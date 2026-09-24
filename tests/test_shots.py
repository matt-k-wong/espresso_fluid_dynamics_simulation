"""GS3 prompt 03 part B: importer/metric correctness on synthetic fixtures only."""

import json

import pytest

from espresso_m1.cli import main as cli_main
from espresso_m1.shots import (brew_ratio, compare_shots, cup_ey_pct,
                               derive_flow, load_shot_file, parse_shot,
                               repeatability, ShotError)


def test_linear_trace_known_derivative():
    shot = load_shot_file("examples/shots/synthetic_linear.json")
    assert shot.metadata.source == "synthetic"
    flows, info = derive_flow(shot.samples, window_s=2.0)
    assert info["window_s"] == 2.0
    assert info["method"] == "local_linear_regression"
    assert info["endpoint_treatment"] == "one_sided_window"
    for q in flows:
        assert q == pytest.approx(1.2, rel=1e-9)
    # Raw retained, no monotonicity forcing.
    assert shot.raw_samples_retained is True
    assert [s.beverage_mass_g for s in shot.samples][:3] == [0.0, 1.2, 2.4]


def test_irregular_trace_and_missing_tds():
    shot = load_shot_file("examples/shots/synthetic_irregular.json")
    assert len(shot.samples) == 16
    flows, _ = derive_flow(shot.samples, window_s=3.0)
    for q in flows[2:-2]:
        assert q == pytest.approx(1.2, rel=0.05)
    assert shot.metadata.tds_pct is None
    assert cup_ey_pct(36.0, None, 18.0) is None
    # Explicit missing TDS fixture never infers EY.
    nod = load_shot_file("examples/shots/synthetic_no_tds.json")
    assert nod.metadata.tds_pct is None
    assert cup_ey_pct(nod.samples[-1].beverage_mass_g, nod.metadata.tds_pct,
                      nod.metadata.dose_g) is None


def test_cup_ey_units_independently_calculated():
    shot = load_shot_file("examples/shots/synthetic_ey.json")
    ey = cup_ey_pct(shot.samples[-1].beverage_mass_g, shot.metadata.tds_pct,
                    shot.metadata.dose_g)
    # 36.0 g * 8.5 % / 18.0 g = 17.0 % by hand.
    assert ey == pytest.approx(36.0 * 8.5 / 18.0, rel=1e-12)
    assert ey == pytest.approx(17.0, rel=1e-12)
    assert brew_ratio(36.0, 18.0) == pytest.approx(2.0)


def test_csv_series_import():
    shot = load_shot_file("examples/shots/synthetic_linear_meta.json",
                          "examples/shots/synthetic_linear_series.csv")
    assert len(shot.samples) == 11
    assert shot.samples[-1].beverage_mass_g == pytest.approx(12.0)


def test_validation_rejections():
    base = json.loads(open("examples/shots/synthetic_linear.json").read())
    # Negative dose.
    bad = dict(base, dose_g=-1.0)
    with pytest.raises(ShotError):
        parse_shot(bad)
    # Nonfinite mass.
    bad2 = json.loads(json.dumps(base))
    bad2["samples"][3]["beverage_mass_g"] = float("inf")
    # json round-trip turns inf into Infinity (invalid JSON); inject directly.
    bad2["samples"][3]["beverage_mass_g"] = float("inf")
    with pytest.raises(ShotError):
        parse_shot(bad2)
    # Decreasing times.
    bad3 = json.loads(json.dumps(base))
    bad3["samples"][5]["time_s"] = 0.5
    with pytest.raises(ShotError, match="decreasing"):
        parse_shot(bad3)
    # Duplicates without rule.
    bad4 = json.loads(json.dumps(base))
    bad4["samples"][2] = dict(bad4["samples"][1])
    bad4["samples"][2]["time_s"] = bad4["samples"][1]["time_s"]
    with pytest.raises(ShotError, match="duplicate"):
        parse_shot(bad4)
    # Duplicates with rule aggregate.
    shot = parse_shot(bad4, aggregate="mean")
    assert len(shot.samples) == len(base["samples"]) - 1
    # Inconsistent units: dose in kg slip (18000 g) or pressure in bar.
    bad5 = json.loads(json.dumps(base))
    bad5["dose_g"] = 18000.0
    with pytest.raises(ShotError, match="inconsistent|plausible"):
        parse_shot(bad5)


def test_negative_increments_flagged_not_removed():
    shot = load_shot_file("examples/shots/synthetic_linear.json")
    # Inject one noisy dip in memory (raw retained by API contract).
    shot.samples[5].beverage_mass_g -= 1.5
    flows, info = derive_flow(shot.samples, window_s=2.0)
    assert info["n_negative_raw_increments"] == 1
    # Raw still dips (not rewritten to monotonic).
    assert shot.samples[5].beverage_mass_g < shot.samples[4].beverage_mass_g
    assert "not directly measured" in info["note"]


def test_compare_holds_out_whole_shots(tmp_path):
    a = load_shot_file("examples/shots/synthetic_linear.json")
    b = load_shot_file("examples/shots/synthetic_irregular.json")
    # Train on a, hold out whole b (never split rows of one shot).
    result = compare_shots([a], [b], t0_s=0.0)
    assert result["fit"]["q_g_s"] == pytest.approx(1.2, rel=1e-9)
    assert result["heldout_residuals"][0]["shot_id"] == b.metadata.shot_id
    assert result["heldout_residuals"][0]["rmse_g"] < 1e-9
    assert "DISABLED" in result["calibration_status"]
    assert "constant Q" in result["limitations"] or "constant-flow" in result["limitations"]


def test_shot_cli_validate_metrics_plot_compare(tmp_path):
    assert cli_main(["shot-validate", "--shot",
                     "examples/shots/synthetic_linear.json"]) == 0
    assert cli_main(["shot-validate", "--shot",
                     "examples/shots/synthetic_linear.json",
                     "--series", "examples/shots/synthetic_linear_series.csv"]) != 0
    out_json = tmp_path / "metrics.json"
    assert cli_main(["shot-metrics", "--shot",
                     "examples/shots/synthetic_ey.json",
                     "--out", str(out_json)]) == 0
    summary = json.loads(out_json.read_text())
    assert summary["cup_EY_pct"] == pytest.approx(17.0, rel=1e-12)
    assert summary["cup_EY_status"] == "measured"
    nod_json = tmp_path / "nod.json"
    assert cli_main(["shot-metrics", "--shot",
                     "examples/shots/synthetic_no_tds.json",
                     "--out", str(nod_json)]) == 0
    assert json.loads(nod_json.read_text())["cup_EY_pct"] is None
    mass_svg = tmp_path / "mass.svg"
    flow_svg = tmp_path / "flow.svg"
    assert cli_main(["shot-plot", "--shot",
                     "examples/shots/synthetic_linear.json",
                     "--out", str(mass_svg), "--flow-out", str(flow_svg),
                     "--baseline-q", "1.2", "--baseline-t0", "0.0",
                     "--baseline-m0", "0.0"]) == 0
    mass_text = mass_svg.read_text()
    assert "time_s" in mass_text or "mass" in mass_text
    assert "observed" in mass_text and "simulated" in mass_text
    cmp_dir = tmp_path / "cmp"
    assert cli_main(["shot-compare", "--shots",
                     "examples/shots/synthetic_linear.json",
                     "examples/shots/synthetic_irregular.json",
                     "--already-wet-t0", "0.0",
                     "--out-dir", str(cmp_dir)]) == 0
    comp = json.loads((cmp_dir / "comparison.json").read_text())
    assert comp["fit"]["q_g_s"] == pytest.approx(1.2, rel=1e-6)
    assert (cmp_dir / "mass.svg").exists()
