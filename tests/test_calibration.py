"""M4 harness on synthetic data only: fit, uncertainty, verdicts, sim bias."""

import json

import pytest

from espresso_m1.calibration import calibrate
from espresso_m1.cli import main as cli_main
from espresso_m1.shots import parse_shot


def make_shot(shot_id, slope, t_end=30.0, step=1.0):
    n = int(t_end / step)
    return parse_shot({
        "schema_version": "gs3_shot_v1",
        "shot_id": shot_id,
        "source": "synthetic",
        "machine": {"profile_id": "x", "feed": "AV reservoir"},
        "dose_g": 18.0,
        "time_zero": "already_wet",
        "cup": {"TDS_pct": None},
        "samples": [{"time_s": i * step, "beverage_mass_g": slope * i * step}
                    for i in range(n + 1)],
    })


def test_pass_and_sim_bias():
    train = [make_shot("a", 1.20), make_shot("b", 1.24)]
    held = [make_shot("c", 1.22)]
    r = calibrate(train, held, 0.0, sim_q_g_s=1.20)
    assert r["status"] == "pass" and r["passed"] is True
    assert r["fit_q_g_s"] == pytest.approx(1.22, rel=1e-12)
    assert r["ci95_mean_g_s"][0] < 1.22 < r["ci95_mean_g_s"][1]
    assert r["verdicts"][0]["verdict"] == "pass"
    assert r["simulator_comparison"]["bias_empirical_minus_sim_g_s"] == pytest.approx(0.02)
    assert r["simulator_comparison"]["sim_inside_95ci"] is True
    assert r["all_synthetic"] is True


def test_failure_disclosed():
    train = [make_shot("a", 1.20), make_shot("b", 1.24)]
    held = [make_shot("z", 2.50)]
    r = calibrate(train, held, 0.0)
    assert r["status"] == "FAIL" and r["passed"] is False
    assert r["verdicts"][0]["verdict"] == "FAIL"


def test_single_train_indeterminate():
    r = calibrate([make_shot("a", 1.20)], [make_shot("c", 1.20)], 0.0)
    assert r["status"] == "indeterminate_single_train"
    assert r["verdicts"][0]["verdict"] == "unassessed_single_train"
    assert r["ci95_mean_g_s"] is None


def test_cli_calibrate_synthetic(tmp_path):
    out = tmp_path / "cal.json"
    rc = cli_main(["shot-calibrate",
                   "--shots", "examples/shots/synthetic_linear.json",
                   "examples/shots/synthetic_irregular.json",
                   "--already-wet-t0", "0.0", "--sim-q", "1.2",
                   "--out", str(out)])
    assert rc == 0  # single train shot -> indeterminate, harness ran fine
    r = json.loads(out.read_text())
    assert r["status"] == "indeterminate_single_train"
    assert r["simulator_comparison"]["bias_empirical_minus_sim_g_s"] == pytest.approx(0.0)
    out2 = tmp_path / "cal2.json"
    rc = cli_main(["shot-calibrate",
                   "--shots", "examples/shots/synthetic_linear.json",
                   "examples/shots/synthetic_irregular.json",
                   "examples/shots/synthetic_no_tds.json",
                   "--holdout", "synthetic-no-tds-01",
                   "--already-wet-t0", "0.0", "--out", str(out2)])
    # Identical train slopes -> degenerate point interval; 1.2067 held-out
    # slope honestly FAILs and exit code discloses it.
    assert rc == 1
    r2 = json.loads(out2.read_text())
    assert r2["n_train"] == 2 and r2["status"] == "FAIL"
