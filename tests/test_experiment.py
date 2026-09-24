"""M6: save/reproduce roundtrip on synthetic fixtures only."""

import json
from pathlib import Path

from espresso_m1.cli import main as cli_main


LIN = "examples/shots/synthetic_linear.json"
IRR = "examples/shots/synthetic_irregular.json"


def _strict(path: Path) -> dict:
    def reject(v):
        raise ValueError(v)
    return json.loads(path.read_text(), parse_constant=reject)


def test_save_reproduce_roundtrip(tmp_path):
    bundle = tmp_path / "exp"
    rc = cli_main(["experiment-save", "--shots", LIN, IRR,
                   "--already-wet-t0", "0.0", "--label", "demo",
                   "--out-dir", str(bundle)])
    assert rc == 0
    manifest = _strict(bundle / "manifest.json")
    assert manifest["experiment_schema_version"] == "experiment_v1"
    assert manifest["train_ids"] == ["synthetic-linear-01"]
    assert manifest["holdout_ids"] == ["synthetic-irregular-01"]
    prov = {e["shot_id"]: e["provenance"] for e in manifest["inputs"]}
    assert prov == {"synthetic-linear-01": "synthetic",
                    "synthetic-irregular-01": "synthetic"}
    assert "fitted_empirical" in manifest["inputs_classification"]
    assert "DISABLED" in manifest["calibration_status"]
    assert (bundle / "comparison.json").exists()
    assert (bundle / "mass.svg").exists() and (bundle / "flow.svg").exists()
    assert (bundle / "export" / "combined_series.csv").exists()
    assert (bundle / "export" / "metrics.csv").exists()
    assert (bundle / "README.txt").exists()
    mass_svg = (bundle / "mass.svg").read_text()
    assert "time_s (s)" in mass_svg and "beverage_mass_g (g)" in mass_svg
    assert "observed" in mass_svg and "simulated" in mass_svg
    rc2 = cli_main(["experiment-reproduce", "--dir", str(bundle)])
    assert rc2 == 0
    report = _strict(bundle / "reproduce_report.json")
    assert report["passed"] is True
    assert report["checks"]["hash_ok"] is True
    assert report["checks"]["q_match"] is True


def test_reproduce_detects_tamper(tmp_path):
    bundle = tmp_path / "exp"
    assert cli_main(["experiment-save", "--shots", LIN, IRR,
                     "--already-wet-t0", "0.0",
                     "--out-dir", str(bundle)]) == 0
    target = bundle / "inputs" / "shots" / Path(IRR).name
    text = target.read_text().replace("36.0", "99.0")
    target.write_text(text)
    rc = cli_main(["experiment-reproduce", "--dir", str(bundle)])
    assert rc == 1
    report = _strict(bundle / "reproduce_report.json")
    assert report["passed"] is False
    assert report["checks"]["hash_ok"] is False


def test_save_rejects_duplicate_and_unknown_holdout(tmp_path):
    rc = cli_main(["experiment-save", "--shots", LIN, LIN,
                   "--already-wet-t0", "0.0",
                   "--out-dir", str(tmp_path / "a")])
    assert rc == 2
    rc = cli_main(["experiment-save", "--shots", LIN, IRR,
                   "--already-wet-t0", "0.0", "--holdout", "nope",
                   "--out-dir", str(tmp_path / "b")])
    assert rc == 2
