"""Transport parameters, catalog coverage, and CLI lifecycle."""

import json

import pytest

from espresso_m1.cli import main as cli_main
from espresso_m1.params import SimParams
from espresso_m1.tparams import (
    TransportParams,
    apply_overrides,
    from_dict,
    inlet_concentration,
    load_catalog,
    to_dict,
)


def test_transport_validation():
    with pytest.raises(ValueError):
        from_dict({"mode": "espresso"})
    with pytest.raises(ValueError):
        from_dict({"M_dose": -1.0})
    with pytest.raises(ValueError):
        from_dict({"f_sol": 1.5})
    with pytest.raises(ValueError):
        from_dict({"D": float("nan")})
    with pytest.raises(TypeError):
        from_dict({"K": True})
    with pytest.raises(ValueError):
        from_dict({"C_in_schedule": [[5.0, 1.0], [2.0, 1.0]]})  # unsorted
    with pytest.raises(ValueError):
        from_dict({"C_in_schedule": [[-1.0, 1.0]]})
    # Extraction demands zero liquid solute so the cup is coffee-attributed.
    with pytest.raises(ValueError):
        from_dict({"mode": "extraction", "C_init": 1.0})
    with pytest.raises(ValueError):
        from_dict({"mode": "extraction", "C_inlet": 1.0})
    # Tracer is strictly passive.
    with pytest.raises(ValueError):
        from_dict({"mode": "tracer", "f_sol": 0.1})
    with pytest.raises(ValueError):
        from_dict({"mode": "tracer", "exchange_rate": 0.1})
    with pytest.raises(KeyError):
        from_dict({"fines": 1.0})


def test_inlet_schedule_piecewise_constant():
    p = from_dict({"mode": "tracer", "f_sol": 0.0, "exchange_rate": 0.0,
                   "C_inlet": 5.0, "C_in_schedule": [[10.0, 0.0], [20.0, 3.0]]})
    assert inlet_concentration(p, 0.0) == 5.0
    assert inlet_concentration(p, 10.0) == 0.0
    assert inlet_concentration(p, 19.999) == 0.0
    assert inlet_concentration(p, 20.0) == 3.0
    assert inlet_concentration(p, 1e6) == 3.0


def test_transport_catalog_entries_all_illustrative():
    catalog = load_catalog()["parameters"]
    for name in to_dict(TransportParams()):
        entry = catalog[name]
        assert entry["provenance"] == "illustrative assumption", name
        assert entry["admissible_range"], name
        assert entry["provenance_detail"], name
    assert catalog["mode"]["default"] == "extraction"


def test_transport_overrides():
    p = apply_overrides(TransportParams(), ["t_end=30", "dt=0.25"])
    assert (p.t_end, p.dt) == (30.0, 0.25)
    with pytest.raises(KeyError):
        apply_overrides(TransportParams(), ["bogus=1"])


def test_transport_cli_success_failure_success(tmp_path):
    out = tmp_path / "t"
    common = ["transport", "--config", "configs/uniform.json",
              "--tconfig", "configs/transport_tracer.json",
              "--tset", "t_end=2", "--out-dir", str(out)]
    assert cli_main(common) == 0
    assert (out / "transport_summary.json").exists()
    assert (out / "transport_timeseries.csv").exists()
    assert (out / "transport_fields.npz").exists()

    assert cli_main(common + ["--tset", "t_end=2", "--set", "max_iter=1"]) == 1
    left = sorted(p.name for p in out.iterdir())
    assert "transport_summary.json" in left  # failure status, not stale data
    failure = json.loads((out / "transport_summary.json").read_text())
    assert failure["converged"] is False and "transport_params" in failure
    assert not (out / "transport_timeseries.csv").exists()
    assert not (out / "transport_fields.npz").exists()

    assert cli_main(common) == 0
    summary = json.loads((out / "transport_summary.json").read_text())
    assert summary["converged"] is True


def test_transport_cli_strict_json_and_zero_cup(tmp_path):
    # Sealed hydraulics: cup stays empty; TDS null in strict JSON.
    out = tmp_path / "sealed"
    rc = cli_main(["transport", "--config", "configs/uniform.json",
                   "--tconfig", "configs/transport_extraction.json",
                   "--set", "R_perf=0", "--tset", "t_end=2",
                   "--out-dir", str(out)])
    assert rc == 0

    def reject_constant(value):
        raise ValueError(value)

    summary = json.loads((out / "transport_summary.json").read_text(),
                         parse_constant=reject_constant)
    assert summary["TDS_pct"] is None
    assert summary["TDS_status"] == "undefined_empty_cup"
    assert summary["EY_pct"] == 0.0


def test_transport_cli_rejects_bad_input(tmp_path):
    out = tmp_path / "bad"
    assert cli_main(["transport", "--config", "configs/uniform.json",
                     "--tconfig", "configs/transport_tracer.json",
                     "--tset", "t_end=-5",
                     "--out-dir", str(out)]) == 2
