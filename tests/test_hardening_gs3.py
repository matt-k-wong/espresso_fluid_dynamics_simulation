"""GS3 prompt 03 part A regressions: invalid-input lifecycle, horizons, boundaries."""

import json

import numpy as np
import pytest

from espresso_m1.cli import main as cli_main
from espresso_m1.driver import run_case
from espresso_m1.params import SimParams
from espresso_m1.tparams import TransportParams, from_dict as t_from_dict
from espresso_m1.transport import (TransportError, assemble_transport,
                                   freeze_snapshot, run_transport)


def _strict_load(path):
    def reject(v):
        raise ValueError(v)
    return json.loads(path.read_text(), parse_constant=reject)


def test_run_success_invalid_success_malformed(tmp_path):
    out = tmp_path / "case"
    assert cli_main(["run", "--config", "configs/uniform.json",
                     "--out-dir", str(out)]) == 0
    (out / "notes.txt").write_text("user file")
    bad = tmp_path / "bad.json"
    bad.write_text("{ not json")
    assert cli_main(["run", "--config", str(bad),
                     "--out-dir", str(out)]) == 2
    assert not (out / "summary.json").exists()
    assert not (out / "fields.npz").exists()
    assert (out / "notes.txt").read_text() == "user file"
    failure = _strict_load(out / "failure.json")
    assert failure["converged"] is False
    assert failure["error_type"] == "JSONDecodeError"
    assert failure["params"] is None
    assert failure["config"] == str(bad)
    assert cli_main(["run", "--config", "configs/uniform.json",
                     "--out-dir", str(out)]) == 0
    assert (out / "summary.json").exists()
    assert not (out / "failure.json").exists()
    assert (out / "notes.txt").read_text() == "user file"


def test_run_success_invalid_success_validation_and_nonconvergence(tmp_path):
    out = tmp_path / "case"
    assert cli_main(["run", "--config", "configs/uniform.json",
                     "--out-dir", str(out)]) == 0
    # Physical validation failure (negative R_perf rejected before solving).
    assert cli_main(["run", "--config", "configs/uniform.json",
                     "--set", "R_perf=-0.01",
                     "--out-dir", str(out)]) == 2
    failure = _strict_load(out / "failure.json")
    assert failure["converged"] is False
    assert failure["error_type"] == "ValueError"
    # Current input context, not old parameters.
    assert failure["set"] == ["R_perf=-0.01"]
    assert failure["params"] is None
    # Solver nonconvergence is a distinct failed run with current params.
    assert cli_main(["run", "--config", "configs/uniform.json",
                     "--set", "max_iter=1",
                     "--out-dir", str(out)]) == 1
    failure = _strict_load(out / "failure.json")
    assert failure["converged"] is False
    assert failure["error_type"] == "ConvergenceError"
    assert failure["params"]["max_iter"] == 1
    assert cli_main(["run", "--config", "configs/uniform.json",
                     "--out-dir", str(out)]) == 0
    assert _strict_load(out / "summary.json")["converged"] is True


def test_run_boolean_json_is_invalid_with_artifact(tmp_path):
    out = tmp_path / "case"
    cfg = tmp_path / "bool.json"
    cfg.write_text(json.dumps({"Nr": True}))
    assert cli_main(["run", "--config", "configs/uniform.json",
                     "--out-dir", str(out)]) == 0
    assert cli_main(["run", "--config", str(cfg),
                     "--out-dir", str(out)]) == 2
    failure = _strict_load(out / "failure.json")
    assert failure["error_type"] == "TypeError"
    assert not (out / "summary.json").exists()


def test_transport_success_invalid_success(tmp_path):
    out = tmp_path / "t"
    base = ["transport", "--config", "configs/uniform.json",
            "--tconfig", "configs/transport_tracer.json",
            "--tset", "t_end=2", "--out-dir", str(out)]
    assert cli_main(base) == 0
    (out / "keep.txt").write_text("keep")
    # Validation failure before solving.
    assert cli_main(base + ["--tset", "t_end=-1"]) == 2
    assert (out / "keep.txt").read_text() == "keep"
    assert not (out / "transport_timeseries.csv").exists()
    assert not (out / "transport_fields.npz").exists()
    failure = _strict_load(out / "transport_summary.json")
    assert failure["converged"] is False
    assert failure["error_type"] == "ValueError"
    assert failure["transport_params"] is None
    assert "t_end=-1" in str(failure["transport_overrides"])
    # Solver nonconvergence.
    assert cli_main(base + ["--set", "max_iter=1"]) == 1
    failure = _strict_load(out / "transport_summary.json")
    assert failure["error_type"] == "ConvergenceError"
    assert failure["hydraulic_params"]["max_iter"] == 1
    assert cli_main(base) == 0
    assert _strict_load(out / "transport_summary.json")["converged"] is True


def test_small_horizon_advances_to_t_end():
    h = SimParams(Nr=4, Nz=8)
    case = run_transport(h, TransportParams(t_end=1e-13, dt=0.1))
    assert case.converged, case.error
    assert case.n_steps == 1
    assert case.result.steps[-1].t == 1e-13


def test_noninteger_step_count_reaches_endpoint():
    h = SimParams(Nr=4, Nz=8)
    case = run_transport(h, TransportParams(t_end=1.0, dt=0.3))
    assert case.converged, case.error
    assert case.n_steps == 4
    assert case.result.steps[-1].t == 1.0
    dts = [s.dt for s in case.result.steps]
    assert dts == pytest.approx([0.3, 0.3, 0.3, 0.1])


def test_schedule_boundaries_hit_exactly_no_pulse_skipped():
    h = SimParams(Nr=4, Nz=8, permeability_mode="uniform")
    tp = t_from_dict({"mode": "tracer", "f_sol": 0.0, "exchange_rate": 0.0,
                      "C_inlet": 5.0,
                      "C_in_schedule": [[0.7, 0.0], [0.9, 3.0]],
                      "t_end": 1.0, "dt": 0.5})
    case = run_transport(h, tp)
    assert case.converged, case.error
    times = [s.t for s in case.result.steps]
    assert 0.7 in times and 0.9 in times and times[-1] == 1.0
    # No step crosses a breakpoint: each step lies within one inlet level.
    for s in case.result.steps:
        assert s.t <= 1.0


def test_mixed_sign_boundary_rejected_despite_positive_aggregate():
    case = run_case(SimParams(Nr=4, Nz=8, permeability_mode="uniform"))
    assert case.converged
    saved = np.array(case.sol.Fz)
    try:
        case.sol.Fz[0, 0] = -1e-7
        case.sol.Fz[1, 0] += 2e-7
        assert float(np.sum(case.sol.Fz[:, 0])) > 0.0
        with pytest.raises(TransportError, match="face"):
            freeze_snapshot(case, 0.018)
    finally:
        case.sol.Fz[:, :] = saved


def test_assemble_rejects_reversed_top_and_preserves_interior():
    from espresso_m1.grid import make_grid
    from espresso_m1.transport import zero_flow_snapshot
    grid = make_grid(2, 1, 0.029, 0.02, 0.0)
    snap = zero_flow_snapshot(grid)
    snap.Fz[0, 0] = -1e-9
    with pytest.raises(TransportError, match="reversed top"):
        assemble_transport(grid, snap.phi, 0.0, snap.Fr, snap.Fz)
    # Valid signed interior flows still assemble (two-cell check).
    snap2 = zero_flow_snapshot(grid)
    snap2.Fr[1, 0] = -2.0e-10
    asm = assemble_transport(grid, snap2.phi, 1e-9, snap2.Fr, snap2.Fz)
    assert asm.L.nnz > 0


def test_sealed_and_forward_freeze_still_succeed():
    sealed = run_case(SimParams(Nr=4, Nz=8, R_perf=0.0))
    assert sealed.converged
    snap = freeze_snapshot(sealed, 0.018)
    assert snap.Q_in == 0.0 and snap.Q_out == 0.0
    fwd = run_case(SimParams(Nr=4, Nz=8, permeability_mode="uniform"))
    snap2 = freeze_snapshot(fwd, 0.018)
    assert snap2.Q_in > 0.0 and snap2.Q_out > 0.0
