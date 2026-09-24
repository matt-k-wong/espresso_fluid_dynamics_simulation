"""Acceptance checks 4-5: linearity without compaction; compaction behavior."""

import numpy as np
import pytest

from espresso_m1.driver import run_case, run_sweep
from espresso_m1.params import BAR_PA

from .conftest import compaction_params, uniform_params


def test_linearity_without_compaction():
    """With compaction disabled, Q must increase linearly with pressure."""
    qs = []
    for bar in (3.0, 6.0, 9.0):
        case = run_case(uniform_params(p_in=bar * BAR_PA))
        assert case.converged
        qs.append(case.diag.Q_out_m3s)
    assert qs[1] / qs[0] == pytest.approx(2.0, rel=1e-9)
    assert qs[2] / qs[0] == pytest.approx(3.0, rel=1e-9)


def test_compaction_bounds_and_positivity():
    case = run_case(compaction_params())
    assert case.converged, case.error
    sol = case.sol
    assert sol is not None
    assert np.all(sol.phi >= 0.20 - 1e-12)
    assert np.all(sol.phi <= 0.35 + 1e-12)
    assert np.all(sol.k > 0.0)
    assert np.all(np.isfinite(sol.p))


def test_compaction_only_reduces_flow():
    q_uni = run_case(uniform_params()).diag.Q_out_m3s
    q_cmp = run_case(compaction_params()).diag.Q_out_m3s
    assert 0.0 < q_cmp < q_uni


def test_uniform_converges_with_finite_cost():
    # Uniform mode is linear: one solve plus one confirmation solve. The
    # repaired solver returns exactly the checked state (no extra polish
    # solve), so this encodes finite cost, not a fixed solve count.
    case = run_case(uniform_params())
    assert case.converged
    assert case.iters == 2


def test_nonconvergence_is_reported_not_silent():
    from espresso_m1.solver import ConvergenceError, solve_pressure
    from espresso_m1.grid import make_grid

    params = compaction_params(max_iter=1)  # cannot confirm an update in 1 pass
    grid = make_grid(params.Nr, params.Nz, params.R, params.L, params.R_perf)
    with pytest.raises(ConvergenceError):
        solve_pressure(params, grid)
    # The driver must surface the failure instead of returning a solution.
    case = run_case(params)
    assert not case.converged
    assert case.sol is None and case.diag is None
    assert case.error != ""


def test_sweep_records_all_points():
    rows = run_sweep(uniform_params(), [3.0 * BAR_PA, 9.0 * BAR_PA])
    assert len(rows) == 2
    assert all(r.converged for r in rows)
    assert rows[1].Q_m3s / rows[0].Q_m3s == pytest.approx(3.0, rel=1e-9)
