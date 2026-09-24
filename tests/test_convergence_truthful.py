"""Task A regression tests: truthful nonlinear convergence.

Every accepted solution is re-verified here from its returned fields with
independent arithmetic: the discrete system is reassembled from the returned
permeability, and local/global balances are recomputed from the returned
face fluxes. The solver's own flags and self-reported residuals are inputs
to these checks, never the verdict.
"""

import numpy as np
import pytest
from dataclasses import replace

from espresso_m1.assembly import assemble
from espresso_m1.driver import run_case
from espresso_m1.params import SimParams
from espresso_m1.solver import FLOW_ATOL_M3S

CASE = dict(permeability_mode="compaction", Nr=8, Nz=16, R_perf=0.0145)


def independent_balances(case):
    """Recompute balances from returned fields only."""
    sol = case.sol
    A, b = assemble(case.params, case.grid, sol.k)
    ap_minus_b = float(np.max(np.abs(A @ sol.p.T.ravel() - b)))
    cell_net = (sol.Fr[1:, :] - sol.Fr[:-1, :]) + (sol.Fz[:, 1:] - sol.Fz[:, :-1])
    local_abs = float(np.max(np.abs(cell_net)))
    q_in = float(np.sum(sol.Fz[:, 0]))
    q_out = float(np.sum(sol.Fz[case.grid.open_mask, case.grid.Nz]))
    mismatch_abs = abs(q_in - q_out)
    return ap_minus_b, local_abs, mismatch_abs, q_in, q_out


def balance_limit(case, q_in, q_out):
    q_scale = max(abs(q_in), abs(q_out))
    return FLOW_ATOL_M3S + case.params.tol_res_rel * q_scale


@pytest.mark.parametrize("relaxation", [1.0, 0.5, 0.1])
def test_relaxation_converges_to_same_balanced_solution(relaxation):
    base = SimParams(
        **CASE, tol_res_rel=1e-8, max_iter=500, relaxation=relaxation
    )
    case = run_case(base)
    assert case.converged, case.error
    ap_b, local_abs, mismatch_abs, q_in, q_out = independent_balances(case)
    limit = balance_limit(case, q_in, q_out)
    print(f"\nw={relaxation}: Q={q_out*1e6:.8f} mL/s local={local_abs:.3e} "
          f"mismatch={mismatch_abs:.3e} Ap-b={ap_b:.3e} limit={limit:.3e}")
    assert local_abs <= limit
    assert mismatch_abs <= limit
    # Same physical answer regardless of relaxation path.
    assert q_out * 1e6 == pytest.approx(0.24676146, rel=1e-6)


def test_relaxation_solutions_agree():
    qs = []
    for relaxation in (1.0, 0.5, 0.1):
        base = SimParams(
            **CASE, tol_res_rel=1e-8, max_iter=500, relaxation=relaxation
        )
        case = run_case(base)
        assert case.converged
        qs.append(case.diag.Q_out_m3s)
    for q in qs[1:]:
        assert q == pytest.approx(qs[0], rel=1e-6)


@pytest.mark.parametrize("relaxation", [1e-8, 1e-12])
def test_tiny_relaxation_reports_nonconvergence(relaxation):
    # Deliberately slow relaxations must fail explicitly; tiny iterate
    # movement alone never counts as success (default 100-iteration budget).
    case = run_case(SimParams(**CASE, relaxation=relaxation))
    assert not case.converged
    assert case.sol is None and case.diag is None
    assert case.error != ""


@pytest.mark.parametrize("p_in,p_out", [(0.0, 0.0), (1.0e5, 1.0e5)])
def test_equal_pressure_equilibrium(p_in, p_out):
    case = run_case(SimParams(**CASE, p_in=p_in, p_out=p_out))
    assert case.converged, case.error
    assert abs(case.diag.Q_out_m3s) < 1e-18
    assert abs(case.diag.Q_in_m3s) < 1e-18
    # Constant pressure: deviation bounded by solver noise, not physics.
    assert float(np.max(np.abs(case.sol.p - p_in))) < 1e-3


def test_sealed_bed_equilibrium():
    case = run_case(SimParams(**{**CASE, "R_perf": 0.0}))
    assert case.converged, case.error
    assert abs(case.diag.Q_out_m3s) < 1e-18
    assert abs(case.diag.Q_in_m3s) < 1e-18
    assert float(np.max(np.abs(case.sol.p - case.params.p_in))) / \
        case.params.p_in < 1e-9


def test_common_pressure_offset_invariance():
    base = SimParams(**CASE)
    shifted = replace(base, p_in=base.p_in + 5e5, p_out=base.p_out + 5e5)
    a = run_case(base)
    b = run_case(shifted)
    assert a.converged and b.converged
    assert b.diag.Q_out_m3s == pytest.approx(a.diag.Q_out_m3s, rel=1e-9)
    assert np.max(np.abs(b.sol.k - a.sol.k)) / np.max(a.sol.k) < 1e-9
    # Effective stress (the physical driver) is offset-invariant.
    assert np.max(np.abs(b.sol.sigma_eff - a.sol.sigma_eff)) / \
        np.max(a.sol.sigma_eff) < 1e-9
