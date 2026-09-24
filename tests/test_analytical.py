"""Acceptance check 1: uniform-k, fully open bottom vs analytical Darcy result."""

import numpy as np
import pytest

from espresso_m1.driver import analytical_cylindrical_Q, run_case

from .conftest import uniform_params


@pytest.mark.parametrize("Nr,Nz", [(10, 20), (20, 40), (40, 80)])
def test_uniform_matches_analytical(Nr, Nz):
    params = uniform_params(Nr=Nr, Nz=Nz)
    case = run_case(params)
    assert case.converged, case.error
    assert case.diag is not None
    expected = analytical_cylindrical_Q(
        params.R, params.k0, params.p_in - params.p_out, params.mu, params.L
    )
    rel = abs(case.diag.Q_out_m3s - expected) / expected
    print(f"\nNr={Nr} Nz={Nz}: Q={case.diag.Q_out_m3s:.6e} "
          f"expected={expected:.6e} rel_err={rel:.3e}")
    assert rel < 1.0e-6


def test_refinement_does_not_degrade():
    errs = []
    for Nr, Nz in [(10, 20), (20, 40), (40, 80)]:
        params = uniform_params(Nr=Nr, Nz=Nz)
        case = run_case(params)
        assert case.converged
        expected = analytical_cylindrical_Q(
            params.R, params.k0, params.p_in - params.p_out, params.mu, params.L
        )
        errs.append(abs(case.diag.Q_out_m3s - expected) / expected)
    # The 1D profile is represented exactly; error must stay at solver noise.
    assert max(errs) < 1.0e-6


def test_uniform_flow_is_one_dimensional():
    params = uniform_params(Nr=20, Nz=40)
    case = run_case(params)
    assert case.converged
    assert case.sol is not None
    # No radial driving force: radial fluxes must vanish to solver precision.
    axial_scale = np.max(np.abs(case.sol.Fz))
    assert np.max(np.abs(case.sol.Fr)) / axial_scale < 1.0e-9
