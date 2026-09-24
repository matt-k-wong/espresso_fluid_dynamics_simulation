"""Acceptance check 3: boundary behavior (axis, sidewall, solid rim)."""

import numpy as np

from espresso_m1.driver import run_case

from .conftest import compaction_params, uniform_params


def test_zero_flux_boundaries():
    params = uniform_params(Nr=20, Nz=40, R_perf=0.0145)
    case = run_case(params)
    assert case.converged, case.error
    sol = case.sol
    assert sol is not None
    # Axis (r = 0) and sidewall (r = R): identically zero radial flux.
    assert np.all(sol.Fr[0, :] == 0.0)
    assert np.all(sol.Fr[params.Nr, :] == 0.0)
    # Solid bottom rim: identically zero axial flux.
    closed = ~case.grid.open_mask
    assert np.any(closed), "test needs a closed rim"
    assert np.all(sol.Fz[closed, params.Nz] == 0.0)
    # Open perforations carry positive downward outflow; top flows inward.
    assert np.all(sol.Fz[case.grid.open_mask, params.Nz] > 0.0)
    assert np.all(sol.Fz[:, 0] > 0.0)


def test_zero_flux_boundaries_compaction():
    params = compaction_params(Nr=12, Nz=24, R_perf=0.0145)
    case = run_case(params)
    assert case.converged, case.error
    assert case.sol is not None
    assert np.all(case.sol.Fr[0, :] == 0.0)
    assert np.all(case.sol.Fr[params.Nr, :] == 0.0)
    assert np.all(case.sol.Fz[~case.grid.open_mask, params.Nz] == 0.0)


def test_partial_perforation_reduces_flow():
    full = run_case(uniform_params(R_perf=0.029))
    part = run_case(uniform_params(R_perf=0.0145))
    assert full.converged and part.converged
    assert part.diag.Q_out_m3s < full.diag.Q_out_m3s
    # Discrete open area approximates pi*R_perf^2 to O(dr).
    assert abs(part.grid.open_area - part.grid.open_area_ideal) / \
        part.grid.open_area_ideal < part.grid.dr / part.grid.R_perf + 1e-12
