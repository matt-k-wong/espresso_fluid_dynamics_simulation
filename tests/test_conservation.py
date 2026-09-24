"""Acceptance check 2: conservation (global mismatch + local residuals)."""

import pytest

from espresso_m1.driver import run_case

from .conftest import compaction_params, uniform_params


@pytest.mark.parametrize("mode", ["uniform", "compaction"])
@pytest.mark.parametrize("R_perf", [0.029, 0.0145])
def test_conservation(mode, R_perf):
    maker = uniform_params if mode == "uniform" else compaction_params
    case = run_case(maker(R_perf=R_perf))
    assert case.converged, case.error
    d = case.diag
    assert d is not None
    print(f"\nmode={mode} R_perf={R_perf}: mismatch_rel={d.mismatch_rel:.3e} "
          f"max_res_rel={d.max_residual_rel:.3e}")
    assert d.mismatch_rel < 1.0e-6
    assert d.max_residual_rel < 1.0e-6
