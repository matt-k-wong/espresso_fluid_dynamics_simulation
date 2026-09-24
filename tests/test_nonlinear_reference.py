"""Task D.1 regression test: nonlinear fully-open flow vs independent integral.

The reference is computed by direct quadrature of the stated constitutive
law, written out here without calling the production closure. For the
submitted defaults at 9 bar the continuum outlet flow is
0.4687763667883215 mL/s (independently established). Mesh errors of the
normalized compaction/uniform flow ratio must shrink at approximately
second order (consecutive error ratios > 3.5) on 10x20, 20x40, 40x80.
A solver tolerance is not a discretization error: this compares the
discrete scheme against the continuum answer.
"""

import numpy as np
import pytest
from dataclasses import replace
from scipy.integrate import quad

from espresso_m1.driver import run_case
from espresso_m1.params import SimParams

Q_EXACT_MLS = 0.4687763667883215


def independent_ratio(params: SimParams) -> float:
    """Integral of k(s)/k0 over the stress interval, by quadrature.

    sigma(t) = sigma_top + t * (p_in - p_out), t in [0, 1]; phi and the
    Carman-Kozeny factor are written out directly from the documented law.
    """
    f0 = params.phi0**3 / (1.0 - params.phi0) ** 2

    def normalized_k(t):
        sigma = params.sigma_eff_top + t * (params.p_in - params.p_out)
        phi = params.phi_min + (params.phi0 - params.phi_min) * np.exp(
            -sigma / params.sigma_c
        )
        return (phi**3 / (1.0 - phi) ** 2) / f0

    return quad(normalized_k, 0.0, 1.0, epsabs=1e-12, epsrel=1e-12)[0]


def test_quadrature_reproduces_independent_continuum_value():
    params = SimParams(permeability_mode="compaction")
    q_exact = (
        np.pi * params.R**2 * params.k0 * (params.p_in - params.p_out)
        / (params.mu * params.L)
        * independent_ratio(params)
        * 1e6
    )
    assert q_exact == pytest.approx(Q_EXACT_MLS, rel=1e-9)


def test_fully_open_mesh_errors_shrink_second_order():
    params = SimParams(permeability_mode="compaction")
    ratio = independent_ratio(params)
    errors = []
    for n in (10, 20, 40):
        sub = replace(params, Nr=n, Nz=2 * n, R_perf=params.R)
        comp = run_case(sub)
        uni = run_case(replace(sub, permeability_mode="uniform"))
        assert comp.converged and uni.converged
        err = abs(comp.diag.Q_out_m3s / uni.diag.Q_out_m3s / ratio - 1.0)
        errors.append(err)
        print(f"\nn={n}: normalized ratio error={err:.4e}")
    assert errors[-1] < 3e-4
    for coarse, fine in zip(errors, errors[1:]):
        assert coarse / fine > 3.5


def test_perforated_ratio_approaches_transform_reference():
    # Kirchhoff-transform consequence: the normalized compaction/uniform
    # flow ratio tends to the same constitutive integral; absolute Q
    # converges more slowly because of the junction regularity.
    params = SimParams(permeability_mode="compaction")
    ratio = independent_ratio(params)
    errors = []
    for n in (10, 20, 40):
        sub = replace(params, Nr=n, Nz=2 * n, R_perf=params.R / 2)
        comp = run_case(sub)
        uni = run_case(replace(sub, permeability_mode="uniform"))
        assert comp.converged and uni.converged
        errors.append(abs(comp.diag.Q_out_m3s / uni.diag.Q_out_m3s / ratio - 1.0))
    assert errors[-1] < 2e-3
    assert all(a > b for a, b in zip(errors, errors[1:]))
