"""Task E4: passive residence time under advection only.

Uniform fully open hydraulics, D=0, lambda=0. Pore transit time is
tau = sum(W)/Q. A finite inlet pulse with steps aligned to its boundaries
must satisfy input = outlet + stored tracer exactly, and the outlet
centroid must approach tau + Tp/2 under grid and dt refinement. Upwind
spreading of the sharp front is expected; no second-order claim is made.
"""

import numpy as np
import pytest

from espresso_m1.params import SimParams
from espresso_m1.tparams import TransportParams
from espresso_m1.transport import run_transport

T_PULSE = 2.0
C_PULSE = 5.0


def pulse_case(Nr, dt, cfl=0.5, t_end=30.0):
    h = SimParams(permeability_mode="uniform", Nr=Nr, Nz=2 * Nr)
    t = TransportParams(
        mode="tracer", f_sol=0.0, exchange_rate=0.0, C_init=0.0,
        C_inlet=C_PULSE, C_in_schedule=[[T_PULSE, 0.0]], D=0.0,
        cfl_advective=cfl, t_end=t_end, dt=dt,
    )
    return run_transport(h, t)


def centroid(case):
    steps = case.result.steps
    ts = np.array([s.t for s in steps])
    dts = np.array([s.dt for s in steps])
    oc = np.array([s.outlet_C for s in steps])
    rate = case.Q_out * oc
    return float(np.sum((ts - dts / 2) * rate * dts) / np.sum(rate * dts))


def test_pulse_mass_balance_exact():
    case = pulse_case(20, 0.5)
    assert case.converged
    r = case.result
    injected = C_PULSE * case.Q_out * T_PULSE
    assert r.M_cup + r.steps[-1].sum_WC == pytest.approx(injected, rel=1e-9)
    assert r.budget_max_abs <= r.budget_tol


def test_centroid_approaches_transit_time_dt_refinement():
    errs = []
    for dt in (0.5, 0.25, 0.125):
        case = pulse_case(10, dt, cfl=0.0)
        assert case.converged
        err = abs(centroid(case) - (case.tau_s + T_PULSE / 2))
        errs.append(err)
        print(f"\ndt={dt}: centroid err={err:.3e} (tau={case.tau_s:.4f})")
    assert errs[-1] < 1e-6
    assert all(a > b for a, b in zip(errs, errs[1:]))


def test_centroid_approaches_transit_time_grid_refinement():
    errs = []
    for Nr in (10, 20, 40):
        case = pulse_case(Nr, 0.5)
        assert case.converged
        err = abs(centroid(case) - (case.tau_s + T_PULSE / 2))
        errs.append(err)
        print(f"\n{Nr}x{2*Nr}: centroid err={err:.3e} (tau={case.tau_s:.4f})")
    assert errs[-1] / case.tau_s < 0.01
    assert all(a > b for a, b in zip(errs, errs[1:]))
