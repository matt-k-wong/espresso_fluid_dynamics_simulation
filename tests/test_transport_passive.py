"""Task E1/E2: passive tracer conservation and two-cell flux directions.

E1: zero initial/inlet solute with B=0 stays zero; a uniform initial and
inlet C is preserved by divergence-free flow up to hydraulic residual
tolerance. E2: positive and negative interior F plus a D>0 gradient move
exactly the mass the donor loses to the receiver; assembly entries are
checked against hand-computed F+H values on two different grids so a
double multiplication by area cannot hide in a shared implementation.
"""

import numpy as np
import pytest

from espresso_m1.grid import make_grid
from espresso_m1.params import SimParams
from espresso_m1.tparams import TransportParams
from espresso_m1.transport import (
    assemble_transport,
    run_transport,
    transport_step,
    zero_flow_snapshot,
)


def test_zero_tracer_stays_zero():
    h = SimParams(permeability_mode="compaction", Nr=8, Nz=16, R_perf=0.0145)
    t = TransportParams(mode="tracer", f_sol=0.0, exchange_rate=0.0,
                        C_init=0.0, C_inlet=0.0, D=1e-9, t_end=5.0, dt=0.5)
    case = run_transport(h, t)
    assert case.converged
    assert np.max(np.abs(case.result.C)) == 0.0
    assert case.result.M_cup == 0.0
    assert case.result.budget_max_abs == 0.0


def test_uniform_tracer_preserved_to_hydraulic_tolerance():
    h = SimParams(permeability_mode="compaction", Nr=8, Nz=16, R_perf=0.0145)
    t = TransportParams(mode="tracer", f_sol=0.0, exchange_rate=0.0,
                        C_init=2.0, C_inlet=2.0, D=1e-9, t_end=5.0, dt=0.5)
    case = run_transport(h, t)
    assert case.converged
    # Divergence-free only to the hydraulic balance tolerance.
    assert np.max(np.abs(case.result.C - 2.0)) / 2.0 < 1e-7


def two_cell_setup(R, F, D=1e-9, c_left=3.0):
    grid = make_grid(2, 1, R, 0.02, 0.0)  # sealed: no boundary fluxes
    snap = zero_flow_snapshot(grid, phi_value=0.35)
    snap.Fr[1, 0] = F
    asm = assemble_transport(grid, snap.phi, D, snap.Fr, snap.Fz)
    C = np.array([[c_left], [0.0]])
    return grid, snap, asm, C


def hand_J(R, F, D, c_left):
    dr = R / 2.0
    area = 2.0 * np.pi * (R / 2.0) * 0.02
    H = D * 0.35 * area / dr
    if F >= 0.0:
        return F * c_left + H * c_left, H
    return H * c_left, H  # upwind takes the right (zero) cell


@pytest.mark.parametrize("R", [0.029, 0.058])
@pytest.mark.parametrize("F", [2.0e-10, -2.0e-10])
def test_two_cell_exact_transfer_both_signs(R, F):
    grid, snap, asm, C = two_cell_setup(R, F)
    J, H = hand_J(R, F, 1e-9, 3.0)
    net = asm.L @ C.T.ravel()
    assert net[0] == pytest.approx(J, rel=1e-12)
    assert net[1] == pytest.approx(-J, rel=1e-12)
    # Advective part is F*C with no extra area factor; dispersive part is H.
    L = asm.L.toarray()
    if F >= 0.0:
        assert L[0, 0] - H == pytest.approx(F, rel=1e-12)
        assert L[0, 1] == pytest.approx(-H, rel=1e-12)
        assert L[1, 0] == pytest.approx(-(F + H), rel=1e-12)
        assert L[1, 1] == pytest.approx(H, rel=1e-12)
    else:
        assert L[0, 0] == pytest.approx(H, rel=1e-12)
        assert L[0, 1] == pytest.approx(F - H, rel=1e-12)
        assert L[1, 0] == pytest.approx(-H, rel=1e-12)
        assert L[1, 1] == pytest.approx(H - F, rel=1e-12)
    # One step: donor loss equals receiver gain to linear-solver roundoff.
    B = np.zeros_like(C)
    C1, _, _ = transport_step(C, B, asm, snap.W, snap.m, 0.0, 0.0, 1.0, 0.0)
    d0 = snap.W[0, 0] * (C1[0, 0] - C[0, 0])
    d1 = snap.W[1, 0] * (C1[1, 0] - C[1, 0])
    assert d0 == pytest.approx(-1.0 * J, rel=1e-9)
    assert abs(d0 + d1) <= 1e-9 * abs(J) * 1.0
