"""Task E5: independent dispersion check on a no-flux eigenmode.

Sealed synthetic flow (F = 0), constant phi, lambda = 0. The cylindrical
diffusion eigenmode

    C = Cbar + A*J0(a*r/R)*cos(pi*z/L)*exp(-D*((a/R)^2+(pi/L)^2)*t),

with a the first positive zero of J1 (~3.83170597), is regular at the axis
and has zero normal derivative at every boundary (J1(a) = 0 kills the
radial derivative at r = R; sin kills the axial derivatives at z = 0, L).
Both radial and axial diffusion contribute to the decay rate. The mode and
its decay are written out here from calculus (scipy supplies only J0/J1
evaluation); grids 8x16, 16x32, 32x64 are compared with time error
controlled by a small dt.
"""

import numpy as np
from scipy.special import j0

from espresso_m1.grid import make_grid
from espresso_m1.transport import (
    assemble_transport,
    transport_step,
    zero_flow_snapshot,
)

A_J1 = 3.83170597
D = 1e-9
CBAR = 2.0
AMP = 1.0
T_END = 3600.0
DT = 5.0


def eigenmode_error(n: int) -> float:
    grid = make_grid(n, 2 * n, 0.029, 0.02, 0.0)
    snap = zero_flow_snapshot(grid, phi_value=0.35)
    asm = assemble_transport(grid, snap.phi, D, snap.Fr, snap.Fz)
    rr = grid.r_centers[:, None] / grid.R
    zz = grid.z_centers[None, :]
    radial = j0(A_J1 * rr)
    axial = np.cos(np.pi * zz / grid.L)
    C = CBAR + AMP * radial * axial
    B = np.zeros_like(C)
    gamma = D * ((A_J1 / grid.R) ** 2 + (np.pi / grid.L) ** 2)
    t = 0.0
    while t < T_END - 1e-9:
        C, B, _ = transport_step(C, B, asm, snap.W, snap.m, 0.0, 0.0, DT, 0.0)
        t += DT
    exact = CBAR + AMP * radial * axial * np.exp(-gamma * T_END)
    return float(
        np.sqrt(np.sum(grid.volumes * (C - exact) ** 2) / np.sum(grid.volumes))
        / AMP
    )


def test_diffusion_eigenmode_second_order():
    errors = [eigenmode_error(n) for n in (8, 16, 32)]
    for n, err in zip((8, 16, 32), errors):
        print(f"\n{n}x{2*n}: L2/AMP={err:.4e}")
    for coarse, fine in zip(errors, errors[1:]):
        assert coarse / fine > 3.0
