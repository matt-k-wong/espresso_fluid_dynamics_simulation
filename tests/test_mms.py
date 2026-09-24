"""Task D.2 regression test: nontrivial 2D manufactured solution.

Manufactured pressure field (satisfies every boundary condition exactly:
top/bottom Dirichlet via sin(0)=sin(pi)=0, zero radial gradient at the axis
and sidewall since f'(0)=f'(R)=0):

    p(r,z) = p_in*(1-z/L) + A*f(r)*sin(a*z),  f = (1-(r/R)^2)^2,  a = pi/L.

The cell source is the *exact* integral of the continuum cylindrical
Laplacian over each finite volume -- it is derived below from calculus, not
by applying the code's discrete matrix to sampled pressures:

    s_ij = -(k0/mu) * A * ( Lap_r[f]_ij - a^2 * F_ij ) * S_j,

with  Lap_r[f]_ij = 2*pi*[r*f'] evaluated on the radial faces
      (since integral of (f''+f'/r)*r dr = [r*f']),
      F_ij   = integral of f*2*pi*r dr over the radial cell
             = pi*[r^2 - r^4/R^2 + r^6/(3*R^4)] on faces,
      S_j    = integral of sin(a*z) dz over the axial cell
             = (cos(a*z_w)-cos(a*z_e))/a,
      r*f'   = -4*r^2/R^2 + 4*r^4/R^4.

Volume-weighted L2 pressure error over amplitude must shrink by more than
3.8 per refinement on 8x16, 16x32, 32x64 (second-order scheme).
"""

import numpy as np
from dataclasses import replace

from espresso_m1.assembly import assemble
from espresso_m1.grid import make_grid
from espresso_m1.params import SimParams
from espresso_m1.solver import solve_linear


def manufactured_error(n: int) -> float:
    par = replace(SimParams(), Nr=n, Nz=2 * n, permeability_mode="uniform")
    grid = make_grid(par.Nr, par.Nz, par.R, par.L, par.R_perf)
    a = np.pi / par.L
    amplitude = par.p_in / 5.0
    rf = grid.r_faces
    z_edges = np.linspace(0.0, par.L, par.Nz + 1)

    # Exact radial integrals on faces.
    rfp = -4.0 * rf**2 / par.R**2 + 4.0 * rf**4 / par.R**4  # r*f'(r)
    face_f = np.pi * (
        rf**2 - rf**4 / par.R**2 + rf**6 / (3.0 * par.R**4)
    )  # primitive of f*2*pi*r
    int_sin = (np.cos(a * z_edges[:-1]) - np.cos(a * z_edges[1:])) / a

    lap_term = 2.0 * np.pi * np.diff(rfp)  # integral of (f''+f'/r) dV per dz
    f_term = np.diff(face_f)  # integral of f dV per dz
    source = (
        -par.k0 / par.mu
        * amplitude
        * (lap_term - a * a * f_term)[:, None]
        * int_sin[None, :]
    )

    A, b = assemble(par, grid, np.full((par.Nr, par.Nz), par.k0))
    numeric = solve_linear(A, b + source.T.ravel()).reshape(par.Nz, par.Nr).T
    exact = par.p_in * (1.0 - grid.z_centers[None, :] / par.L) + amplitude * (
        1.0 - (grid.r_centers[:, None] / par.R) ** 2
    ) ** 2 * np.sin(a * grid.z_centers[None, :])
    return float(
        np.sqrt(np.sum(grid.volumes * (numeric - exact) ** 2)
                / np.sum(grid.volumes)) / amplitude
    )


def test_manufactured_2d_second_order_accuracy():
    errors = [manufactured_error(n) for n in (8, 16, 32)]
    for n, err in zip((8, 16, 32), errors):
        print(f"\n{n}x{2*n}: L2/amplitude={err:.4e}")
    for coarse, fine in zip(errors, errors[1:]):
        assert coarse / fine > 3.8
