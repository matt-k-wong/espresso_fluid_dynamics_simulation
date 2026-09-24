"""Finite-volume assembly for steady saturated Darcy flow.

Discrete equations
------------------
For each cell, the sum of outward volumetric face fluxes is zero::

    (Fr_e - Fr_w) + (Fz_s - Fz_n) = 0,

with Darcy fluxes ``F = T * (p_here - p_there)`` and transmissibility
``T = (k_face / mu) * A_face / d``.

Face mobility (harmonic averaging)
----------------------------------
Interior faces use the harmonic average of the neighboring cell
permeabilities::

    k_face = 2 / (1/k_L + 1/k_R).

Rationale: for one-dimensional steady flow through two half-cells in
series, the harmonic average is the *exact* equivalent permeability that
preserves flux continuity (equal pressure at the shared face). It stays
bounded by twice the smaller value, so a compacted low-k cell correctly
throttles the face flux, whereas an arithmetic average would let a
high-k neighbor leak through the low-k cell. Boundary half-faces use the
adjacent cell value (single material between center and boundary).

Linear solver choice
--------------------
The assembled matrix is sparse with a 5-point stencil, symmetric positive
definite whenever at least one Dirichlet face exists (top inlet and any
open bottom column anchor the pressure; pure-Neumann rows are absent
because axis/sidewall/rim fluxes are identically zero rather than
extra unknowns). Milestone-1 grids are small (default 20x40 = 800
unknowns; refinement studies stay below ~10^4), so the sparse direct
solver ``scipy.sparse.linalg.spsolve`` (SuperLU) is used: it is robust,
deterministic, and needs no preconditioner tuning. Larger or 3D problems
would call for CG with an SPD-compatible preconditioner (e.g. algebraic
multigrid or incomplete Cholesky); generic ILU is not symmetric and is not
suitable for CG. That is out of scope.
"""

from __future__ import annotations

import numpy as np
from scipy import sparse

from .grid import Grid, cell_index
from .params import SimParams


def harmonic_face_k(k_left: float, k_right: float) -> float:
    """Harmonic average of two strictly positive cell permeabilities."""
    if not (k_left > 0.0 and k_right > 0.0):
        raise ValueError("Permeabilities must be strictly positive")
    return 2.0 / (1.0 / k_left + 1.0 / k_right)


def radial_face_T(
    grid: Grid, mu: float, k: np.ndarray, i_face: int, j: int
) -> float:
    """Transmissibility of interior radial face ``i_face`` (1..Nr-1)."""
    kf = harmonic_face_k(k[i_face - 1, j], k[i_face, j])
    return kf / mu * grid.rad_areas[i_face] / grid.dr


def axial_face_T(
    grid: Grid, mu: float, k: np.ndarray, i: int, j_face: int
) -> float:
    """Transmissibility of interior axial face ``j_face`` (1..Nz-1)."""
    kf = harmonic_face_k(k[i, j_face - 1], k[i, j_face])
    return kf / mu * grid.col_areas[i] / grid.dz


def top_face_T(grid: Grid, mu: float, k: np.ndarray, i: int) -> float:
    """Transmissibility of the top Dirichlet half-face in column ``i``."""
    return k[i, 0] / mu * grid.col_areas[i] / (grid.dz / 2.0)


def bottom_face_T(grid: Grid, mu: float, k: np.ndarray, i: int) -> float:
    """Transmissibility of the bottom Dirichlet half-face in column ``i``."""
    return k[i, grid.Nz - 1] / mu * grid.col_areas[i] / (grid.dz / 2.0)


def assemble(
    params: SimParams, grid: Grid, k: np.ndarray
) -> tuple[sparse.csr_matrix, np.ndarray]:
    """Assemble the sparse linear system ``A p = b`` for cell pressures.

    ``k`` is the ``(Nr, Nz)`` cell permeability field. Returns ``(A_csr, b)``
    with ``p`` ordered by ``idx = j*Nr + i``.
    """
    if np.any(~np.isfinite(k)) or np.any(k <= 0.0):
        raise ValueError("Permeability field must be finite and positive")
    Nr, Nz = grid.Nr, grid.Nz
    mu = params.mu
    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []
    b = np.zeros(Nr * Nz)
    for j in range(Nz):
        for i in range(Nr):
            r = cell_index(Nr, i, j)
            diag = 0.0
            # West radial neighbor (i = 0 column has no west face: the axis
            # face area is identically zero, hence zero flux).
            if i > 0:
                T = radial_face_T(grid, mu, k, i, j)
                diag += T
                rows.append(r)
                cols.append(cell_index(Nr, i - 1, j))
                data.append(-T)
            # East radial neighbor (sidewall: no flux, nothing added).
            if i < Nr - 1:
                T = radial_face_T(grid, mu, k, i + 1, j)
                diag += T
                rows.append(r)
                cols.append(cell_index(Nr, i + 1, j))
                data.append(-T)
            # North: neighbor above, or top Dirichlet (inlet pressure).
            if j > 0:
                T = axial_face_T(grid, mu, k, i, j)
                diag += T
                rows.append(r)
                cols.append(cell_index(Nr, i, j - 1))
                data.append(-T)
            else:
                T = top_face_T(grid, mu, k, i)
                diag += T
                b[r] += T * params.p_in
            # South: neighbor below, open-bottom Dirichlet, or solid rim.
            if j < Nz - 1:
                T = axial_face_T(grid, mu, k, i, j + 1)
                diag += T
                rows.append(r)
                cols.append(cell_index(Nr, i, j + 1))
                data.append(-T)
            elif grid.open_mask[i]:
                T = bottom_face_T(grid, mu, k, i)
                diag += T
                b[r] += T * params.p_out
            # else: solid rim, zero normal flux -> no contribution.
            rows.append(r)
            cols.append(r)
            data.append(diag)
    A = sparse.coo_matrix((data, (rows, cols)), shape=(Nr * Nz, Nr * Nz)).tocsr()
    return A, b


def compute_face_fluxes(
    params: SimParams, grid: Grid, k: np.ndarray, p: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Volumetric face fluxes in m^3/s from cell pressures.

    Returns ``(Fr, Fz)``:

    - ``Fr`` shape ``(Nr+1, Nz)``: radial fluxes, positive in +r.
      ``Fr[0, :]`` (axis) and ``Fr[Nr, :]`` (sidewall) are identically zero.
    - ``Fz`` shape ``(Nr, Nz+1)``: axial fluxes, positive downward (+z).
      ``Fz[:, 0]`` is the top inflow, ``Fz[:, Nz]`` the bottom outflow
      (zero on the solid rim by construction).
    """
    Nr, Nz = grid.Nr, grid.Nz
    mu = params.mu
    Fr = np.zeros((Nr + 1, Nz))
    Fz = np.zeros((Nr, Nz + 1))
    for j in range(Nz):
        for i_face in range(1, Nr):
            T = radial_face_T(grid, mu, k, i_face, j)
            Fr[i_face, j] = T * (p[i_face - 1, j] - p[i_face, j])
    for i in range(Nr):
        for j_face in range(1, Nz):
            T = axial_face_T(grid, mu, k, i, j_face)
            Fz[i, j_face] = T * (p[i, j_face - 1] - p[i, j_face])
    for i in range(Nr):
        Fz[i, 0] = top_face_T(grid, mu, k, i) * (params.p_in - p[i, 0])
        if grid.open_mask[i]:
            Fz[i, Nz] = bottom_face_T(grid, mu, k, i) * (p[i, Nz - 1] - params.p_out)
        else:
            Fz[i, Nz] = 0.0
    return Fr, Fz
