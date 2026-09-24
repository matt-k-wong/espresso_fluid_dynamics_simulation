"""Axisymmetric cell-centered finite-volume grid.

Coordinates: r in [0, R], z in [0, L], z increasing downward.
All face areas and cell volumes include the axisymmetric 2*pi factor, so
face fluxes are true volumetric rates in m^3/s.

Mixed bottom boundary (perforated plate with solid rim)
------------------------------------------------------
The bottom (z = L) is split into a perforated part held at the outlet gauge
pressure ``p_out`` (Dirichlet) and a solid rim with zero normal flux
(homogeneous Neumann). The split is a stair-step approximation on the
finite-volume grid: a radial column ``i`` is *open* when its cell-center
radius satisfies ``r_center[i] < R_perf`` and *closed* otherwise.

- Open column: the bottom half-face carries a Dirichlet flux
  ``F = T_bot * (p_cell - p_out)`` with ``T_bot = (k_cell/mu) * A/(dz/2)``.
  This is the standard cell-centered treatment of a Dirichlet boundary
  (half-cell distance from center to boundary).
- Closed column: the bottom face flux is set identically to zero.

The open area is ``sum(col_areas[open])`` and approximates ``pi*R_perf^2``
to first order in ``dr``; the geometric error is reported in the
``Grid.open_area`` / ``Grid.open_area_ideal`` attributes and checked in the
test suite. No partial-face weighting is used in milestone 1.

Solid-rim approximation
-----------------------
A closed bottom column imposes zero normal flux on its whole bottom face, so
the Dirichlet/Neumann junction follows column boundaries (a stair-step in
the continuum picture). For ``R_perf = R/2`` with even ``Nr`` the junction
lies exactly on a grid face and the discrete open area is exact; the reduced
solution regularity at the junction still limits convergence there, which is
a discretization effect distinct from any area error. Geometric
(requested vs represented open area) and discretization errors are reported
separately.

A positive ``R_perf`` that resolves to zero open columns is *not* a sealed
bed: it is an under-resolved request and is rejected with an actionable
error (increase ``Nr`` or ``R_perf``). Only ``R_perf = 0`` selects the
intentionally sealed bed (uniform ``p = p_in``, zero flow).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Grid:
    Nr: int
    Nz: int
    R: float
    L: float
    R_perf: float
    dr: float
    dz: float
    r_centers: np.ndarray  # (Nr,) cell-center radii
    z_centers: np.ndarray  # (Nz,) cell-center depths
    r_faces: np.ndarray  # (Nr+1,) radial face positions
    col_areas: np.ndarray  # (Nr,) axial (top/bottom) face area per column [m^2]
    rad_areas: np.ndarray  # (Nr+1,) radial face areas (full 2*pi*r*dz) [m^2]
    volumes: np.ndarray  # (Nr, Nz) cell volumes [m^3]
    open_mask: np.ndarray  # (Nr,) bool, True where the bottom is perforated
    open_area: float  # discrete open area [m^2]
    open_area_ideal: float  # pi * min(R_perf, R)^2 [m^2]
    bed_area: float  # pi * R^2 [m^2]


def make_grid(Nr: int, Nz: int, R: float, L: float, R_perf: float) -> Grid:
    dr = R / Nr
    dz = L / Nz
    r_faces = np.linspace(0.0, R, Nr + 1)
    r_centers = (np.arange(Nr) + 0.5) * dr
    z_centers = (np.arange(Nz) + 0.5) * dz
    # Axial face area of column i: annulus pi*(re^2 - rw^2).
    col_areas = np.pi * (r_faces[1:] ** 2 - r_faces[:-1] ** 2)
    # Radial face areas (axisymmetric shells).
    rad_areas = 2.0 * np.pi * r_faces * dz
    volumes = col_areas[:, None] * np.ones((1, Nz)) * dz
    open_mask = r_centers < R_perf
    if R_perf > 0.0 and not np.any(open_mask):
        raise ValueError(
            f"R_perf={R_perf:.6e} m opens no bottom columns on this grid "
            f"(Nr={Nr}, dr={dr:.6e} m, innermost cell center r={r_centers[0]:.6e} m, "
            f"requested open area pi*R_perf^2={np.pi * R_perf**2:.6e} m^2). "
            "This is an under-resolved opening, not a sealed bed: increase Nr "
            "(finer dr) or increase R_perf so that at least the innermost "
            "cell center satisfies r_center < R_perf. Use R_perf=0 for an "
            "intentionally sealed bed."
        )
    open_area = float(np.sum(col_areas[open_mask]))
    open_area_ideal = float(np.pi * min(R_perf, R) ** 2)
    return Grid(
        Nr=Nr,
        Nz=Nz,
        R=R,
        L=L,
        R_perf=R_perf,
        dr=dr,
        dz=dz,
        r_centers=r_centers,
        z_centers=z_centers,
        r_faces=r_faces,
        col_areas=col_areas,
        rad_areas=rad_areas,
        volumes=volumes,
        open_mask=open_mask,
        open_area=open_area,
        open_area_ideal=open_area_ideal,
        bed_area=float(np.pi * R**2),
    )


def cell_index(Nr: int, i: int, j: int) -> int:
    """Row-major flat index for cell (i, j): idx = j*Nr + i."""
    return j * Nr + i
