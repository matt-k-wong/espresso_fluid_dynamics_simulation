"""Permeability closures (milestone 1).

Two modes are offered:

1. ``uniform`` -- constant permeability ``k = k0`` everywhere. Used for
   verification against the analytical cylindrical Darcy solution.

2. ``compaction`` -- a reduced, pressure-dependent hydraulic closure:

   a. Effective stress from a prescribed top load plus Terzaghi's principle
      with a constant total axial stress::

          sigma_eff(p) = sigma_eff_top + p_in - p.

      Load assumption: the total axial stress is uniform over the bed and
      equal to its value at the top boundary, ``sigma_total = sigma_eff_top
      + p_in`` (Biot coefficient 1, quasi-static, uniaxial, no wall
      friction, no shear). The local effective stress is then the total
      stress minus the local pore pressure. This is a *reduced hydraulic
      closure*, not a mechanical deformation model: the bed geometry is
      fixed and no displacement, strain, or solid mass balance is solved
      (those belong to milestone B, which requires separate approval).

   b. Bounded porosity-effective-stress law::

          phi(sigma) = phi_min + (phi0 - phi_min) * exp(-sigma_+ / sigma_c),

      with ``sigma_+ = max(sigma_eff, 0)`` (no tensile states; the bed
      cannot be looser than its poured reference state). Porosity is hence
      confined to ``[phi_min, phi0]`` with ``0 < phi_min < phi0 < 1``.

   c. Carman-Kozeny relation anchored at the reference state::

          k(phi) = k0 * [phi^3/(1-phi)^2] / [phi0^3/(1-phi0)^2].

      Because ``phi >= phi_min > 0``, permeability is strictly positive.

All compaction coefficients are *illustrative assumptions* (see the
parameter catalog): they are order-of-magnitude placeholders, not measured
coffee properties, and must be calibrated against experiments (see REPORT).
"""

from __future__ import annotations

import numpy as np

from .params import SimParams


def effective_stress(
    p: np.ndarray, p_in: float, sigma_eff_top: float
) -> np.ndarray:
    """Terzaghi effective stress under constant total axial stress.

    sigma_eff = sigma_eff_top + p_in - p, clipped below at zero.
    """
    return np.maximum(sigma_eff_top + p_in - p, 0.0)


def porosity_from_stress(
    sigma: np.ndarray, phi0: float, phi_min: float, sigma_c: float
) -> np.ndarray:
    """Bounded exponential porosity law; output lies in [phi_min, phi0]."""
    phi = phi_min + (phi0 - phi_min) * np.exp(-np.maximum(sigma, 0.0) / sigma_c)
    return np.clip(phi, phi_min, phi0)


def kozeny_ratio(phi: np.ndarray, phi_ref: float) -> np.ndarray:
    """k(phi)/k(phi_ref) from the Carman-Kozeny porosity function."""
    f = phi**3 / (1.0 - phi) ** 2
    f_ref = phi_ref**3 / (1.0 - phi_ref) ** 2
    return f / f_ref


def cell_permeability(
    p: np.ndarray, params: SimParams
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``(k, phi, sigma_eff)`` cell fields for a pressure field.

    Shapes are ``(Nr, Nz)``. In ``uniform`` mode ``k = k0`` and ``phi = phi0``
    everywhere (``sigma_eff`` is still reported from the stress relation for
    transparency but does not feed back into ``k``).
    """
    if params.permeability_mode == "uniform":
        k = np.full_like(p, params.k0, dtype=float)
        phi = np.full_like(p, params.phi0, dtype=float)
        sigma = effective_stress(p, params.p_in, params.sigma_eff_top)
        return k, phi, sigma
    if params.permeability_mode == "compaction":
        sigma = effective_stress(p, params.p_in, params.sigma_eff_top)
        phi = porosity_from_stress(sigma, params.phi0, params.phi_min, params.sigma_c)
        k = params.k0 * kozeny_ratio(phi, params.phi0)
        # Defensive floor: the porosity bound already guarantees k > 0, but
        # never allow a non-positive permeability to reach the assembler.
        k = np.maximum(k, 1.0e-30)
        return k, phi, sigma
    raise ValueError(f"Unknown permeability_mode {params.permeability_mode!r}")
