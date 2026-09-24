"""Nonlinear (Picard) solver for the pressure/permeability coupling.

Iteration contract (repaired)
-----------------------------
1. Pressure is solved with the *current iterate* of ``k``. Relaxation only
   shapes the next iterate; it never enters the definition of the residual.
2. For every pressure candidate, the **unrelaxed** physical fields
   ``k(p)``, ``phi(p)``, ``sigma(p)`` are evaluated.
3. The candidate's actual face fluxes are computed with that ``k(p)``; local
   net outward flux per cell (m^3/s) and the global inlet-minus-outlet
   mismatch (m^3/s) are evaluated from those fluxes.
4. A candidate is accepted only when the pressure-update criterion, the
   flow-update criterion, AND the physical local/global balance criteria all
   pass. The returned state is exactly the checked state (no unchecked polish
   solve follows acceptance).
5. Convergence diagnostics are populated from the returned state and the
   actual consecutive candidates used in its acceptance.

Tolerances (combined absolute + relative, documented units)
-----------------------------------------------------------
- Pressure update: ``dp_abs <= tol_p_atol_Pa + tol_p_rel * |p_in - p_out|``.
  The relative scale is the *imposed pressure difference*, so a common
  pressure offset on both boundaries cannot change acceptance.
  Reported ``dp_rel = dp_abs / max(|p_in - p_out|, tol_p_atol_Pa)``.
- Flow update: ``dq_abs <= FLOW_ATOL_M3S + tol_q_rel * max(|Q_new|, |Q_prev|)``
  with physical outlet flows. Reported ``dq_rel`` uses the same denominator
  with ``FLOW_ATOL_M3S`` as floor.
- Balances: with ``Q_scale = max(|Q_in|, |Q_out|)`` of the returned state,
  ``balance_limit = FLOW_ATOL_M3S + tol_res_rel * Q_scale``; each cell's
  absolute net flux and the absolute global mismatch must be below it.
  Reported ``res_rel = max(local_abs, mismatch_abs) / max(Q_scale,
  FLOW_ATOL_M3S)``.
- ``FLOW_ATOL_M3S = 1e-18`` m^3/s is an explicit *numerical* default (absolute
  floor against division by ~zero flux), not a physical parameter.

Shifted solve
-------------
The linear system is solved for ``u = p - p_out`` (top Dirichlet
``p_in - p_out``, outlet Dirichlet ``0``) and shifted back afterwards. The
matrix is unchanged; the unknown is O(pressure difference) instead of
O(absolute pressure), which reduces floating-point cancellation when both
boundary pressures carry a large common offset. Physics is offset-invariant
(``sigma_eff`` depends on ``p_in - p``), so only numerics change.

In ``uniform`` mode the system is linear: the first solve gives the exact
discrete solution and the second confirms it, so convergence takes 2 linear
solves. In ``compaction`` mode permeability depends on pressure and the
Picard fixed point above is iterated. If the iteration budget is exhausted,
a ``ConvergenceError`` carrying the last indicators is raised -- the last
iterate is *never* silently returned as a solution.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from scipy.sparse.linalg import spsolve

from .assembly import assemble, compute_face_fluxes
from .closures import cell_permeability
from .grid import Grid
from .params import SimParams

FLOW_ATOL_M3S = 1.0e-18  # absolute floor for flow/balance criteria [m^3/s]


class ConvergenceError(RuntimeError):
    """Raised when the Picard iteration fails to converge.

    Attributes: iters, dp_rel, dq_rel, res_rel (last values), message.
    """

    def __init__(
        self, message: str, iters: int, dp_rel: float, dq_rel: float, res_rel: float
    ) -> None:
        super().__init__(message)
        self.iters = iters
        self.dp_rel = dp_rel
        self.dq_rel = dq_rel
        self.res_rel = res_rel


@dataclass
class NonlinearResult:
    p: np.ndarray  # (Nr, Nz) cell pressures [Pa, gauge]
    k: np.ndarray  # (Nr, Nz) cell permeabilities [m^2]
    phi: np.ndarray  # (Nr, Nz) cell porosities [-]
    sigma_eff: np.ndarray  # (Nr, Nz) effective stresses [Pa]
    Fr: np.ndarray  # (Nr+1, Nz) radial face fluxes [m^3/s, +r]
    Fz: np.ndarray  # (Nr, Nz+1) axial face fluxes [m^3/s, +z downward]
    iters: int  # linear solves performed
    dp_rel: float  # pressure-update indicator, scale max(|dp_top|, atol) [-]
    dq_rel: float  # flow-update indicator, scale max(|Q|, atol) [-]
    res_rel: float  # balance indicator, scale max(Q_scale, atol) [-]
    dp_abs_Pa: float  # max |p_candidate - p_previous| [Pa]
    dp_limit_Pa: float  # combined pressure acceptance limit [Pa]
    dq_abs_m3s: float  # |Q_out_candidate - Q_out_previous| [m^3/s]
    local_abs_m3s: float  # max |cell net outward flux| of returned state
    mismatch_abs_m3s: float  # |Q_in - Q_out| of returned state [m^3/s]
    q_scale_m3s: float  # max(|Q_in|, |Q_out|) of returned state [m^3/s]
    balance_limit_m3s: float  # combined balance acceptance limit [m^3/s]


def initial_pressure_guess(params: SimParams, grid: Grid) -> np.ndarray:
    """Linear-in-z profile between p_in and p_out, tiled over radius."""
    profile = params.p_in + (params.p_out - params.p_in) * (
        grid.z_centers / params.L
    )
    return np.tile(profile[None, :], (grid.Nr, 1))


def solve_linear(A, b: np.ndarray) -> np.ndarray:
    """Sparse direct solve (SuperLU via ``spsolve``); see module notes."""
    p = spsolve(A.tocsr(), b)
    p = np.asarray(p, dtype=float)
    if np.any(~np.isfinite(p)):
        raise RuntimeError("Linear solver returned non-finite pressures")
    return p


def solve_pressure(params: SimParams, grid: Grid) -> NonlinearResult:
    """Solve the (possibly nonlinear) discrete Darcy problem."""
    dp_top = params.p_in - params.p_out
    shifted = replace(params, p_in=dp_top, p_out=0.0)

    p_prev = initial_pressure_guess(params, grid)
    k_prev, _, _ = cell_permeability(p_prev, params)
    Q_prev: float | None = None

    dp_abs = float("inf")
    dq_abs = float("inf")
    local_abs = float("inf")
    mismatch_abs = float("inf")
    q_scale = 0.0

    it = 0
    while it < params.max_iter:
        it += 1
        # 1. Solve with the current k iterate (relaxation never enters here).
        A, b = assemble(shifted, grid, k_prev)
        u_new = solve_linear(A, b).reshape(grid.Nz, grid.Nr).T
        p_new = u_new + params.p_out

        # 2-3. Unrelaxed physical fields and the candidate's actual fluxes.
        k_phys, phi_phys, sig_phys = cell_permeability(p_new, params)
        Fr, Fz = compute_face_fluxes(params, grid, k_phys, p_new)
        Q_in = float(np.sum(Fz[:, 0]))
        Q_out = float(np.sum(Fz[grid.open_mask, grid.Nz]))
        cell_net = (Fr[1:, :] - Fr[:-1, :]) + (Fz[:, 1:] - Fz[:, :-1])
        local_abs = float(np.max(np.abs(cell_net)))
        mismatch_abs = abs(Q_in - Q_out)
        q_scale = max(abs(Q_in), abs(Q_out))

        dp_abs = float(np.max(np.abs(p_new - p_prev)))
        dp_limit = params.tol_p_atol_Pa + params.tol_p_rel * abs(dp_top)
        if Q_prev is None:
            dq_abs = float("inf")
        else:
            dq_abs = abs(Q_out - Q_prev)
        dq_denom = max(abs(Q_out), abs(Q_prev) if Q_prev is not None else 0.0)
        dq_limit = FLOW_ATOL_M3S + params.tol_q_rel * dq_denom
        balance_limit = FLOW_ATOL_M3S + params.tol_res_rel * q_scale

        # 4. Accept only the checked state: updates AND physical balances.
        if (
            Q_prev is not None
            and dp_abs <= dp_limit
            and dq_abs <= dq_limit
            and local_abs <= balance_limit
            and mismatch_abs <= balance_limit
        ):
            return NonlinearResult(
                p=p_new,
                k=k_phys,
                phi=phi_phys,
                sigma_eff=sig_phys,
                Fr=Fr,
                Fz=Fz,
                iters=it,
                dp_rel=dp_abs / max(abs(dp_top), params.tol_p_atol_Pa),
                dq_rel=dq_abs / max(dq_denom, FLOW_ATOL_M3S),
                res_rel=max(local_abs, mismatch_abs)
                / max(q_scale, FLOW_ATOL_M3S),
                dp_abs_Pa=dp_abs,
                dp_limit_Pa=dp_limit,
                dq_abs_m3s=dq_abs,
                local_abs_m3s=local_abs,
                mismatch_abs_m3s=mismatch_abs,
                q_scale_m3s=q_scale,
                balance_limit_m3s=balance_limit,
            )

        # 5. Relaxation shapes the next iterate only.
        w = params.relaxation
        k_prev = w * k_phys + (1.0 - w) * k_prev
        p_prev, Q_prev = p_new, Q_out

    raise ConvergenceError(
        f"Picard iteration did not converge in {params.max_iter} iterations "
        f"(dp_rel={dp_abs / max(abs(dp_top), params.tol_p_atol_Pa):.3e}, "
        f"dq_rel={dq_abs / max(dq_denom, FLOW_ATOL_M3S):.3e}, "
        f"res_rel={max(local_abs, mismatch_abs) / max(q_scale, FLOW_ATOL_M3S):.3e})",
        iters=params.max_iter,
        dp_rel=dp_abs / max(abs(dp_top), params.tol_p_atol_Pa),
        dq_rel=dq_abs / max(dq_denom, FLOW_ATOL_M3S),
        res_rel=max(local_abs, mismatch_abs) / max(q_scale, FLOW_ATOL_M3S),
    )
