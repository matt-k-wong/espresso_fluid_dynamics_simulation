"""Post-processing diagnostics: flow rates, balances, resistance."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from .grid import Grid
from .params import SimParams
from .solver import FLOW_ATOL_M3S, NonlinearResult

ML_PER_M3 = 1.0e6  # m^3 -> mL


@dataclass
class Diagnostics:
    Q_in_m3s: float  # total inflow through the top [m^3/s]
    Q_out_m3s: float  # total outflow through the perforations [m^3/s]
    Q_out_mLs: float  # outlet flow [mL/s]
    mismatch_abs_m3s: float  # |Q_in - Q_out| [m^3/s]
    mismatch_rel: float  # |Q_in - Q_out| / max(|Q_out|, 1e-18) [-]
    max_residual_abs_m3s: float  # max |cell net outflow| [m^3/s]
    max_residual_rel: float  # max residual / max(|Q_out|, 1e-18) [-]
    delta_p_Pa: float  # p_in - p_out [Pa]
    R_hyd: Optional[float]  # hydraulic resistance delta_p / Q_out, None if Q<=0
    R_hyd_status: str  # "ok" or "undefined_zero_flow"
    k_min: float
    k_max: float
    phi_min: float
    phi_max: float
    p_min: float
    p_max: float
    # Physical balance of the returned state (acceptance scale).
    q_scale_m3s: float  # max(|Q_in|, |Q_out|) [m^3/s]
    balance_limit_m3s: float  # FLOW_ATOL_M3S + tol_res_rel * q_scale [m^3/s]
    open_area_m2: float  # represented (discrete) open area [m^2]
    open_area_ideal_m2: float  # requested open area pi*R_perf^2 [m^2]


def evaluate(
    params: SimParams, grid: Grid, sol: NonlinearResult
) -> Diagnostics:
    Q_in = float(np.sum(sol.Fz[:, 0]))
    Q_out = float(np.sum(sol.Fz[grid.open_mask, grid.Nz]))
    mismatch_abs = abs(Q_in - Q_out)
    denom = max(abs(Q_out), 1.0e-18)
    mismatch_rel = mismatch_abs / denom

    # Local finite-volume residuals: net outward flux per cell.
    res = (
        (sol.Fr[1:, :] - sol.Fr[:-1, :])
        + (sol.Fz[:, 1:] - sol.Fz[:, :-1])
    )
    max_res_abs = float(np.max(np.abs(res)))
    max_res_rel = max_res_abs / denom

    delta_p = params.p_in - params.p_out
    if Q_out > 0.0:
        R_hyd: Optional[float] = delta_p / Q_out
        R_hyd_status = "ok"
    else:
        R_hyd = None
        R_hyd_status = "undefined_zero_flow"
    q_scale = max(abs(Q_in), abs(Q_out))
    return Diagnostics(
        Q_in_m3s=Q_in,
        Q_out_m3s=Q_out,
        Q_out_mLs=Q_out * ML_PER_M3,
        mismatch_abs_m3s=mismatch_abs,
        mismatch_rel=mismatch_rel,
        max_residual_abs_m3s=max_res_abs,
        max_residual_rel=max_res_rel,
        delta_p_Pa=delta_p,
        R_hyd=R_hyd,
        R_hyd_status=R_hyd_status,
        k_min=float(np.min(sol.k)),
        k_max=float(np.max(sol.k)),
        phi_min=float(np.min(sol.phi)),
        phi_max=float(np.max(sol.phi)),
        p_min=float(np.min(sol.p)),
        p_max=float(np.max(sol.p)),
        q_scale_m3s=q_scale,
        balance_limit_m3s=FLOW_ATOL_M3S + params.tol_res_rel * q_scale,
        open_area_m2=grid.open_area,
        open_area_ideal_m2=grid.open_area_ideal,
    )
