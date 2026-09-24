"""High-level drivers: single cases, pressure sweeps, analytical reference.

This module is free of plotting and file output by design so the numerical
core can later be ported (e.g. to Rust); the CLI layer owns all I/O.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import math

import numpy as np

from .diagnostics import Diagnostics, evaluate
from .grid import Grid, make_grid
from .params import BAR_PA, SimParams, validate
from .solver import ConvergenceError, NonlinearResult, solve_pressure


@dataclass
class CaseResult:
    params: SimParams
    grid: Grid
    sol: NonlinearResult | None
    diag: Diagnostics | None
    converged: bool
    iters: int
    error: str = ""


@dataclass
class SweepRow:
    p_in_Pa: float
    p_in_bar: float
    converged: bool
    iters: int
    Q_m3s: float
    Q_mLs: float
    mismatch_abs_m3s: float
    mismatch_rel: float
    max_residual_abs_m3s: float
    max_residual_rel: float
    delta_p_Pa: float
    R_hyd: float  # NaN when undefined (zero flow); converged flag + error say why
    R_hyd_status: str
    k_min: float
    k_max: float
    phi_min: float
    phi_max: float
    open_area_m2: float
    open_area_ideal_m2: float
    error: str = ""


def analytical_cylindrical_Q(R: float, k: float, delta_p: float, mu: float, L: float) -> float:
    """1D cylindrical Darcy result Q = pi R^2 k (p_in - p_out) / (mu L)."""
    return float(np.pi * R**2 * k * delta_p / (mu * L))


def run_case(params: SimParams) -> CaseResult:
    """Run one steady case; failed convergence is reported, never hidden."""
    validate(params)
    grid = make_grid(params.Nr, params.Nz, params.R, params.L, params.R_perf)
    try:
        sol = solve_pressure(params, grid)
    except ConvergenceError as exc:
        return CaseResult(
            params=params,
            grid=grid,
            sol=None,
            diag=None,
            converged=False,
            iters=exc.iters,
            error=str(exc),
        )
    diag = evaluate(params, grid, sol)
    return CaseResult(
        params=params, grid=grid, sol=sol, diag=diag, converged=True,
        iters=sol.iters,
    )


def run_sweep(base: SimParams, p_in_values: list[float]) -> list[SweepRow]:
    """Run a pressure sweep over inlet pressures (Pa, gauge).

    Every point is attempted independently from the same initial guess;
    non-converged points are recorded with ``converged=False`` and NaN
    flow/diagnostic fields so failures stay visible in the CSV. Only
    ``ConvergenceError`` is treated as an ordinary point failure; programming
    errors propagate. An empty pressure list is rejected (it cannot succeed).
    """
    if len(p_in_values) == 0:
        raise ValueError("Empty sweep: no inlet pressures requested")
    for value in p_in_values:
        if not math.isfinite(value) or value < 0.0:
            raise ValueError(
                f"Invalid sweep pressure {value!r}: inlet gauge pressures "
                "must be finite and non-negative"
            )
    rows: list[SweepRow] = []
    for p_in in p_in_values:
        params = replace(base, p_in=float(p_in))
        case = run_case(params)
        if case.converged:
            assert case.diag is not None
            d = case.diag
            rows.append(
                SweepRow(
                    p_in_Pa=float(p_in),
                    p_in_bar=float(p_in) / BAR_PA,
                    converged=True,
                    iters=case.iters,
                    Q_m3s=d.Q_out_m3s,
                    Q_mLs=d.Q_out_mLs,
                    mismatch_abs_m3s=d.mismatch_abs_m3s,
                    mismatch_rel=d.mismatch_rel,
                    max_residual_abs_m3s=d.max_residual_abs_m3s,
                    max_residual_rel=d.max_residual_rel,
                    delta_p_Pa=d.delta_p_Pa,
                    R_hyd=d.R_hyd if d.R_hyd is not None else float("nan"),
                    R_hyd_status=d.R_hyd_status,
                    k_min=d.k_min,
                    k_max=d.k_max,
                    phi_min=d.phi_min,
                    phi_max=d.phi_max,
                    open_area_m2=d.open_area_m2,
                    open_area_ideal_m2=d.open_area_ideal_m2,
                )
            )
        else:
            rows.append(
                SweepRow(
                    p_in_Pa=float(p_in),
                    p_in_bar=float(p_in) / BAR_PA,
                    converged=False,
                    iters=case.iters,
                    Q_m3s=float("nan"),
                    Q_mLs=float("nan"),
                    mismatch_abs_m3s=float("nan"),
                    mismatch_rel=float("nan"),
                    max_residual_abs_m3s=float("nan"),
                    max_residual_rel=float("nan"),
                    delta_p_Pa=float(p_in) - base.p_out,
                    R_hyd=float("nan"),
                    R_hyd_status="no_solution",
                    k_min=float("nan"),
                    k_max=float("nan"),
                    phi_min=float("nan"),
                    phi_max=float("nan"),
                    open_area_m2=case.grid.open_area,
                    open_area_ideal_m2=case.grid.open_area_ideal,
                    error=case.error,
                )
            )
    return rows
