"""Transient solute transport and lumped extraction on frozen hydraulics.

One hydraulically converged, saturated, isothermal snapshot is frozen per
run: grid, porosity, permeability, and face volume fluxes never change
during the transport run. Time zero means an already wetted bed.

Semidiscrete equations per cell i (W = phi*V fixed storage, passive
approximation: dissolving mass does not alter geometry or properties):

    W_i dC_i/dt + sum_outward(J_if) = R_i
    dB_i/dt = -R_i
    R_i = lambda * (B_i - K*m_i*C_i)                 [kg/s]

Reversible partition exchange; local equilibrium is B_i = K*m_i*C_i.
Inflowing liquid can raise a local B_i above its initial value; that is
redistribution, not a conservation bug, and B is never clipped. K = 0 is
the irreversible first-order release limit.

Interior face i -> j (F signed, contains face area already; never divided
by phi -- pore velocity u/phi is a travel-time reference, not a flux):

    C_up = C_i if F >= 0 else C_j
    J = F*C_up + H*(C_i - C_j),  H = D*harm(phi_i,phi_j)*A/d    [m^3/s]

Top: total inward solute flux F_in*C_in(t) per face (Danckwerts-type; no
separate Dirichlet diffusion flux). Outlet: outward advection F*C_cell with
zero dispersive gradient. Axis/sidewall/rim: zero total flux. Hydraulic flow
must be forward or zero; reversal is rejected, not accommodated.

Time update is backward Euler with first-order upwind and centered
dispersion, one sparse direct solve per step after algebraic elimination of
the solid pool (see ``transport_step``). Every step checks finite values,
linear residual, positivity (no corrective clipping), and the global solute
budget; failures are reported explicitly via ``TransportError``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

from .diagnostics import ML_PER_M3  # noqa: F401  (re-exported for CLI/tests)
from .driver import CaseResult, run_case
from .grid import Grid, cell_index
from .params import SimParams
from .solver import FLOW_ATOL_M3S
from .tparams import (
    TransportParams,
    inlet_concentration,
    schedule_breakpoints,
    validate,
)

LIN_REL_TOL = 1e-8  # linear-solve residual acceptance (relative, inf-norm)
POS_REL_TOL = 1e-9  # positivity acceptance relative to concentration scale
BUDGET_ATOL_KG = 1e-12  # global solute budget absolute floor [kg]
BUDGET_RTOL = 1e-9  # global solute budget relative factor [-]


class TransportError(RuntimeError):
    """Explicit failure of a transport run (unsupported state or bad step)."""


def _hmean(a: float, b: float) -> float:
    if not (a > 0.0 and b > 0.0):
        raise ValueError("Harmonic mean needs strictly positive values")
    return 2.0 / (1.0 / a + 1.0 / b)


@dataclass
class HydraulicSnapshot:
    """Frozen hydraulic input to a transport run."""

    grid: Grid
    phi: np.ndarray  # (Nr, Nz) frozen porosity [-]
    m: np.ndarray  # (Nr, Nz) dry coffee mass per cell [kg], sums to M_dose
    W: np.ndarray  # (Nr, Nz) liquid storage volume [m^3]
    Fr: np.ndarray  # frozen radial face volume fluxes [m^3/s]
    Fz: np.ndarray  # frozen axial face volume fluxes [m^3/s]
    Q_in: float  # frozen top inflow [m^3/s]
    Q_out: float  # frozen outlet outflow [m^3/s]


def freeze_snapshot(case: CaseResult, M_dose: float) -> HydraulicSnapshot:
    """Freeze a converged hydraulic case; reject reversed boundary flow.

    Dry mass is distributed uniformly by bulk volume,
    ``m_i = M_dose * V_i / sum(V)``, with the last cell corrected so the
    masses sum to the requested dose to floating-point precision. The dose
    is never inferred from the pressure-dependent porosity field.

    Boundary enforcement is face by face from the actual flux arrays (not
    aggregate diagnostics): every top face ``Fz[:, 0]`` must satisfy
    ``F >= -tol`` (forward/downward or zero) and every open bottom face
    ``Fz[open, Nz]`` must satisfy ``Fb >= -tol``; closed-rim bottom faces
    must satisfy ``|Fb| <= tol``. Rate sums ``Q_in``/``Q_out`` are
    recomputed from the (clamped) arrays. Faces in ``[-tol, 0)`` are
    clamped to exactly zero; the neglected rate is at most ``tol`` per
    face (numerical noise scale, bounded by ``FLOW_ATOL + tol_res*Q``)
    and is covered by the per-step global budget check, so no silent
    mass source is introduced. Interior radial/axial faces keep their
    signed values (counter-flow there is physical). A deliberately
    mixed-sign boundary with positive aggregate flow is rejected.
    """
    if not case.converged or case.sol is None or case.diag is None:
        raise TransportError("Cannot freeze a non-converged hydraulic case")
    tol = max(FLOW_ATOL_M3S, case.sol.balance_limit_m3s)
    grid = case.grid
    Fr = np.array(case.sol.Fr, dtype=float)
    Fz = np.array(case.sol.Fz, dtype=float)
    # Face-by-face boundary gate on the actual arrays.
    for i in range(grid.Nr):
        f_top = float(Fz[i, 0])
        if f_top < -tol:
            raise TransportError(
                "Unsupported boundary: reversed top face "
                f"i={i} F={f_top:.3e} m^3/s < -tol ({tol:.3e}). "
                "Transport requires forward-or-zero boundary flow face by face."
            )
        if f_top < 0.0:
            Fz[i, 0] = 0.0
    for i in range(grid.Nr):
        f_bot = float(Fz[i, grid.Nz])
        if grid.open_mask[i]:
            if f_bot < -tol:
                raise TransportError(
                    "Unsupported boundary: reversed outlet face "
                    f"i={i} F={f_bot:.3e} m^3/s < -tol ({tol:.3e}). "
                    "Transport requires forward-or-zero boundary flow face by face."
                )
            if f_bot < 0.0:
                Fz[i, grid.Nz] = 0.0
        else:
            if abs(f_bot) > tol:
                raise TransportError(
                    "Unsupported boundary: non-zero closed-rim face "
                    f"i={i} F={f_bot:.3e} m^3/s (|F| > tol {tol:.3e})."
                )
            Fz[i, grid.Nz] = 0.0
    # Recompute rate sums from the actual (clamped) arrays.
    q_in = float(np.sum(Fz[:, 0]))
    q_out = float(np.sum(Fz[grid.open_mask, grid.Nz])) if np.any(grid.open_mask) else 0.0
    volumes = grid.volumes
    m = M_dose * volumes / np.sum(volumes)
    m.flat[-1] += M_dose - np.sum(m)
    return HydraulicSnapshot(
        grid=grid,
        phi=np.array(case.sol.phi, dtype=float),
        m=m,
        W=np.array(case.sol.phi, dtype=float) * volumes,
        Fr=Fr,
        Fz=Fz,
        Q_in=q_in,
        Q_out=q_out,
    )


def zero_flow_snapshot(
    grid: Grid, phi_value: float = 0.35, M_dose: float = 0.018
) -> HydraulicSnapshot:
    """Synthetic no-flow snapshot for core tests (not a pressure feature).

    Zero face fluxes, uniform porosity, same dry-mass rule as production.
    """
    volumes = grid.volumes
    m = M_dose * volumes / np.sum(volumes)
    m.flat[-1] += M_dose - np.sum(m)
    phi = np.full((grid.Nr, grid.Nz), phi_value)
    return HydraulicSnapshot(
        grid=grid,
        phi=phi,
        m=m,
        W=phi * volumes,
        Fr=np.zeros((grid.Nr + 1, grid.Nz)),
        Fz=np.zeros((grid.Nr, grid.Nz + 1)),
        Q_in=0.0,
        Q_out=0.0,
    )


@dataclass
class TransportAssembly:
    """Discrete transport operator: net outward solute flux = L*C - b_in."""

    L: sparse.csr_matrix  # (N, N) in m^3/s, idx = j*Nr + i
    top_inflow: np.ndarray  # (Nr, Nz) inward top-face volume rates [m^3/s]
    outflow: np.ndarray  # (Nr, Nz) positive outward volume rates [m^3/s]


def assemble_transport(
    grid: Grid, phi: np.ndarray, D: float, Fr: np.ndarray, Fz: np.ndarray
) -> TransportAssembly:
    """Assemble the solute transport operator from frozen volume fluxes."""
    Nr, Nz = grid.Nr, grid.Nz
    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []
    top_inflow = np.zeros((Nr, Nz))
    outflow = np.zeros((Nr, Nz))

    def add(r: int, c: int, v: float) -> None:
        if v != 0.0:
            rows.append(r)
            cols.append(c)
            data.append(v)

    def interior(r_l: int, c_l: tuple[int, int], r_r: int, c_r: tuple[int, int],
                 F: float, H: float) -> None:
        i_l, j_l = c_l
        i_r, j_r = c_r
        if F >= 0.0:
            # J = F*C_L + H*(C_L - C_R); +J outward from L, -J from R.
            add(r_l, r_l, F + H)
            add(r_l, r_r, -H)
            add(r_r, r_l, -(F + H))
            add(r_r, r_r, H)
            outflow[i_l, j_l] += F
        else:
            # J = F*C_R + H*(C_L - C_R).
            add(r_l, r_l, H)
            add(r_l, r_r, F - H)
            add(r_r, r_l, -H)
            add(r_r, r_r, H - F)
            outflow[i_r, j_r] += -F

    for j in range(Nz):
        for i_face in range(1, Nr):
            i_l, i_r = i_face - 1, i_face
            F = float(Fr[i_face, j])
            H = D * _hmean(phi[i_l, j], phi[i_r, j]) * grid.rad_areas[i_face] / grid.dr
            interior(cell_index(Nr, i_l, j), (i_l, j),
                     cell_index(Nr, i_r, j), (i_r, j), F, H)
    for i in range(Nr):
        for j_face in range(1, Nz):
            j_u, j_d = j_face - 1, j_face
            F = float(Fz[i, j_face])
            H = D * _hmean(phi[i, j_u], phi[i, j_d]) * grid.col_areas[i] / grid.dz
            interior(cell_index(Nr, i, j_u), (i, j_u),
                     cell_index(Nr, i, j_d), (i, j_d), F, H)
    for i in range(Nr):
        r = cell_index(Nr, i, 0)
        F = float(Fz[i, 0])
        if F < 0.0:
            raise TransportError(
                "Unsupported boundary: reversed top flow "
                f"i={i} F={F:.3e} m^3/s. Transport requires forward-or-zero "
                "boundary flow; reversed boundaries are rejected, not accommodated."
            )
        if F > 0.0:
            top_inflow[i, 0] = F  # inward total flux F*C_in(t) via b_in
        # F == 0 contributes nothing.
        if grid.open_mask[i]:
            r_bot = cell_index(Nr, i, Nz - 1)
            Fb = float(Fz[i, Nz])
            if Fb < 0.0:
                raise TransportError(
                    "Unsupported boundary: reversed outlet flow "
                    f"i={i} F={Fb:.3e} m^3/s. Transport requires forward-or-zero "
                    "boundary flow; reversed boundaries are rejected, not accommodated."
                )
            add(r_bot, r_bot, Fb)  # outward advection (forward flow: Fb >= 0)
            if Fb > 0.0:
                outflow[i, Nz - 1] += Fb
        else:
            Fb = float(Fz[i, Nz])
            if Fb != 0.0:
                raise TransportError(
                    "Unsupported boundary: non-zero closed-rim flux "
                    f"i={i} F={Fb:.3e} m^3/s. Rim faces must carry zero flow."
                )
    L = sparse.coo_matrix(
        (data, (rows, cols)), shape=(Nr * Nz, Nr * Nz)
    ).tocsr()
    return TransportAssembly(L=L, top_inflow=top_inflow, outflow=outflow)


def transport_step(
    C_old: np.ndarray,
    B_old: np.ndarray,
    asm: TransportAssembly,
    W: np.ndarray,
    m: np.ndarray,
    K: float,
    lam: float,
    dt: float,
    C_in: float,
) -> tuple[np.ndarray, np.ndarray, float]:
    """One backward-Euler step with the solid pool eliminated algebraically.

    ``beta = dt*lam/(1+dt*lam)``; solves
    ``[diag(W + beta*K*m) + dt*L] C_new = W*C_old + beta*B_old + dt*b_in``
    with ``b_in = top_inflow*C_in``, then
    ``B_new = (B_old + dt*lam*K*m*C_new)/(1+dt*lam)``.
    Returns ``(C_new, B_new, lin_rel)`` with the relative inf-norm linear
    residual. No second release source is added after the solve.
    """
    if not dt > 0.0:
        raise TransportError(f"Step requires dt > 0, got {dt!r}")
    beta = dt * lam / (1.0 + dt * lam)
    Nr, Nz = W.shape
    A = sparse.diags((W + beta * K * m).T.ravel()) + dt * asm.L
    b_in = asm.top_inflow * C_in
    rhs = (W * C_old + beta * B_old + dt * b_in).T.ravel()
    A = A.tocsr()
    C_flat = np.asarray(spsolve(A, rhs), dtype=float)
    if np.any(~np.isfinite(C_flat)):
        raise TransportError("Transport linear solve returned non-finite C")
    rhs_scale = max(float(np.max(np.abs(rhs))), 1.0e-300)
    lin_rel = float(np.max(np.abs(A @ C_flat - rhs)) / rhs_scale)
    C_new = C_flat.reshape(Nz, Nr).T
    B_new = (B_old + dt * lam * K * m * C_new) / (1.0 + dt * lam)
    return C_new, B_new, lin_rel


def advective_dt_limit(
    W: np.ndarray, outflow: np.ndarray, cfl: float
) -> float | None:
    """Conservative accuracy limit ``cfl*min(W/outflow)`` over outflow cells.

    Returns None when the limit is disabled (cfl <= 0) or no cell has
    positive outflow. Backward Euler is stable beyond it; the limit guards
    accuracy, and resolution is assessed by refinement, not by stability.
    """
    if cfl <= 0.0:
        return None
    mask = outflow > 0.0
    if not np.any(mask):
        return None
    return float(cfl * np.min(W[mask] / outflow[mask]))


@dataclass
class StepRecord:
    t: float
    dt: float
    C_in: float
    M_cup: float
    V_cup: float
    outlet_C: float  # flux-weighted instantaneous outlet conc [kg/m^3]
    sum_B: float
    sum_WC: float
    cumulative_inlet: float
    ledger_signed: float
    lin_rel: float
    min_C: float
    min_B: float


@dataclass
class TransportResult:
    C: np.ndarray  # final liquid concentration [kg/m^3]
    B: np.ndarray  # final solid inventory [kg]
    m: np.ndarray  # dry mass per cell [kg]
    W: np.ndarray  # storage volume per cell [m^3]
    e_local: np.ndarray  # local extraction fraction (B0-B)/m [-]
    e_mean: float  # dry-mass-weighted mean [-]
    e_var: float  # dry-mass-weighted variance [-]
    steps: list[StepRecord]
    n_steps: int
    dt_requested: float
    dt_advective_limit: float | None
    budget_max_abs: float
    budget_tol: float
    lin_rel_max: float
    M_cup: float
    V_cup: float
    TDS_pct: float | None
    TDS_status: str
    EY_pct: float | None  # None in tracer mode: no coffee EY for injected tracer
    EY_status: str
    M_dose: float
    f_sol: float


@dataclass
class TransportCase:
    converged: bool
    result: TransportResult | None = None
    n_steps: int = 0
    error: str = ""
    error_type: str = ""
    Q_out: float = float("nan")
    tau_s: float | None = None  # None when Q_out == 0 (no residence time)


def cup_metrics(
    M_cup: float, V_cup: float, rho_solution: float, M_dose: float
) -> tuple[float | None, str, float]:
    """TDS percent (None when the cup is empty), status, and EY percent.

    Raises TransportError if solute mass exceeds solution mass instead of
    clipping TDS.
    """
    M_solution = rho_solution * V_cup
    if V_cup <= 0.0:
        return None, "undefined_empty_cup", 0.0
    if M_cup > M_solution:
        raise TransportError(
            f"Invalid cup state: solute mass {M_cup:.6e} kg exceeds solution "
            f"mass {M_solution:.6e} kg; TDS is not reported."
        )
    return 100.0 * M_cup / M_solution, "ok", 100.0 * M_cup / M_dose


def run_transport(hparams: SimParams, tparams: TransportParams) -> TransportCase:
    """Run passive/exchange transport on one frozen hydraulic snapshot.

    Time stepping hits every inlet-schedule breakpoint exactly (no pulse
    is skipped silently) and always advances to the requested ``t_end`` on
    success. Representability rule: ``t_end`` and every breakpoint must be
    finite with ``0 < t <= t_end`` reachable; a step with ``t+dt == t``
    (nonadvancing) fails explicitly via ``TransportError`` instead of
    stalling or skipping. Small but representable horizons (e.g.
    ``t_end=1e-13``) take explicit small steps.
    """
    validate(tparams)
    case = run_case(hparams)
    if not case.converged:
        return TransportCase(converged=False, error=f"Hydraulics: {case.error}",
                             error_type="ConvergenceError")
    try:
        snap = freeze_snapshot(case, tparams.M_dose)
        grid = snap.grid
        asm = assemble_transport(grid, snap.phi, tparams.D, snap.Fr, snap.Fz)
    except TransportError as exc:
        q = float("nan")
        try:
            q = float(case.diag.Q_out_m3s) if case.diag is not None else float("nan")
        except Exception:
            pass
        return TransportCase(converged=False, error=str(exc),
                             error_type=type(exc).__name__, Q_out=q)

    if tparams.mode == "extraction":
        B = tparams.f_sol * snap.m
        C = np.zeros_like(snap.W)
    else:
        B = np.zeros_like(snap.W)
        C = np.full_like(snap.W, tparams.C_init)
    B0 = np.array(B)

    initial_total = float(np.sum(B) + np.sum(snap.W * C))
    dt_limit = advective_dt_limit(snap.W, asm.outflow, tparams.cfl_advective)
    dt_req = tparams.dt if dt_limit is None else min(tparams.dt, dt_limit)
    if not (math.isfinite(dt_req) and dt_req > 0.0):
        return TransportCase(converged=False,
                             error=f"Unsupported time step: dt_req={dt_req!r} "
                             "is not finite-positive.",
                             error_type="TransportError", Q_out=snap.Q_out)
    t_end = float(tparams.t_end)
    if not (math.isfinite(t_end) and t_end > 0.0):
        # validate() already rejects this; explicit guard before solving.
        return TransportCase(converged=False,
                             error=f"Unsupported horizon: t_end={t_end!r}.",
                             error_type="TransportError", Q_out=snap.Q_out)
    targets = sorted(
        {float(x) for x in schedule_breakpoints(tparams)
         if math.isfinite(float(x)) and 0.0 < float(x) < t_end}
        | {t_end}
    )
    # Explicit representability pre-check: the first step must advance.
    _first_nxt = min(targets)
    _first_dt = dt_req if dt_req < _first_nxt else _first_nxt
    if not (math.isfinite(_first_dt) and _first_dt > 0.0 and (0.0 + _first_dt) > 0.0):
        return TransportCase(
            converged=False,
            error=f"Unsupported horizon: t_end={t_end!r} with dt_req={dt_req!r} "
            "does not advance from t=0.",
            error_type="TransportError", Q_out=snap.Q_out)
    M_cup = 0.0
    V_cup = 0.0
    cumulative_inlet = 0.0
    steps: list[StepRecord] = []
    budget_max_abs = 0.0
    lin_rel_max = 0.0
    t = 0.0
    try:
        while t < t_end:
            candidates = [x for x in targets if x > t]
            if not candidates:
                raise TransportError(
                    f"No forward schedule target from t={t:.6e} s "
                    f"toward t_end={t_end:.6e} s."
                )
            nxt = min(candidates)
            if not (math.isfinite(nxt) and nxt > t and nxt <= t_end):
                raise TransportError(
                    f"Invalid schedule target nxt={nxt!r} from t={t:.6e} s."
                )
            gap = nxt - t
            if not (math.isfinite(gap) and gap > 0.0):
                raise TransportError(
                    f"Nonadvancing schedule gap nxt-t={gap!r} at t={t:.6e} s."
                )
            dt_step = dt_req if dt_req < gap else gap
            if not (math.isfinite(dt_step) and dt_step > 0.0):
                raise TransportError(
                    f"Non-positive step dt={dt_step!r} at t={t:.6e} s."
                )
            # Exact boundary hit avoids floating accumulation error.
            if dt_step == gap:
                t_next = nxt
            else:
                t_next = t + dt_step
            if not (math.isfinite(t_next) and t_next > t):
                raise TransportError(
                    f"Nonadvancing t+dt: t={t:.6e} s dt={dt_step:.6e} s "
                    "does not advance; failing explicitly."
                )
            C_in = inlet_concentration(tparams, t)
            C, B, lin_rel = transport_step(
                C, B, asm, snap.W, snap.m, tparams.K,
                tparams.exchange_rate, dt_step, C_in,
            )
            lin_rel_max = max(lin_rel_max, lin_rel)
            if lin_rel > LIN_REL_TOL:
                raise TransportError(
                    f"Linear residual {lin_rel:.3e} exceeds {LIN_REL_TOL:.1e} "
                    f"at t={t:.6f} s"
                )
            c_scale = max(float(np.max(C)), C_in, 1.0e-30)
            if float(np.min(C)) < -POS_REL_TOL * c_scale:
                raise TransportError(
                    f"Negative concentration {float(np.min(C)):.3e} kg/m^3 "
                    f"beyond tolerance at t={t:.6f} s; no clipping applied"
                )
            b_scale = max(float(np.max(B)), 1.0e-30)
            if float(np.min(B)) < -POS_REL_TOL * b_scale:
                raise TransportError(
                    f"Negative solid inventory {float(np.min(B)):.3e} kg "
                    f"beyond tolerance at t={t:.6f} s; no clipping applied"
                )
            # Cup collection uses the same numerical outlet flux as transport.
            dM = dt_step * float(
                np.sum(snap.Fz[snap.grid.open_mask, grid.Nz] * C[snap.grid.open_mask, grid.Nz - 1])
            )
            dV = dt_step * snap.Q_out
            dIn = dt_step * float(np.sum(asm.top_inflow)) * C_in
            M_cup += dM
            V_cup += dV
            cumulative_inlet += dIn
            total = float(np.sum(B) + np.sum(snap.W * C)) + M_cup
            expected = initial_total + cumulative_inlet
            signed = total - expected
            budget_max_abs = max(budget_max_abs, abs(signed))
            budget_tol = BUDGET_ATOL_KG + BUDGET_RTOL * max(
                initial_total, cumulative_inlet
            )
            if abs(signed) > budget_tol:
                raise TransportError(
                    f"Global budget error {signed:.3e} kg exceeds "
                    f"{budget_tol:.3e} kg at t={t + dt_step:.6f} s"
                )
            if snap.Q_out > 0.0:
                outlet_C = float(
                    np.sum(snap.Fz[snap.grid.open_mask, grid.Nz]
                           * C[snap.grid.open_mask, grid.Nz - 1]) / snap.Q_out
                )
            else:
                outlet_C = float("nan")
            t = t_next
            steps.append(StepRecord(
                t=t, dt=dt_step, C_in=C_in, M_cup=M_cup, V_cup=V_cup,
                outlet_C=outlet_C, sum_B=float(np.sum(B)),
                sum_WC=float(np.sum(snap.W * C)),
                cumulative_inlet=cumulative_inlet, ledger_signed=signed,
                lin_rel=lin_rel, min_C=float(np.min(C)),
                min_B=float(np.min(B)),
            ))
        TDS_pct, TDS_status, EY_pct = cup_metrics(
            M_cup, V_cup, tparams.rho_solution, tparams.M_dose
        )
        if tparams.mode == "tracer":
            # Injected tracer is external solute: no coffee EY is reported.
            EY_pct, EY_status = None, "not_coffee_tracer"
        else:
            EY_status = "coffee"
    except TransportError as exc:
        return TransportCase(
            converged=False, n_steps=len(steps), error=str(exc),
            error_type=type(exc).__name__,
            Q_out=snap.Q_out,
        )
        # Successful positive horizons must reach the declared endpoint.
    if not steps or steps[-1].t != t_end:
        return TransportCase(
            converged=False, n_steps=len(steps),
            error=f"Did not reach t_end={t_end!r}: final t="
            f"{steps[-1].t if steps else 0.0!r}.",
            error_type="TransportError", Q_out=snap.Q_out)

    e_local = (B0 - B) / snap.m
    e_mean = float(np.sum(snap.m * e_local) / tparams.M_dose)
    e_var = float(np.sum(snap.m * (e_local - e_mean) ** 2) / tparams.M_dose)
    budget_tol = BUDGET_ATOL_KG + BUDGET_RTOL * max(initial_total, cumulative_inlet)
    tau_s = float(np.sum(snap.W) / snap.Q_out) if snap.Q_out > 0.0 else None
    return TransportCase(
        converged=True,
        result=TransportResult(
            C=C, B=B, m=snap.m, W=snap.W, e_local=e_local,
            e_mean=e_mean, e_var=e_var, steps=steps, n_steps=len(steps),
            dt_requested=tparams.dt, dt_advective_limit=dt_limit,
            budget_max_abs=budget_max_abs, budget_tol=budget_tol,
            lin_rel_max=lin_rel_max, M_cup=M_cup, V_cup=V_cup,
            TDS_pct=TDS_pct, TDS_status=TDS_status, EY_pct=EY_pct,
            EY_status=EY_status,
            M_dose=tparams.M_dose, f_sol=tparams.f_sol,
        ),
        n_steps=len(steps),
        Q_out=snap.Q_out,
        tau_s=tau_s,
    )


ASSUMPTIONS = [
    "One frozen saturated isothermal hydraulic snapshot per run; grid, porosity, permeability, and face fluxes fixed.",
    "Time zero is an already wetted bed; no wetting, gas, deformation, fines, heat, or particle diffusion.",
    "Lumped reversible solid/liquid exchange; passive storage (dissolution does not alter geometry or properties).",
    "Backward Euler, first-order upwind advection, centered dispersion; no TVD limiter.",
    "All transport/exchange coefficients are illustrative assumptions; no experimental calibration is supplied.",
]
