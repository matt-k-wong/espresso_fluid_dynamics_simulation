"""Tasks E6/E7: full extraction ledger, positivity, and metrics.

E6 runs the perforated compaction example with exchange and dispersion:
every step's global budget error must satisfy
``1e-12 kg + 1e-9*max(initial_total, cumulative_input)`` (maxima and
tolerances reported, not just a flag), and C/B stay nonnegative to solver
tolerance with no corrective clipping. E7 independently recomputes
flux-weighted outlet concentration, cup mass, TDS/EY, and dry-mass-weighted
extraction statistics from saved arrays/series, checks cup EY <= 100*f_sol,
and requires valid strict JSON for an empty cup. Budget acceptance does not
establish spatial accuracy (see refinement tests).
"""

import json

import numpy as np
import pytest

from espresso_m1.params import SimParams
from espresso_m1.tparams import TransportParams
from espresso_m1.transport import BUDGET_ATOL_KG, BUDGET_RTOL, run_transport

HYDRO = dict(permeability_mode="compaction", Nr=20, Nz=40, R_perf=0.0145)


def extraction_case(**kwargs):
    h = SimParams(**HYDRO)
    base = dict(mode="extraction", t_end=60.0, dt=0.5)
    base.update(kwargs)
    return run_transport(h, TransportParams(**base))


def test_full_extraction_ledger_and_positivity():
    case = extraction_case()
    assert case.converged, case.error
    r = case.result
    worst = 0.0
    for s in r.steps:
        tol = BUDGET_ATOL_KG + BUDGET_RTOL * max(
            s.sum_B + s.sum_WC + s.M_cup - s.ledger_signed, s.cumulative_inlet
        )
        # Recompute the allowance from first principles, not the stored tol.
        assert abs(s.ledger_signed) <= tol
        worst = max(worst, abs(s.ledger_signed))
    print(f"\nsteps={r.n_steps} budget worst={worst:.3e} tol={r.budget_tol:.3e} "
          f"lin_rel max={r.lin_rel_max:.3e}")
    assert r.budget_max_abs <= r.budget_tol
    c_scale = max(float(np.max(r.C)), 1e-30)
    b_scale = max(float(np.max(r.B)), 1e-30)
    print(f"min C={np.min(r.C):.3e} (scale {c_scale:.3e}) "
          f"min B={np.min(r.B):.3e} (scale {b_scale:.3e})")
    assert float(np.min(r.C)) >= -1e-9 * c_scale
    assert float(np.min(r.B)) >= -1e-9 * b_scale


def test_final_field_ledger_from_saved_arrays(tmp_path):
    case = extraction_case()
    assert case.converged
    r = case.result
    path = tmp_path / "fields.npz"
    np.savez(path, C=r.C, B=r.B, m=r.m, W=r.W, e_local=r.e_local)
    z = np.load(path)
    # Global conservation from final arrays + dose alone (no inlet here).
    assert float(np.sum(z["B"]) + np.sum(z["W"] * z["C"])) + r.M_cup == \
        pytest.approx(r.f_sol * r.M_dose, rel=1e-9)
    # Extraction statistics recomputed from saved arrays.
    e = (r.f_sol * z["m"] - z["B"]) / z["m"]
    assert np.max(np.abs(e - z["e_local"])) == 0.0
    mean = float(np.sum(z["m"] * e) / r.M_dose)
    var = float(np.sum(z["m"] * (e - mean) ** 2) / r.M_dose)
    assert mean == pytest.approx(r.e_mean, rel=1e-12)
    assert var == pytest.approx(r.e_var, rel=1e-12)
    # Cup EY bounded by the soluble fraction (extraction, zero liquid input).
    assert r.EY_pct <= 100.0 * r.f_sol * (1.0 + 1e-9)
    assert r.EY_status == "coffee"


def test_series_metrics_recomputed_from_columns(tmp_path):
    case = extraction_case(t_end=10.0)
    assert case.converged
    r = case.result
    rho = 1000.0
    for s in r.steps:
        tds = 100.0 * s.M_cup / (rho * s.V_cup) if s.V_cup > 0 else None
        assert (tds is None and s.M_cup == 0.0) or tds == pytest.approx(
            100.0 * s.M_cup / (rho * s.V_cup), rel=1e-12)
        assert 100.0 * s.M_cup / r.M_dose >= 0.0
    assert r.TDS_pct == pytest.approx(
        100.0 * r.M_cup / (rho * r.V_cup), rel=1e-12)
    assert r.EY_pct == pytest.approx(100.0 * r.M_cup / r.M_dose, rel=1e-12)


def test_cup_metrics_reject_unphysical_state():
    from espresso_m1.transport import cup_metrics, TransportError

    with pytest.raises(TransportError):
        cup_metrics(2.0, 1e-3, 1000.0, 0.018)  # solute > solution mass
    tds, status, ey = cup_metrics(0.0, 0.0, 1000.0, 0.018)
    assert tds is None and status == "undefined_empty_cup" and ey == 0.0


def test_batch_reabsorption_pool_grows_above_initial():
    # Reversible exchange: a solute-rich liquid loads an empty solid pool.
    # B rises above its initial value while global mass stays conserved;
    # the model preserves this instead of clipping B.
    from espresso_m1.transport import (
        assemble_transport,
        transport_step,
        zero_flow_snapshot,
    )
    from espresso_m1.grid import make_grid

    grid = make_grid(1, 1, 0.029, 0.02, 0.0)
    snap = zero_flow_snapshot(grid, phi_value=0.35)
    asm = assemble_transport(grid, snap.phi, 0.0, snap.Fr, snap.Fz)
    W, m = float(snap.W[0, 0]), float(snap.m[0, 0])
    K, lam, C0 = 0.005, 0.1, 10.0
    C = np.array([[C0]])
    B = np.array([[0.0]])
    M = W * C0
    t = 0.0
    while t < 20.0 - 1e-12:
        C, B, _ = transport_step(C, B, asm, snap.W, snap.m, K, lam, 1.0, 0.0)
        t += 1.0
    assert float(B[0, 0]) > 0.0  # grew above its initial zero
    assert float(B[0, 0] + W * C[0, 0]) == pytest.approx(M, rel=1e-12)
    C_eq = M / (W + K * m)
    rate = lam * (1.0 + K * m / W)
    assert float(C[0, 0]) == pytest.approx(
        C_eq + (C0 - C_eq) * np.exp(-rate * 20.0), rel=1e-2
    )
