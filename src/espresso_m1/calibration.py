"""M4 calibration harness: empirical flow fit with uncertainty, held-out verdicts.

Fits a single identifiable quantity — the constant-flow slope ``q`` over a
declared already-wet interval — on TRAIN whole shots only, reports a 95 %
t-interval (mean) and 95 % prediction interval (new shot), and verdicts each
HELD-OUT whole shot against the prediction interval. Optionally compares the
fitted ``q`` against a simulator constant ``Q`` (empirical-vs-simulated bias;
permeability is never varied to fit). Failures to predict are disclosed, not
hidden.

Status: harness verified on synthetic fixtures only. The first run involving
``source: measured`` data requires scientific review; see ``real_data_note``.
Cup TDS alone still cannot determine f_sol, K, lambda, D, k0.
"""

from __future__ import annotations

import math
import statistics

from .shots import (CALIBRATION_STATUS, LIMITATIONS_NOTE, Shot,
                    ShotError, fit_constant_flow_train)

T_LABEL = "95 % Student-t"


def _slope(shot: Shot, t0_s: float) -> float:
    ts = [p.time_s for p in shot.samples if p.time_s >= t0_s]
    ms = [p.beverage_mass_g for p, t in zip(shot.samples, [p.time_s for p in shot.samples]) if t >= t0_s]
    if len(ts) < 2 or (ts[-1] - ts[0]) <= 0:
        raise ShotError(f"shot {shot.metadata.shot_id}: no usable interval at/after t0={t0_s!r}")
    return (ms[-1] - ms[0]) / (ts[-1] - ts[0])


def calibrate(train: list[Shot], heldout: list[Shot], t0_s: float,
              sim_q_g_s: float | None = None) -> dict:
    """Fit q on train shots; verdict held-out shots; compare with simulator Q."""
    if not train:
        raise ShotError("calibrate needs >= 1 training shot")
    if not math.isfinite(t0_s):
        raise ShotError(f"t0 must be finite, got {t0_s!r}")
    if sim_q_g_s is not None and not math.isfinite(sim_q_g_s):
        raise ShotError(f"sim-q must be finite, got {sim_q_g_s!r}")
    fit = fit_constant_flow_train(train, t0_s)
    slopes: list[float] = fit["train_slopes"]
    n = len(slopes)
    q = fit["q_g_s"]
    ci95 = pi95 = None
    t95 = None
    if n >= 2:
        from scipy.stats import t as _t
        s = statistics.stdev(slopes)
        t95 = float(_t.ppf(0.975, n - 1))
        se = s / math.sqrt(n)
        ci95 = [q - t95 * se, q + t95 * se]
        half = t95 * s * math.sqrt(1.0 + 1.0 / n)
        pi95 = [q - half, q + half]
    verdicts: list[dict] = []
    for s in heldout:
        h = _slope(s, t0_s)
        if pi95 is None:
            verdicts.append({"shot_id": s.metadata.shot_id, "slope_g_s": h,
                             "verdict": "unassessed_single_train",
                             "detail": "n_train < 2: no prediction interval; fit reported without verdict"})
        else:
            ok = pi95[0] <= h <= pi95[1]
            verdicts.append({"shot_id": s.metadata.shot_id, "slope_g_s": h,
                             "verdict": "pass" if ok else "FAIL",
                             "interval_g_s": list(pi95),
                             "detail": ("held-out slope inside 95 % prediction interval"
                                        if ok else "held-out slope OUTSIDE 95 % prediction interval")})
    assessed = [v for v in verdicts if v["verdict"] in ("pass", "FAIL")]
    if not assessed:
        status = "indeterminate_single_train" if verdicts else "no_heldout"
        passed = False if verdicts else True
    else:
        status = "pass" if all(v["verdict"] == "pass" for v in assessed) else "FAIL"
        passed = status == "pass"
    sim_cmp: dict | None = None
    if sim_q_g_s is not None:
        bias = q - float(sim_q_g_s)
        inside = (ci95[0] <= float(sim_q_g_s) <= ci95[1]) if ci95 is not None else None
        sim_cmp = {"sim_q_g_s": float(sim_q_g_s), "empirical_q_g_s": q,
                   "bias_empirical_minus_sim_g_s": bias,
                   "rel_bias": (bias / float(sim_q_g_s) if float(sim_q_g_s) != 0 else None),
                   "sim_inside_95ci": inside,
                   "note": ("Positive bias = measured flow faster than simulator constant Q. "
                            "Simulator Q is overlaid as-is; k is NOT refitted.")}
    sources = {s.metadata.source for s in list(train) + list(heldout)}
    return {"fit_q_g_s": q, "t0_s": t0_s, "n_train": n,
            "train_slopes_g_s": slopes,
            "t95": t95, "interval_label": T_LABEL if n >= 2 else None,
            "ci95_mean_g_s": ci95, "pi95_new_shot_g_s": pi95,
            "verdicts": verdicts, "status": status, "passed": passed,
            "simulator_comparison": sim_cmp,
            "all_synthetic": sources <= {"synthetic"},
            "real_data_note": ("Harness verified on synthetic fixtures only. "
                               "First run with source=measured data requires scientific review; "
                               "a pass here is code verification, not coffee validation."
                               if "measured" in sources else
                               "Synthetic-only run: importer/harness correctness, no coffee claims."),
            "limitations": LIMITATIONS_NOTE + " Fit uses whole training shots only.",
            "calibration_status": CALIBRATION_STATUS}
