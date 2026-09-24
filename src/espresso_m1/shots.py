"""Versioned GS3 shot observations: schema, importer, metrics, comparison.

Schema version ``gs3_shot_v1``. No measured GS3 data is invented here; public
fixtures are labeled ``synthetic``. This module never connects to or controls
the machine.

Shot metadata (JSON) covers: shot ID, source (measured/synthetic), machine
profile ID, AV/reservoir feed, bean batch/roast date if known, grinder name
and ordinal setting (not particle diameter), dose_g, basket ID, beverage
target_g, user temperature setting with measurement location, pre-brew
settings only if confirmed, time-zero convention, notes, and optional cup
TDS with method/uncertainty. No exact machine serial is stored in public
examples.

Time series (inline JSON ``samples`` or external CSV) covers: time_s,
beverage_mass_g, optional pressure_Pa + sensor location, optional
temperature_K + sensor location. Missing data stays missing (null/empty).
``time_zero`` records whether time zero is pump activation, first liquid, or
the modeled already-wet state.

Validation rejects nonfinite values, decreasing times, duplicate timestamps
unless an explicit aggregation rule is selected, negative dose, and
inconsistent units. Raw mass samples are retained; small negative mass
increments from balance noise are flagged, and flow is derived by a
documented optional smoothing method (window + endpoint treatment reported).
Raw data is never rewritten to force monotonicity, and a noisy derivative is
reported as derived (smoothed) flow, never as directly measured flow.

Cup EY is only computed when TDS (or another actual solute measurement) is
present; it is never inferred from mass alone.

The current simulator predicts constant Q in a transport run. Comparisons
use an explicit constant-flow baseline over a declared already-wet interval
and disclose this limitation; permeability is never secretly varied to fit
a trace. Effective conductance versus k + pressure-loss decomposition, and
f_sol/K/lambda/D/k0 from TDS alone, are documented as non-identifiable.
Calibration against real coffee is disabled until real data is supplied
(see ``calibration_status``).
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass, field

SHOT_SCHEMA_VERSION = "gs3_shot_v1"

SOURCES = ("measured", "synthetic")
TIME_ZERO_KINDS = ("pump_activation", "first_liquid", "already_wet")
AGGREGATIONS = ("mean", "first", "last", "median")

CALIBRATION_STATUS = (
    "Calibration against real coffee is DISABLED until real measured shot "
    "data is supplied. Synthetic fixtures validate importer/metric "
    "correctness only and must never be labeled measured."
)

LIMITATIONS_NOTE = (
    "Simulator limitation: a transport run predicts constant Q (frozen "
    "hydraulics, already-wet bed). Any comparison overlays this constant-flow "
    "baseline and reports residuals; it does not vary permeability to fit the "
    "trace. Without puck-pressure measurement, effective conductance may be "
    "estimated but k and pressure losses cannot both be claimed identified. "
    "Cup TDS alone does not uniquely determine f_sol, K, lambda, D and k0."
)


class ShotError(ValueError):
    """Invalid shot observation (schema/validation failure)."""


@dataclass
class ShotSample:
    time_s: float
    beverage_mass_g: float
    pressure_Pa: float | None = None
    pressure_location: str | None = None
    temperature_K: float | None = None
    temperature_location: str | None = None


@dataclass
class ShotMetadata:
    shot_id: str
    source: str
    machine_profile_id: str
    feed: str
    dose_g: float
    beverage_target_g: float | None = None
    basket_id: str | None = None
    bean_batch: str | None = None
    roast_date: str | None = None
    grinder: str | None = None
    grinder_setting: object = None
    grinder_setting_is_ordinal: bool = True
    temp_setting_C: float | None = None
    temp_location: str | None = None
    prebrew: dict | None = None
    time_zero: str = "pump_activation"
    notes: str = ""
    tds_pct: float | None = None
    tds_method: str | None = None
    tds_uncertainty_pct: float | None = None


@dataclass
class Shot:
    metadata: ShotMetadata
    samples: list[ShotSample] = field(default_factory=list)
    raw_samples_retained: bool = True


def _is_finite_number(x) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(float(x))


def _require_finite(name: str, value, errors: list[str]) -> float | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        errors.append(f"{name} must be a number or null, got {value!r}")
        return None
    if not math.isfinite(float(value)):
        errors.append(f"{name} must be finite, got {value!r}")
        return None
    return float(value)


def load_series_csv(path: str) -> list[dict]:
    """Load a timeseries CSV with header time_s,beverage_mass_g,... .

    Empty fields stay missing (None). Non-numeric non-empty fields raise
    ShotError. Header must contain time_s and beverage_mass_g.
    """
    with open(path, "r", newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            raise ShotError(f"{path}: missing CSV header")
        needed = {"time_s", "beverage_mass_g"}
        missing = needed - set(reader.fieldnames)
        if missing:
            raise ShotError(f"{path}: missing columns {sorted(missing)}; "
                            f"header={reader.fieldnames}")
        rows: list[dict] = []
        for lineno, row in enumerate(reader, start=2):
            out: dict = {}
            for key in ("time_s", "beverage_mass_g", "pressure_Pa", "temperature_K"):
                raw = row.get(key, "")
                if raw is None or (isinstance(raw, str) and raw.strip() == ""):
                    out[key] = None
                else:
                    try:
                        out[key] = float(str(raw).strip())
                    except ValueError:
                        raise ShotError(f"{path}:{lineno}: {key}={raw!r} is not numeric")
            for key in ("pressure_location", "temperature_location"):
                raw = row.get(key, "")
                if raw is None or (isinstance(raw, str) and raw.strip() == ""):
                    out[key] = None
                else:
                    out[key] = str(raw).strip()
            out["_lineno"] = lineno
            rows.append(out)
    return rows


def _aggregate_duplicates(rows: list[dict], rule: str) -> list[dict]:
    """Aggregate rows sharing identical time_s per rule; locations use first non-null."""
    from collections import defaultdict
    groups: dict[float, list[dict]] = defaultdict(list)
    for r in rows:
        groups[float(r["time_s"])].append(r)
    out: list[dict] = []
    for t in sorted(groups):
        g = groups[t]
        if len(g) == 1:
            out.append(g[0])
            continue
        vals: dict = {"time_s": t}
        for key in ("beverage_mass_g", "pressure_Pa", "temperature_K"):
            nums = [r[key] for r in g if r[key] is not None]
            if not nums:
                vals[key] = None
            elif rule == "first":
                vals[key] = next(r[key] for r in g if r[key] is not None)
            elif rule == "last":
                vals[key] = next(r[key] for r in reversed(g) if r[key] is not None)
            elif rule == "median":
                s = sorted(nums)
                n = len(s)
                vals[key] = s[n // 2] if n % 2 == 1 else 0.5 * (s[n // 2 - 1] + s[n // 2])
            else:  # mean
                vals[key] = sum(nums) / len(nums)
        for key in ("pressure_location", "temperature_location"):
            vals[key] = next((r[key] for r in g if r[key]), None)
        vals["_lineno"] = g[0].get("_lineno")
        vals["_aggregated_n"] = len(g)
        out.append(vals)
    return out


def parse_shot(data: dict, series_rows: list[dict] | None = None,
               aggregate: str | None = None) -> Shot:
    """Validate a shot JSON dict (plus optional external CSV rows) into a Shot.

    ``aggregate`` selects duplicate-timestamp handling (mean/first/last/median);
    without it duplicates are rejected. Raises ShotError on any violation.
    Raw samples are retained unchanged (no monotonicity forcing).
    """
    errors: list[str] = []
    if not isinstance(data, dict):
        raise ShotError("Shot JSON must be an object")
    if data.get("schema_version") != SHOT_SCHEMA_VERSION:
        raise ShotError(f"schema_version must be {SHOT_SCHEMA_VERSION!r}, "
                        f"got {data.get('schema_version')!r}")
    shot_id = data.get("shot_id")
    if not isinstance(shot_id, str) or not shot_id.strip():
        errors.append("shot_id must be a non-empty string")
    source = data.get("source")
    if source not in SOURCES:
        errors.append(f"source must be one of {SOURCES}, got {source!r}")
    machine = data.get("machine") or {}
    profile_id = machine.get("profile_id", data.get("machine_profile_id", ""))
    if not isinstance(profile_id, str) or not profile_id.strip():
        errors.append("machine.profile_id must be a non-empty string")
    feed = machine.get("feed", data.get("feed", ""))
    if not isinstance(feed, str) or not feed.strip():
        errors.append("machine.feed must be a non-empty string (e.g. 'AV reservoir')")
    dose_g = _require_finite("dose_g", data.get("dose_g"), errors)
    if dose_g is not None and not dose_g > 0:
        errors.append(f"dose_g must be positive, got {dose_g!r}")
    bev_target = _require_finite("beverage_target_g", data.get("beverage_target_g"), errors)
    if bev_target is not None and not bev_target > 0:
        errors.append("beverage_target_g must be positive when present")
    time_zero = data.get("time_zero", "pump_activation")
    if time_zero not in TIME_ZERO_KINDS:
        errors.append(f"time_zero must be one of {TIME_ZERO_KINDS}, got {time_zero!r}")
    # Grinder: ordinal setting, never particle diameter.
    grinder = data.get("grinder") or {}
    if isinstance(grinder, str):
        grinder = {"name": grinder}
    gname = grinder.get("name")
    gsetting = grinder.get("setting", data.get("grinder_setting"))
    if gsetting is not None and isinstance(gsetting, bool):
        errors.append("grinder.setting must be an ordinal number/string, not boolean")
    if isinstance(gsetting, str) and gsetting.strip() == "":
        gsetting = None
    # Temperature setting with location.
    temp = data.get("temperature") or {}
    tset = _require_finite("temperature.setting_C", temp.get("setting_C",
                           data.get("temp_setting_C")), errors)
    tloc = temp.get("location", data.get("temp_location"))
    if tset is not None and (not isinstance(tloc, str) or not tloc.strip()):
        errors.append("temperature.location is required when temperature.setting_C is given")
    # Pre-brew only if confirmed.
    prebrew = data.get("prebrew")
    if prebrew is not None and not isinstance(prebrew, dict):
        errors.append("prebrew must be an object or null (only if confirmed)")
    # Cup TDS with method/uncertainty; never infer EY without it.
    cup = data.get("cup") or {}
    tds = _require_finite("cup.TDS_pct", cup.get("TDS_pct", data.get("TDS_pct")), errors)
    tds_method = cup.get("method", data.get("TDS_method"))
    tds_unc = _require_finite("cup.TDS_uncertainty_pct",
                              cup.get("uncertainty_pct", data.get("TDS_uncertainty_pct")), errors)
    if tds is not None:
        if not 0.0 <= tds <= 30.0:
            errors.append(f"cup.TDS_pct={tds!r} outside plausible 0-30 % (inconsistent units?)")
    if tds_unc is not None and not tds_unc >= 0:
        errors.append("cup TDS uncertainty must be non-negative when present")
    # Samples: inline or external CSV rows.
    inline = data.get("samples")
    if series_rows is None:
        if inline is None:
            errors.append("no samples: provide inline 'samples' or external CSV rows")
            rows: list[dict] = []
        elif not isinstance(inline, list) or len(inline) == 0:
            errors.append("samples must be a non-empty list")
            rows = []
        else:
            rows = []
            for idx, s in enumerate(inline):
                if not isinstance(s, dict):
                    errors.append(f"samples[{idx}] must be an object")
                    continue
                rows.append({
                    "time_s": s.get("time_s"),
                    "beverage_mass_g": s.get("beverage_mass_g"),
                    "pressure_Pa": s.get("pressure_Pa"),
                    "pressure_location": s.get("pressure_location"),
                    "temperature_K": s.get("temperature_K"),
                    "temperature_location": s.get("temperature_location"),
                    "_lineno": idx + 1,
                })
    else:
        rows = [dict(r) for r in series_rows]
        if inline is not None:
            errors.append("provide either inline samples or external CSV, not both")
    # Per-row finite checks (before ordering checks).
    clean: list[dict] = []
    for r in rows:
        loc = r.get("_lineno", "?")
        t = r.get("time_s")
        m = r.get("beverage_mass_g")
        if isinstance(t, bool) or not isinstance(t, (int, float)) or not math.isfinite(float(t)):
            errors.append(f"row {loc}: time_s must be finite, got {t!r}")
            continue
        if isinstance(m, bool) or not isinstance(m, (int, float)) or not math.isfinite(float(m)):
            errors.append(f"row {loc}: beverage_mass_g must be finite, got {m!r}")
            continue
        t_f, m_f = float(t), float(m)
        if m_f < 0.0:
            errors.append(f"row {loc}: beverage_mass_g must be >= 0, got {m_f!r}")
        p = r.get("pressure_Pa")
        if p is not None:
            if isinstance(p, bool) or not isinstance(p, (int, float)) or not math.isfinite(float(p)):
                errors.append(f"row {loc}: pressure_Pa must be finite or missing, got {p!r}")
            else:
                p_f = float(p)
                if not 0.0 <= p_f <= 2.0e6:
                    errors.append(f"row {loc}: pressure_Pa={p_f!r} outside 0-2e6 Pa "
                                  "(inconsistent units? use Pa, not bar)")
                r["pressure_Pa"] = p_f
        tk = r.get("temperature_K")
        if tk is not None:
            if isinstance(tk, bool) or not isinstance(tk, (int, float)) or not math.isfinite(float(tk)):
                errors.append(f"row {loc}: temperature_K must be finite or missing, got {tk!r}")
            else:
                tk_f = float(tk)
                if not 250.0 <= tk_f <= 373.0:
                    errors.append(f"row {loc}: temperature_K={tk_f!r} outside 250-373 K "
                                  "(inconsistent units? use K, not C)")
                r["temperature_K"] = tk_f
        # Plausible mass/time ranges catch g/kg or ms/s slips.
        if not 0.0 <= t_f <= 600.0:
            errors.append(f"row {loc}: time_s={t_f!r} outside 0-600 s (inconsistent units?)")
        if not 0.0 <= m_f <= 1000.0:
            errors.append(f"row {loc}: beverage_mass_g={m_f!r} outside 0-1000 g")
        r["time_s"], r["beverage_mass_g"] = t_f, m_f
        clean.append(r)
    # Ordering / duplicates.
    if aggregate is not None and aggregate not in AGGREGATIONS:
        errors.append(f"aggregate must be one of {AGGREGATIONS}, got {aggregate!r}")
    if not errors and clean:
        times = [r["time_s"] for r in clean]
        seen: dict[float, int] = {}
        dups: set[float] = set()
        for t in times:
            if t in seen:
                dups.add(t)
            seen[t] = seen.get(t, 0) + 1
        if dups and aggregate is None:
            errors.append(f"duplicate timestamps {sorted(dups)} without explicit "
                          "aggregation rule (pass --aggregate mean|first|last|median)")
        elif dups:
            assert aggregate is not None
            clean = _aggregate_duplicates(clean, aggregate)
            times = [r["time_s"] for r in clean]
        for a, b in zip(times, times[1:]):
            if b < a:
                errors.append(f"decreasing times: {a!r} followed by {b!r}")
                break
            if b == a and aggregate is None:
                # Already reported as duplicate; avoid double count.
                pass
    # Dose plausibility (g, not kg/mg).
    if dose_g is not None and not 2.0 <= dose_g <= 40.0:
        errors.append(f"dose_g={dose_g!r} outside plausible 2-40 g (inconsistent units? use g)")
    if errors:
        raise ShotError("Invalid shot " + (str(shot_id) if isinstance(shot_id, str) else "?")
                        + ": " + "; ".join(errors))
    assert dose_g is not None
    meta = ShotMetadata(
        shot_id=str(shot_id).strip(),
        source=str(source),
        machine_profile_id=str(profile_id).strip(),
        feed=str(feed).strip(),
        dose_g=float(dose_g),
        beverage_target_g=float(bev_target) if bev_target is not None else None,
        basket_id=data.get("basket_id") or (data.get("basket") or {}).get("id"),
        bean_batch=(data.get("beans") or {}).get("batch", data.get("bean_batch")),
        roast_date=(data.get("beans") or {}).get("roast_date", data.get("roast_date")),
        grinder=gname if isinstance(gname, str) else data.get("grinder_name"),
        grinder_setting=gsetting,
        grinder_setting_is_ordinal=True,
        temp_setting_C=float(tset) if tset is not None else None,
        temp_location=str(tloc).strip() if isinstance(tloc, str) and tloc.strip() else None,
        prebrew=prebrew,
        time_zero=str(time_zero),
        notes=str(data.get("notes", "")),
        tds_pct=float(tds) if tds is not None else None,
        tds_method=str(tds_method) if isinstance(tds_method, str) and tds_method.strip() else None,
        tds_uncertainty_pct=float(tds_unc) if tds_unc is not None else None,
    )
    samples = [ShotSample(
        time_s=float(r["time_s"]),
        beverage_mass_g=float(r["beverage_mass_g"]),
        pressure_Pa=float(r["pressure_Pa"]) if r.get("pressure_Pa") is not None else None,
        pressure_location=r.get("pressure_location"),
        temperature_K=float(r["temperature_K"]) if r.get("temperature_K") is not None else None,
        temperature_location=r.get("temperature_location"),
    ) for r in clean]
    return Shot(metadata=meta, samples=samples)


def load_shot_file(json_path: str, series_csv: str | None = None,
                   aggregate: str | None = None) -> Shot:
    """Load and validate one shot JSON file, optionally joined with a CSV series."""
    with open(json_path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    rows = load_series_csv(series_csv) if series_csv is not None else None
    return parse_shot(data, rows, aggregate)


def brew_ratio(beverage_mass_g: float, dose_g: float) -> float:
    """Beverage-to-dose ratio [-]; both in grams."""
    if not (math.isfinite(beverage_mass_g) and math.isfinite(dose_g) and dose_g > 0):
        raise ShotError(f"brew_ratio needs finite masses and positive dose, "
                        f"got beverage={beverage_mass_g!r} dose={dose_g!r}")
    return float(beverage_mass_g) / float(dose_g)


def cup_ey_pct(beverage_mass_g: float, tds_pct: float | None, dose_g: float) -> float | None:
    """Measured cup EY %; None when TDS is missing (never inferred)."""
    if tds_pct is None:
        return None
    if not (math.isfinite(beverage_mass_g) and math.isfinite(tds_pct) and math.isfinite(dose_g)):
        raise ShotError("cup EY needs finite inputs")
    if not dose_g > 0:
        raise ShotError("cup EY needs positive dose")
    if not 0.0 <= tds_pct <= 30.0:
        raise ShotError(f"TDS {tds_pct!r} % outside 0-30")
    return 100.0 * float(beverage_mass_g) * (float(tds_pct) / 100.0) / float(dose_g)


def negative_increment_flags(samples: list[ShotSample]) -> list[int]:
    """Indices i>0 where raw mass decreases (balance noise); raw is retained."""
    flags: list[int] = []
    for i in range(1, len(samples)):
        if samples[i].beverage_mass_g < samples[i - 1].beverage_mass_g:
            flags.append(i)
    return flags


def derive_flow(samples: list[ShotSample], window_s: float = 1.0) -> tuple[list[float], dict]:
    """Derive beverage flow (g/s) by local linear regression in a time window.

    For each sample i, all samples with ``|t - t_i| <= window_s/2`` are fit by
    least squares ``m = a + q*t``; ``q`` is the derived (smoothed) flow at i.
    Endpoints use one-sided windows (only available side). Returns
    ``(flows, info)`` where info reports ``window_s``, method
    ``local_linear_regression``, endpoint treatment ``one_sided``, the number
    of flagged raw negative increments, and a reminder that this is derived,
    not directly measured, flow. Raw samples are never modified.
    """
    if not math.isfinite(window_s) or not window_s > 0:
        raise ShotError(f"flow window_s must be finite-positive, got {window_s!r}")
    if len(samples) < 2:
        raise ShotError("need >= 2 samples to derive flow")
    times = [s.time_s for s in samples]
    masses = [s.beverage_mass_g for s in samples]
    half = window_s / 2.0
    flows: list[float] = []
    for i, ti in enumerate(times):
        idx = [j for j, tj in enumerate(times) if abs(tj - ti) <= half + 1e-12]
        if len(idx) < 2:
            idx = [max(0, i - 1), i] if i > 0 else [i, min(len(times) - 1, i + 1)]
            if idx[0] == idx[1]:
                flows.append(0.0)
                continue
        tj = [times[j] for j in idx]
        mj = [masses[j] for j in idx]
        tbar = sum(tj) / len(tj)
        mbar = sum(mj) / len(mj)
        denom = sum((t - tbar) ** 2 for t in tj)
        if denom == 0.0:
            flows.append(0.0)
        else:
            flows.append(sum((t - tbar) * (m - mbar) for t, m in zip(tj, mj)) / denom)
    info = {
        "method": "local_linear_regression",
        "window_s": float(window_s),
        "endpoint_treatment": "one_sided_window",
        "n_negative_raw_increments": len(negative_increment_flags(samples)),
        "note": ("Derived (smoothed) flow, not directly measured flow. "
                 "Raw mass samples retained; small negative increments are "
                 "balance noise and are flagged, not removed."),
    }
    return flows, info


def constant_flow_baseline(times: list[float], q_g_s: float,
                           t0_s: float, m0_g: float) -> list[float]:
    """Constant-flow mass trace m(t) = m0 + q*(t - t0) for t >= t0 (else m0)."""
    if not (math.isfinite(q_g_s) and math.isfinite(t0_s) and math.isfinite(m0_g)):
        raise ShotError("baseline needs finite q, t0, m0")
    return [float(m0_g) + float(q_g_s) * max(0.0, t - float(t0_s)) for t in times]


def shot_summary(shot: Shot, flow_window_s: float = 1.0) -> dict:
    """Brew ratio, measured cup EY (when TDS present), flows, and flags."""
    final_mass = shot.samples[-1].beverage_mass_g if shot.samples else 0.0
    ratio = brew_ratio(final_mass, shot.metadata.dose_g)
    ey = cup_ey_pct(final_mass, shot.metadata.tds_pct, shot.metadata.dose_g)
    flows, finfo = derive_flow(shot.samples, flow_window_s)
    return {
        "schema_version": SHOT_SCHEMA_VERSION,
        "shot_id": shot.metadata.shot_id,
        "source": shot.metadata.source,
        "time_zero": shot.metadata.time_zero,
        "dose_g": shot.metadata.dose_g,
        "beverage_target_g": shot.metadata.beverage_target_g,
        "final_beverage_mass_g": final_mass,
        "brew_ratio": ratio,
        "cup_TDS_pct": shot.metadata.tds_pct,
        "cup_TDS_method": shot.metadata.tds_method,
        "cup_TDS_uncertainty_pct": shot.metadata.tds_uncertainty_pct,
        "cup_EY_pct": ey,
        "cup_EY_status": ("measured" if ey is not None
                          else "missing_no_TDS_no_inference"),
        "n_samples": len(shot.samples),
        "negative_raw_increments": negative_increment_flags(shot.samples),
        "derived_flow_g_s": flows,
        "derived_flow_info": finfo,
        "calibration_status": CALIBRATION_STATUS,
        "limitations": LIMITATIONS_NOTE,
    }


def repeatability(shots: list[Shot]) -> dict:
    """Mean/std of final mass and brew ratio across repeated shots."""
    import statistics
    finals = [s.samples[-1].beverage_mass_g for s in shots if s.samples]
    ratios = [brew_ratio(f, s.metadata.dose_g) for f, s in zip(finals, shots)]
    eys = [cup_ey_pct(f, s.metadata.tds_pct, s.metadata.dose_g)
           for f, s in zip(finals, shots)]
    eys_known = [e for e in eys if e is not None]
    return {
        "n_shots": len(shots),
        "shot_ids": [s.metadata.shot_id for s in shots],
        "final_mass_g": {
            "mean": statistics.fmean(finals) if finals else None,
            "stdev": statistics.stdev(finals) if len(finals) >= 2 else 0.0,
            "values": finals,
        },
        "brew_ratio": {
            "mean": statistics.fmean(ratios) if ratios else None,
            "stdev": statistics.stdev(ratios) if len(ratios) >= 2 else 0.0,
            "values": ratios,
        },
        "cup_EY_pct": {
            "mean": statistics.fmean(eys_known) if eys_known else None,
            "stdev": statistics.stdev(eys_known) if len(eys_known) >= 2 else 0.0,
            "values": eys,
            "n_missing_TDS": sum(1 for e in eys if e is None),
        },
    }


def fit_constant_flow_train(shots: list[Shot], t0_s: float) -> dict:
    """Fit constant flow q on training shots only (whole shots, not rows).

    q is the mean slope (final-initial)/(t_final-t0) over the already-wet
    interval [t0, end] across training shots. Hold out whole shots for testing.
    """
    slopes: list[float] = []
    for s in shots:
        ts = [p.time_s for p in s.samples if p.time_s >= t0_s]
        ms = [p.beverage_mass_g for p, q in zip(s.samples, [p.time_s for p in s.samples]) if q >= t0_s]
        if len(ts) >= 2 and (ts[-1] - ts[0]) > 0:
            slopes.append((ms[-1] - ms[0]) / (ts[-1] - ts[0]))
    if not slopes:
        raise ShotError(f"no training data at/after t0={t0_s!r}")
    import statistics
    q = statistics.fmean(slopes)
    return {"q_g_s": q, "t0_s": t0_s, "n_train": len(shots),
            "train_slopes": slopes,
            "stdev": statistics.stdev(slopes) if len(slopes) >= 2 else 0.0}


def residuals_vs_baseline(shot: Shot, q_g_s: float, t0_s: float, m0_g: float) -> dict:
    """Residuals of a shot's raw mass trace versus a constant-flow baseline."""
    times = [p.time_s for p in shot.samples]
    masses = [p.beverage_mass_g for p in shot.samples]
    pred = constant_flow_baseline(times, q_g_s, t0_s, m0_g)
    resid = [m - p for m, p in zip(masses, pred)]
    import statistics
    return {
        "shot_id": shot.metadata.shot_id,
        "q_g_s": q_g_s, "t0_s": t0_s, "m0_g": m0_g,
        "max_abs_resid_g": max(abs(r) for r in resid) if resid else 0.0,
        "rmse_g": (sum(r * r for r in resid) / len(resid)) ** 0.5 if resid else 0.0,
        "mean_resid_g": statistics.fmean(resid) if resid else 0.0,
        "residuals_g": resid,
    }


def compare_shots(train: list[Shot], heldout: list[Shot], t0_s: float) -> dict:
    """Train constant-flow baseline on train shots; evaluate on held-out whole shots."""
    if not train:
        raise ShotError("compare needs >= 1 training shot")
    fit = fit_constant_flow_train(train, t0_s)
    q, t0 = fit["q_g_s"], fit["t0_s"]
    results: list[dict] = []
    for s in heldout:
        # m0: mass at/after t0 (first sample >= t0), never fitted per-row.
        after = [(p.time_s, p.beverage_mass_g) for p in s.samples if p.time_s >= t0]
        if not after:
            raise ShotError(f"held-out shot {s.metadata.shot_id} has no data >= t0={t0}")
        m0 = after[0][1]
        results.append(residuals_vs_baseline(s, q, t0, m0))
    return {
        "schema_version": SHOT_SCHEMA_VERSION,
        "fit": fit,
        "train_repeatability": repeatability(train),
        "heldout_repeatability": repeatability(heldout) if heldout else None,
        "heldout_residuals": results,
        "calibration_status": CALIBRATION_STATUS,
        "limitations": LIMITATIONS_NOTE + " Baseline fit uses whole training shots only.",
    }


def write_svg_plot(path: str, series: list[tuple[str, list[float], list[float], str]],
                   xlabel: str = "time_s (s)",
                   ylabel: str = "value (see title)",
                   title: str = "plot") -> None:
    """Write a minimal SVG line plot (stdlib only, no matplotlib).

    ``series`` is a list of (label, xs, ys, kind) where kind is
    'observed' (solid) or 'simulated' (dashed). Axes carry explicit
    units via ``xlabel``/``ylabel``; observations vs simulations are
    distinguished by style + legend.
    """
    import math as _math
    allx = [x for _, xs, _, _ in series for x in xs]
    ally = [y for _, _, ys, _ in series for y in ys if _math.isfinite(y)]
    if not allx or not ally:
        raise ShotError("no plottable data")
    x0, x1 = min(allx), max(allx)
    y0, y1 = min(ally), max(ally)
    if x1 == x0:
        x1 = x0 + 1.0
    if y1 == y0:
        y1 = y0 + 1.0
    padx = 0.05 * (x1 - x0)
    pady = 0.10 * (y1 - y0)
    x0, x1, y0, y1 = x0 - padx, x1 + padx, y0 - pady, y1 + pady
    W, H, ML, MR, MT, MB = 640, 420, 64, 16, 28, 52
    pw, ph = W - ML - MR, H - MT - MB

    def sx(x: float) -> float:
        return ML + (x - x0) / (x1 - x0) * pw

    def sy(y: float) -> float:
        return MT + ph - (y - y0) / (y1 - y0) * ph

    import html as _html
    colors = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e"]
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" role="img">']
    parts.append(f'<title>{_html.escape(title)}</title>')
    parts.append(f'<rect x="0" y="0" width="{W}" height="{H}" fill="white"/>')
    parts.append(f'<text x="{W/2:.0f}" y="16" font-size="12" text-anchor="middle">{_html.escape(title)}</text>')
    parts.append(f'<rect x="{ML}" y="{MT}" width="{pw}" height="{ph}" fill="none" stroke="black"/>')
    parts.append(f'<text x="{ML + pw/2:.0f}" y="{H - 8}" font-size="11" text-anchor="middle">{_html.escape(xlabel)}</text>')
    parts.append(f'<text x="14" y="{MT + ph/2:.0f}" font-size="11" text-anchor="middle" '
                 f'transform="rotate(-90 14,{MT + ph/2:.0f})">{_html.escape(ylabel)}</text>')
    # gridlines + ticks (3 x, 3 y)
    for k in range(4):
        xv = x0 + (x1 - x0) * k / 3.0
        parts.append(f'<line x1="{sx(xv):.1f}" y1="{MT}" x2="{sx(xv):.1f}" y2="{MT + ph}" stroke="#ddd"/>')
        parts.append(f'<text x="{sx(xv):.1f}" y="{MT + ph + 16}" font-size="10" text-anchor="middle">{xv:.2g}</text>')
        yv = y0 + (y1 - y0) * k / 3.0
        parts.append(f'<line x1="{ML}" y1="{sy(yv):.1f}" x2="{ML + pw}" y2="{sy(yv):.1f}" stroke="#ddd"/>')
        parts.append(f'<text x="{ML - 6}" y="{sy(yv) + 3:.1f}" font-size="10" text-anchor="end">{yv:.2g}</text>')
    for idx, (label, xs, ys, kind) in enumerate(series):
        col = colors[idx % len(colors)]
        dash = ' stroke-dasharray="6,3"' if kind == "simulated" else ""
        pts = " ".join(f"{sx(x):.1f},{sy(y):.1f}" for x, y in zip(xs, ys) if _math.isfinite(y))
        parts.append(f'<polyline points="{pts}" fill="none" stroke="{col}" stroke-width="2"{dash}/>')
        parts.append(f'<text x="{ML + pw - 4}" y="{MT + 16 + idx * 16}" font-size="11" '
                     f'text-anchor="end" fill="{col}">{label} [{kind}]</text>')
    parts.append("</svg>")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(parts) + "\n")
