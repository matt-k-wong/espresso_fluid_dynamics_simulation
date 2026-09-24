"""M6 local experiment explorer: save / reproduce a shot comparison bundle.

A bundle is a self-contained directory that reproduces one
``shot-compare`` result from its own ``inputs/`` copies, distinguishes
measured / synthetic / illustrative / fitted-empirical inputs, exports tidy
CSVs, and states limitations. No new fluid physics, no machine control, no
permeability fitting. Synthetic fixtures only validate the harness.
"""

from __future__ import annotations

import csv
import datetime
import hashlib
import json
import math
from pathlib import Path

from .shots import (CALIBRATION_STATUS, LIMITATIONS_NOTE, SHOT_SCHEMA_VERSION,
                    Shot, ShotError, compare_shots, constant_flow_baseline,
                    derive_flow, load_shot_file, repeatability, shot_summary,
                    write_svg_plot)

EXPERIMENT_SCHEMA_VERSION = "experiment_v1"


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _write_json_strict(path: Path, obj: object) -> None:
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n")


def _classify_shot(shot: Shot) -> str:
    return "measured" if shot.metadata.source == "measured" else "synthetic"


def save_experiment(shot_paths: list[str], series_paths: list[str] | None,
                    already_wet_t0: float, holdout: list[str],
                    flow_window: float, out_dir: str,
                    label: str = "") -> dict:
    """Validate shots, run whole-shot holdout comparison, write bundle."""
    if series_paths and len(series_paths) != len(shot_paths):
        raise ShotError("--series length must match --shots length")
    if not math.isfinite(already_wet_t0):
        raise ShotError(f"already-wet-t0 must be finite, got {already_wet_t0!r}")
    if not math.isfinite(flow_window) or not flow_window > 0:
        raise ShotError(f"flow-window must be finite-positive, got {flow_window!r}")
    shots: list[Shot] = []
    for i, sp in enumerate(shot_paths):
        cp = series_paths[i] if series_paths else None
        shots.append(load_shot_file(sp, cp, None))
    ids = [s.metadata.shot_id for s in shots]
    if len(set(ids)) != len(ids):
        raise ShotError(f"duplicate shot_id in inputs: {ids}")
    hold = list(holdout or [])
    if hold:
        unknown = sorted(set(hold) - set(ids))
        if unknown:
            raise ShotError(f"holdout ids match no loaded shots: {unknown}")
        train = [s for s in shots if s.metadata.shot_id not in set(hold)]
        held = [s for s in shots if s.metadata.shot_id in set(hold)]
    else:
        if len(shots) < 2:
            raise ShotError("need >= 2 shots to compare (or pass --holdout)")
        train, held = shots[:-1], shots[-1:]
        hold = [held[0].metadata.shot_id]
    result = compare_shots(train, held, already_wet_t0)
    q, t0 = result["fit"]["q_g_s"], result["fit"]["t0_s"]

    out = Path(out_dir)
    inputs_shots = out / "inputs" / "shots"
    inputs_series = out / "inputs" / "series"
    export = out / "export"
    for d in (inputs_shots, inputs_series, export):
        d.mkdir(parents=True, exist_ok=True)
    # Copy inputs verbatim; record hashes for reproduce check.
    input_files: list[dict] = []
    for i, sp in enumerate(shot_paths):
        src = Path(sp)
        dst = inputs_shots / src.name
        dst.write_bytes(src.read_bytes())
        entry: dict = {"role": "shot_json", "shot_id": ids[i],
                       "src": str(sp), "bundle_path": f"inputs/shots/{src.name}",
                       "sha256": _sha256_file(dst),
                       "provenance": _classify_shot(shots[i])}
        if series_paths:
            csrc = Path(series_paths[i])
            cdst = inputs_series / csrc.name
            cdst.write_bytes(csrc.read_bytes())
            entry["series_src"] = str(series_paths[i])
            entry["series_bundle_path"] = f"inputs/series/{csrc.name}"
            entry["series_sha256"] = _sha256_file(cdst)
        input_files.append(entry)

    _write_json_strict(out / "comparison.json", result)
    mass_series, flow_series = [], []
    for s in shots:
        times = [p.time_s for p in s.samples]
        masses = [p.beverage_mass_g for p in s.samples]
        mass_series.append((f"{s.metadata.shot_id} mass (observed)",
                            times, masses, "observed"))
        after = [(p.time_s, p.beverage_mass_g) for p in s.samples if p.time_s >= t0]
        m0 = after[0][1] if after else masses[0]
        mass_series.append((f"{s.metadata.shot_id} baseline (simulated)",
                            times, constant_flow_baseline(times, q, t0, m0),
                            "simulated"))
        flows, _ = derive_flow(s.samples, flow_window)
        flow_series.append((f"{s.metadata.shot_id} derived flow (observed)",
                            times, flows, "observed"))
    write_svg_plot(str(out / "mass.svg"), mass_series,
                   xlabel="time_s (s)", ylabel="beverage_mass_g (g)",
                   title=f"{label or 'experiment'} mass vs time")
    write_svg_plot(str(out / "flow.svg"), flow_series,
                   xlabel="time_s (s)", ylabel="derived_flow_g_s (g/s)",
                   title=f"{label or 'experiment'} derived flow")
    # Tidy exports (raw retained; derived flow labeled derived).
    with open(export / "combined_series.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["shot_id", "provenance", "time_s",
                                           "beverage_mass_g"])
        w.writeheader()
        for s in shots:
            prov = _classify_shot(s)
            for p in s.samples:
                w.writerow({"shot_id": s.metadata.shot_id, "provenance": prov,
                            "time_s": p.time_s,
                            "beverage_mass_g": p.beverage_mass_g})
    with open(export / "metrics.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["shot_id", "provenance", "dose_g",
                                           "final_mass_g", "brew_ratio",
                                           "cup_TDS_pct", "cup_EY_pct"])
        w.writeheader()
        for s in shots:
            summ = shot_summary(s, flow_window)
            w.writerow({"shot_id": s.metadata.shot_id,
                        "provenance": _classify_shot(s),
                        "dose_g": s.metadata.dose_g,
                        "final_mass_g": summ["final_beverage_mass_g"],
                        "brew_ratio": summ["brew_ratio"],
                        "cup_TDS_pct": ("" if summ["cup_TDS_pct"] is None
                                        else summ["cup_TDS_pct"]),
                        "cup_EY_pct": ("" if summ["cup_EY_pct"] is None
                                       else summ["cup_EY_pct"])})
    manifest = {
        "experiment_schema_version": EXPERIMENT_SCHEMA_VERSION,
        "shot_schema_version": SHOT_SCHEMA_VERSION,
        "label": label,
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "already_wet_t0_s": already_wet_t0,
        "flow_window_s": flow_window,
        "train_ids": [s.metadata.shot_id for s in train],
        "holdout_ids": [s.metadata.shot_id for s in held],
        "inputs": input_files,
        "inputs_classification": {
            "measured": "shot traces with source=measured (observations).",
            "synthetic": "shot traces with source=synthetic (harness tests only).",
            "illustrative": "simulator constants/assumptions (not fitted).",
            "fitted_empirical": ("constant-flow q baseline fitted on TRAIN whole "
                                 "shots only; not a permeability or physics fit."),
        },
        "outputs": ["comparison.json", "mass.svg", "flow.svg",
                    "export/combined_series.csv", "export/metrics.csv"],
        "limitations": LIMITATIONS_NOTE + " Baseline fit uses whole training shots only.",
        "calibration_status": CALIBRATION_STATUS,
        "reproduce_cmd": "experiment-reproduce --dir <bundle>",
    }
    _write_json_strict(out / "manifest.json", manifest)
    readme = (
        f"Experiment: {label or '(unlabeled)'}\n"
        f"Bundle reproduces one whole-shot-holdout constant-flow comparison.\n"
        f"Train: {', '.join(manifest['train_ids'])} | "
        f"Holdout: {', '.join(manifest['holdout_ids'])} | "
        f"t0={already_wet_t0:g} s q={q:.4f} g/s\n"
        "Provenance: per-shot measured/synthetic in manifest inputs; "
        "q is fitted-empirical (train shots only), simulator constants are "
        "illustrative. Observed=solid, simulated=dashed.\n"
        f"Limitations: {manifest['limitations']}\n"
        f"Calibration: {CALIBRATION_STATUS}\n"
        "Reproduce: PYTHONPATH=src ./.venv/bin/python -m espresso_m1.cli "
        "experiment-reproduce --dir <bundle>\n"
        "Export: export/combined_series.csv (raw) + export/metrics.csv.\n"
    )
    (out / "README.txt").write_text(readme)
    return {"manifest": manifest, "comparison": result, "out_dir": str(out)}


def _close(a: float | None, b: float | None, rel=1e-9, abs_=1e-12) -> bool:
    if a is None or b is None:
        return a is b
    if not (math.isfinite(a) and math.isfinite(b)):
        return False
    return abs(a - b) <= max(abs_ , rel * max(abs(a), abs(b)))


def reproduce_experiment(bundle_dir: str) -> dict:
    """Re-hash inputs, recompute comparison from bundle inputs/, check match."""
    root = Path(bundle_dir)
    try:
        manifest = json.loads((root / "manifest.json").read_text())
        saved = json.loads((root / "comparison.json").read_text())
    except (FileNotFoundError, json.JSONDecodeError) as exc:
        raise ShotError(f"bundle unreadable: {exc}") from exc
    if manifest.get("experiment_schema_version") != EXPERIMENT_SCHEMA_VERSION:
        raise ShotError("unsupported experiment_schema_version")
    # Verify stored hashes before trusting anything.
    hash_ok = True
    hash_notes: list[str] = []
    for entry in manifest.get("inputs", []):
        for key, relp in (("sha256", entry.get("bundle_path")),
                          ("series_sha256", entry.get("series_bundle_path"))):
            if relp is None:
                continue
            p = root / relp
            if not p.exists():
                hash_ok = False
                hash_notes.append(f"missing {relp}")
                continue
            if _sha256_file(p) != entry.get(key):
                hash_ok = False
                hash_notes.append(f"hash mismatch {relp}")
    # Reload from bundle copies only (self-contained).
    shots: list[Shot] = []
    for entry in manifest.get("inputs", []):
        jp = root / entry["bundle_path"]
        cp = root / entry["series_bundle_path"] if entry.get("series_bundle_path") else None
        shots.append(load_shot_file(str(jp), str(cp) if cp else None, None))
    by_id = {s.metadata.shot_id: s for s in shots}
    train = [by_id[i] for i in manifest["train_ids"]]
    held = [by_id[i] for i in manifest["holdout_ids"]]
    recomputed = compare_shots(train, held, manifest["already_wet_t0_s"])
    checks: dict = {}
    checks["q_match"] = _close(recomputed["fit"]["q_g_s"],
                               saved["fit"]["q_g_s"])
    r1 = (recomputed["heldout_residuals"][0]["rmse_g"]
          if recomputed["heldout_residuals"] else 0.0)
    r0 = (saved["heldout_residuals"][0]["rmse_g"]
          if saved.get("heldout_residuals") else 0.0)
    checks["heldout_rmse_match"] = _close(r1, r0)
    t1 = recomputed["train_repeatability"]["final_mass_g"]["mean"]
    t0m = saved["train_repeatability"]["final_mass_g"]["mean"]
    checks["train_mean_match"] = _close(t1, t0m)
    checks["hash_ok"] = hash_ok
    passed = all(checks.values())
    report = {"passed": passed, "checks": checks, "hash_notes": hash_notes,
              "recomputed_fit_q": recomputed["fit"]["q_g_s"],
              "saved_fit_q": saved["fit"]["q_g_s"],
              "limitations": saved.get("limitations", LIMITATIONS_NOTE),
              "calibration_status": CALIBRATION_STATUS}
    return report
