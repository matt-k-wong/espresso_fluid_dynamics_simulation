"""Minimal CLI: run one case or a pressure sweep (all file output lives here)."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
from dataclasses import asdict
from pathlib import Path

import numpy as np

from .driver import SweepRow, run_case, run_sweep
from .params import BAR_PA, SimParams, apply_overrides, from_json_file
from .tparams import TransportParams, apply_overrides as apply_toverrides
from .tparams import from_dict as tparams_from_dict
from .transport import ASSUMPTIONS, run_transport

SWEEP_FIELDS = [
    "p_in_Pa",
    "p_in_bar",
    "Q_m3s",
    "Q_mLs",
    "converged",
    "iters",
    "mismatch_abs_m3s",
    "mismatch_rel",
    "max_residual_abs_m3s",
    "max_residual_rel",
    "delta_p_Pa",
    "R_hyd",
    "R_hyd_status",
    "k_min",
    "k_max",
    "phi_min",
    "phi_max",
    "open_area_m2",
    "open_area_ideal_m2",
    "error",
]

# Artifacts owned by `run` inside --out-dir. Only these names are ever
# created or deleted by this command; nothing else in the directory is touched.
OWNED_RUN_FILES = ("summary.json", "fields.npz", "failure.json")

# Artifacts owned by `transport` inside --out-dir.
OWNED_TRANSPORT_FILES = (
    "transport_summary.json",
    "transport_timeseries.csv",
    "transport_fields.npz",
)

TRANSPORT_SERIES_FIELDS = [
    "t_s", "dt_s", "C_inlet_kgm3", "Q_out_m3s", "M_solute_cup_kg",
    "V_cup_m3", "TDS_pct", "EY_pct", "outlet_C_kgm3", "sum_B_kg",
    "sum_WC_kg", "cumulative_inlet_kg", "ledger_signed_kg",
    "lin_rel", "min_C_kgm3", "min_B_kg",
]


def parse_bar_list(spec: str) -> list[float]:
    """Parse ``'a,b,c'`` or ``'start:stop:step'`` in bar; return Pa values.

    Ranges include the start, advance by ``step``, and never pass ``stop``
    (``1:2:0.6`` means [1, 1.6] bar). Descending ranges need a negative step
    (``3:1:-1`` means [3, 2, 1]). Zero steps, direction-inconsistent steps,
    non-finite values, and empty results are rejected with ``ValueError``.
    """
    text = spec.strip()
    if not text:
        raise ValueError("Empty pressure specification")
    if ":" in text:
        chunks = text.split(":")
        if len(chunks) != 3:
            raise ValueError(f"Range {text!r} must have the form start:stop:step")
        try:
            start, stop, step = (float(x) for x in chunks)
        except ValueError:
            raise ValueError(f"Range {text!r} contains a non-numeric value")
        for value in (start, stop, step):
            if not math.isfinite(value):
                raise ValueError(f"Range {text!r} contains a non-finite value")
        if step == 0.0:
            raise ValueError(f"Range {text!r} has a zero step")
        if (stop - start) * step < 0.0:
            raise ValueError(
                f"Range {text!r}: step sign is inconsistent with start->stop "
                "(use a negative step for descending ranges)"
            )
        eps = 1e-9 * max(1.0, abs(start), abs(stop), abs(step))
        bars: list[float] = []
        n = 0
        while n <= 1_000_000:
            current = start + n * step
            if step > 0.0 and current > stop + eps:
                break
            if step < 0.0 and current < stop - eps:
                break
            bars.append(current)
            n += 1
        if n > 1_000_000:
            raise ValueError(f"Range {text!r} is too long")
    else:
        try:
            bars = [float(x) for x in text.split(",")]
        except ValueError:
            raise ValueError(f"List {text!r} contains a non-numeric value")
        for value in bars:
            if not math.isfinite(value):
                raise ValueError(f"List {text!r} contains a non-finite value")
    if not bars:
        raise ValueError(f"Specification {text!r} yields no pressures")
    return [float(v) * BAR_PA for v in bars]


def _stage_publish(out_dir: Path, payloads: list[tuple[str, object]]) -> None:
    """Write staged files and publish them in order (status/summary last).

    Each ``(name, writer)`` payload is written to ``name.tmp`` then atomically
    moved to ``name`` via ``os.replace``. Callers order payloads so the
    summary/status file is published last.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, writer in payloads:
        tmp = out_dir / (name + ".tmp")
        try:
            writer(tmp)
            os.replace(tmp, out_dir / name)
        finally:
            try:
                os.unlink(tmp)
            except OSError:
                pass


def _write_json_strict(path: Path, obj: object) -> None:
    """Write JSON that strict parsers accept (NaN/Infinity raise, never emit)."""
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n")


def _remove_owned(out_dir: Path, names: tuple[str, ...]) -> None:
    for name in names:
        try:
            (out_dir / name).unlink()
        except FileNotFoundError:
            pass


def _publish_run_invalid(out_dir: Path, error: Exception,
                         args: argparse.Namespace,
                         params_dict: object = None) -> None:
    """Record an invalid-input run as a failed new run.

    Removes owned stale success artifacts (summary/fields) and writes a
    strict-JSON failure record carrying the current input context (config
    path and overrides), never old parameters. Unrelated files are untouched.
    Only expected input/config failures reach here; programming errors
    propagate instead.
    """
    payload = {
        "converged": False,
        "error": str(error),
        "error_type": type(error).__name__,
        "iters": 0,
        "params": params_dict,
        "config": getattr(args, "config", None),
        "set": getattr(args, "set", None),
        "p_in_bar": getattr(args, "p_in_bar", None),
    }
    _remove_owned(out_dir, ("summary.json", "fields.npz"))
    _stage_publish(out_dir, [(
        "failure.json",
        lambda tmp: _write_json_strict(tmp, payload),
    )])


def _publish_transport_invalid(out_dir: Path, error: Exception,
                               args: argparse.Namespace,
                               hparams_dict: object = None,
                               tparams_dict: object = None) -> None:
    """Record an invalid transport invocation as a failed new run.

    Removes owned stale success artifacts and writes a strict-JSON failure
    status with the current hydraulic/transport input context. Unrelated
    files are preserved. Only expected input/config failures reach here.
    """
    from .transport import ASSUMPTIONS as _ASSUMPTIONS

    payload = {
        "converged": False,
        "error": str(error),
        "error_type": type(error).__name__,
        "n_steps": 0,
        "hydraulic_params": hparams_dict,
        "transport_params": tparams_dict,
        "hydraulic_config": getattr(args, "config", None),
        "transport_config": getattr(args, "tconfig", None),
        "hydraulic_overrides": getattr(args, "set", None),
        "transport_overrides": getattr(args, "tset", None),
        "assumptions": list(_ASSUMPTIONS),
    }
    _remove_owned(out_dir, ("transport_summary.json",
                            "transport_timeseries.csv",
                            "transport_fields.npz"))
    _stage_publish(out_dir, [(
        "transport_summary.json",
        lambda tmp: _write_json_strict(tmp, payload),
    )])


def _case_params(args: argparse.Namespace) -> SimParams:
    params = from_json_file(args.config)
    if args.set:
        params = apply_overrides(params, args.set)
    if getattr(args, "p_in_bar", None) is not None:
        params = apply_overrides(params, [f"p_in={args.p_in_bar * BAR_PA}"])
    return params


def cmd_run(args: argparse.Namespace) -> int:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    try:
        params = _case_params(args)
    except (ValueError, KeyError, TypeError, OSError) as exc:
        # Invalid input is a failed new run in the requested directory:
        # drop stale success, keep strict JSON + current input context.
        _publish_run_invalid(out_dir, exc, args, None)
        print(f"INVALID INPUT: {exc}", file=sys.stderr)
        return 2
    try:
        case = run_case(params)
    except (ValueError, KeyError) as exc:
        # Physical validation failure after parsing (e.g. under-resolved
        # opening): same failed-run contract with current params context.
        try:
            params_dict = asdict(params)
        except Exception:
            params_dict = None
        _publish_run_invalid(out_dir, exc, args, params_dict)
        print(f"INVALID INPUT: {exc}", file=sys.stderr)
        return 2
    if not case.converged:
        # No stale success output may masquerade as this run's result.
        _remove_owned(out_dir, ("summary.json", "fields.npz"))
        _stage_publish(out_dir, [(
            "failure.json",
            lambda tmp: _write_json_strict(tmp, {
                "converged": False,
                "iters": case.iters,
                "error": case.error,
                "error_type": "ConvergenceError",
                "params": asdict(params),
                "config": getattr(args, "config", None),
                "set": getattr(args, "set", None),
                "p_in_bar": getattr(args, "p_in_bar", None),
            }),
        )])
        print(f"FAILED TO CONVERGE: {case.error}", file=sys.stderr)
        return 1
    assert case.sol is not None and case.diag is not None
    d = case.diag
    summary = {
        "converged": True,
        "iters": case.iters,
        "Q_m3s": d.Q_out_m3s,
        "Q_mLs": d.Q_out_mLs,
        "mismatch_abs_m3s": d.mismatch_abs_m3s,
        "mismatch_rel": d.mismatch_rel,
        "max_residual_abs_m3s": d.max_residual_abs_m3s,
        "max_residual_rel": d.max_residual_rel,
        "delta_p_Pa": d.delta_p_Pa,
        "R_hyd": d.R_hyd,  # null when undefined (zero flow); see R_hyd_status
        "R_hyd_status": d.R_hyd_status,
        "k_min": d.k_min,
        "k_max": d.k_max,
        "phi_min": d.phi_min,
        "phi_max": d.phi_max,
        "p_min": d.p_min,
        "p_max": d.p_max,
        "q_scale_m3s": d.q_scale_m3s,
        "balance_limit_m3s": d.balance_limit_m3s,
        "open_area_m2": d.open_area_m2,
        "open_area_ideal_m2": d.open_area_ideal_m2,
        "dp_abs_Pa": case.sol.dp_abs_Pa,
        "dp_limit_Pa": case.sol.dp_limit_Pa,
        "dq_abs_m3s": case.sol.dq_abs_m3s,
        "local_abs_m3s": case.sol.local_abs_m3s,
        "params": asdict(params),
    }
    fields = {
        "p": case.sol.p, "k": case.sol.k, "phi": case.sol.phi,
        "sigma_eff": case.sol.sigma_eff, "Fr": case.sol.Fr, "Fz": case.sol.Fz,
        "r_centers": case.grid.r_centers, "z_centers": case.grid.z_centers,
    }
    def _write_fields(tmp: Path) -> None:
        # np.savez_compressed appends ".npz" to plain path names, so it gets
        # an open binary handle to honor the staged ".tmp" name exactly.
        with open(tmp, "wb") as fh:
            np.savez_compressed(fh, **fields)

    # No stale failure record may survive a success; summary goes last.
    _remove_owned(out_dir, ("failure.json",))
    _stage_publish(out_dir, [
        ("fields.npz", _write_fields),
        ("summary.json", lambda tmp: _write_json_strict(tmp, summary)),
    ])
    print(f"converged in {case.iters} iteration(s)")
    print(f"Q_out = {d.Q_out_m3s:.6e} m^3/s = {d.Q_out_mLs:.6f} mL/s")
    print(f"inlet/outlet mismatch = {d.mismatch_abs_m3s:.3e} m^3/s "
          f"(rel {d.mismatch_rel:.3e})")
    print(f"max cell residual = {d.max_residual_abs_m3s:.3e} m^3/s "
          f"(rel {d.max_residual_rel:.3e})")
    print(f"pressure drop = {d.delta_p_Pa:.6e} Pa, "
          f"R_hyd = {d.R_hyd} Pa s/m^3 [{d.R_hyd_status}]")
    print(f"k range = [{d.k_min:.3e}, {d.k_max:.3e}] m^2, "
          f"phi range = [{d.phi_min:.5f}, {d.phi_max:.5f}]")
    print(f"wrote {out_dir / 'summary.json'} and {out_dir / 'fields.npz'}")
    return 0


def cmd_sweep(args: argparse.Namespace) -> int:
    params = _case_params_no_pin(args)
    p_vals = parse_bar_list(args.p_in_bar)
    rows = run_sweep(params, p_vals)
    out_csv = Path(args.out)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=SWEEP_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))
    n_ok = sum(r.converged for r in rows)
    print(f"sweep: {n_ok}/{len(rows)} points converged -> {out_csv}")
    for r in rows:
        status = "ok " if r.converged else "FAIL"
        q = f"{r.Q_mLs:10.5f}" if r.converged else "       nan"
        print(f"  [{status}] p_in={r.p_in_bar:6.2f} bar  Q={q} mL/s  "
              f"iters={r.iters}")
    return 0 if n_ok == len(rows) else 2


def _case_params_no_pin(args: argparse.Namespace) -> SimParams:
    params = from_json_file(args.config)
    if args.set:
        params = apply_overrides(params, args.set)
    return params


def _tparams_from_args(args: argparse.Namespace) -> TransportParams:
    import json as _json

    with open(args.tconfig, "r", encoding="utf-8") as fh:
        tparams = tparams_from_dict(_json.load(fh))
    if args.tset:
        tparams = apply_toverrides(tparams, args.tset)
    return tparams


def cmd_transport(args: argparse.Namespace) -> int:
    from dataclasses import asdict as _asdict

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    hparams = None
    tparams = None
    try:
        hparams = _case_params_no_pin(args)
        tparams = _tparams_from_args(args)
    except (ValueError, KeyError, TypeError, OSError) as exc:
        h_dict = _asdict(hparams) if hparams is not None else None
        t_dict = None
        try:
            t_dict = _asdict(tparams) if tparams is not None else None
        except Exception:
            t_dict = None
        _publish_transport_invalid(out_dir, exc, args, h_dict, t_dict)
        print(f"INVALID INPUT: {exc}", file=sys.stderr)
        return 2
    assert hparams is not None and tparams is not None
    try:
        tc = run_transport(hparams, tparams)
    except (ValueError, KeyError) as exc:
        _publish_transport_invalid(
            out_dir, exc, args, _asdict(hparams), _asdict(tparams))
        print(f"INVALID INPUT: {exc}", file=sys.stderr)
        return 2
    if not tc.converged or tc.result is None:
        err = tc.error or ""
        if err.startswith("Hydraulics:"):
            err_type = "ConvergenceError"
        else:
            err_type = getattr(tc, "error_type", None) or "TransportError"
        _remove_owned(out_dir, ("transport_summary.json",
                                "transport_timeseries.csv",
                                "transport_fields.npz"))
        _stage_publish(out_dir, [(
            "transport_summary.json",
            lambda tmp, _tc=tc, _err_type=err_type: _write_json_strict(tmp, {
                "converged": False,
                "error": _tc.error,
                "error_type": _err_type,
                "n_steps": _tc.n_steps,
                "hydraulic_params": _asdict(hparams),
                "transport_params": _asdict(tparams),
                "hydraulic_config": getattr(args, "config", None),
                "transport_config": getattr(args, "tconfig", None),
                "hydraulic_overrides": getattr(args, "set", None),
                "transport_overrides": getattr(args, "tset", None),
                "assumptions": ASSUMPTIONS,
            }),
        )])
        print(f"TRANSPORT FAILED: {tc.error}", file=sys.stderr)
        return 1
    r = tc.result
    assert r is not None
    summary = {
        "converged": True,
        "mode": tparams.mode,
        "n_steps": r.n_steps,
        "t_end": tparams.t_end,
        "dt_requested": r.dt_requested,
        "dt_advective_limit": r.dt_advective_limit,
        "Q_out_m3s": tc.Q_out,
        "residence_time_tau_s": tc.tau_s,
        "M_solute_cup_kg": r.M_cup,
        "V_cup_m3": r.V_cup,
        "TDS_pct": r.TDS_pct,  # null when the cup is empty
        "TDS_status": r.TDS_status,
        "EY_pct": r.EY_pct,  # null in tracer mode (no coffee EY)
        "EY_status": r.EY_status,
        "e_mean": r.e_mean,
        "e_var": r.e_var,
        "budget_max_abs_kg": r.budget_max_abs,
        "budget_tol_kg": r.budget_tol,
        "lin_rel_max": r.lin_rel_max,
        "min_C_kgm3": float(np.min(r.C)),
        "min_B_kg": float(np.min(r.B)),
        "hydraulic_params": _asdict(hparams),
        "transport_params": _asdict(tparams),
        "assumptions": ASSUMPTIONS,
    }

    def _write_series(tmp: Path) -> None:
        with open(tmp, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=TRANSPORT_SERIES_FIELDS)
            writer.writeheader()
            for s in r.steps:
                writer.writerow({
                    "t_s": s.t, "dt_s": s.dt, "C_inlet_kgm3": s.C_in,
                    "Q_out_m3s": tc.Q_out, "M_solute_cup_kg": s.M_cup,
                    "V_cup_m3": s.V_cup,
                    "TDS_pct": (100.0 * s.M_cup
                                / (tparams.rho_solution * s.V_cup)
                                if s.V_cup > 0.0 else ""),
                    "EY_pct": (100.0 * s.M_cup / tparams.M_dose
                               if tparams.mode == "extraction" else ""),
                    "outlet_C_kgm3": s.outlet_C, "sum_B_kg": s.sum_B,
                    "sum_WC_kg": s.sum_WC,
                    "cumulative_inlet_kg": s.cumulative_inlet,
                    "ledger_signed_kg": s.ledger_signed,
                    "lin_rel": s.lin_rel, "min_C_kgm3": s.min_C,
                    "min_B_kg": s.min_B,
                })

    def _write_fields(tmp: Path) -> None:
        from .grid import make_grid as _make_grid

        _grid = _make_grid(hparams.Nr, hparams.Nz, hparams.R, hparams.L,
                           hparams.R_perf)
        with open(tmp, "wb") as fh:
            np.savez_compressed(
                fh, C=r.C, B=r.B, m=r.m, W=r.W, e_local=r.e_local,
                r_centers=_grid.r_centers, z_centers=_grid.z_centers,
            )

    _remove_owned(out_dir, ("transport_summary.json",))
    _stage_publish(out_dir, [
        ("transport_timeseries.csv", _write_series),
        ("transport_fields.npz", _write_fields),
        ("transport_summary.json",
         lambda tmp: _write_json_strict(tmp, summary)),
    ])
    print(f"transport converged in {r.n_steps} step(s) [{tparams.mode}]")
    print(f"cup: M={r.M_cup:.6e} kg V={r.V_cup:.6e} m^3 "
          f"TDS={r.TDS_pct} % [{r.TDS_status}] EY={r.EY_pct} % "
          f"[{r.EY_status}]")
    print(f"budget max |res| = {r.budget_max_abs:.3e} kg "
          f"(tol {r.budget_tol:.3e}); lin_rel max = {r.lin_rel_max:.3e}")
    print(f"wrote {out_dir / 'transport_summary.json'}, "
           f"{out_dir / 'transport_timeseries.csv'}, "
           f"{out_dir / 'transport_fields.npz'}")
    return 0


def cmd_shot_validate(args: argparse.Namespace) -> int:
    from .shots import ShotError, load_shot_file
    try:
        shot = load_shot_file(args.shot, args.series, args.aggregate)
    except (ShotError, FileNotFoundError, OSError, ValueError, KeyError) as exc:
        print(f"INVALID SHOT: {exc}", file=sys.stderr)
        return 2
    print(f"valid shot {shot.metadata.shot_id} [{shot.metadata.source}] "
          f"n={len(shot.samples)} dose={shot.metadata.dose_g:g} g "
          f"time_zero={shot.metadata.time_zero}")
    return 0


def cmd_shot_metrics(args: argparse.Namespace) -> int:
    from .shots import ShotError, load_shot_file, shot_summary
    try:
        shot = load_shot_file(args.shot, args.series, args.aggregate)
        summary = shot_summary(shot, args.flow_window)
    except (ShotError, FileNotFoundError, OSError, ValueError, KeyError) as exc:
        print(f"INVALID SHOT: {exc}", file=sys.stderr)
        return 2
    print(f"shot {summary['shot_id']} [{summary['source']}]")
    print(f"final mass = {summary['final_beverage_mass_g']:.3f} g "
          f"dose = {summary['dose_g']:.3f} g "
          f"brew ratio = {summary['brew_ratio']:.4f}")
    if summary["cup_EY_pct"] is None:
        print("cup EY: missing (no TDS; no inference made)")
    else:
        print(f"cup TDS = {summary['cup_TDS_pct']:.3f} % "
              f"[{summary['cup_TDS_method']}] "
              f"cup EY = {summary['cup_EY_pct']:.3f} % [measured]")
    print(f"derived flow window = {summary['derived_flow_info']['window_s']:g} s "
          f"[{summary['derived_flow_info']['method']}, "
          f"{summary['derived_flow_info']['endpoint_treatment']}] "
          f"neg_raw_increments = "
          f"{summary['derived_flow_info']['n_negative_raw_increments']}")
    if args.out is not None:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        _write_json_strict(out, summary)
        print(f"wrote {out}")
    return 0


def cmd_shot_plot(args: argparse.Namespace) -> int:
    from .shots import (ShotError, constant_flow_baseline, derive_flow,
                        load_shot_file, write_svg_plot)
    try:
        shot = load_shot_file(args.shot, args.series, args.aggregate)
        times = [s.time_s for s in shot.samples]
        masses = [s.beverage_mass_g for s in shot.samples]
        flows, _ = derive_flow(shot.samples, args.flow_window)
    except (ShotError, FileNotFoundError, OSError, ValueError, KeyError) as exc:
        print(f"INVALID SHOT: {exc}", file=sys.stderr)
        return 2
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    mass_series: list[tuple[str, list[float], list[float], str]] = [
        (f"{shot.metadata.shot_id} mass (observed)", times, masses, "observed"),
    ]
    if args.baseline_q is not None:
        if args.baseline_t0 is None or args.baseline_m0 is None:
            print("error: --baseline-q needs --baseline-t0 and --baseline-m0",
                  file=sys.stderr)
            return 2
        pred = constant_flow_baseline(times, args.baseline_q,
                                      args.baseline_t0, args.baseline_m0)
        mass_series.append(("constant-flow baseline (simulated)",
                            times, pred, "simulated"))
    write_svg_plot(str(out), mass_series,
                   xlabel="time_s (s)", ylabel="beverage_mass_g (g)",
                   title=f"{shot.metadata.shot_id} mass vs time")
    print(f"wrote mass-vs-time (time_s vs beverage_mass_g) -> {out} "
          "[observed=solid, simulated=dashed]")
    if args.flow_out is not None:
        fout = Path(args.flow_out)
        fout.parent.mkdir(parents=True, exist_ok=True)
        write_svg_plot(str(fout),
                       [(f"{shot.metadata.shot_id} derived flow (observed)",
                         times, flows, "observed")],
                       xlabel="time_s (s)", ylabel="derived_flow_g_s (g/s)",
                       title=f"{shot.metadata.shot_id} derived flow")
    print(f"wrote derived flow (time_s vs derived_flow_g_s) -> {fout} "
          f"[window {args.flow_window:g} s; derived, not measured]")
    return 0


def cmd_experiment_save(args: argparse.Namespace) -> int:
    from .experiment import save_experiment
    from .shots import ShotError
    try:
        info = save_experiment(args.shots, args.series,
                               args.already_wet_t0, args.holdout or [],
                               args.flow_window, args.out_dir,
                               args.label or "")
    except (ShotError, FileNotFoundError, OSError, ValueError, KeyError) as exc:
        print(f"INVALID EXPERIMENT: {exc}", file=sys.stderr)
        return 2
    comp = info["comparison"]
    print(f"saved experiment '{info['manifest']['label']}' -> {info['out_dir']}")
    print(f"train {info['manifest']['train_ids']} holdout "
          f"{info['manifest']['holdout_ids']} "
          f"q={comp['fit']['q_g_s']:.4f} g/s")
    for r in comp["heldout_residuals"]:
        print(f"held-out {r['shot_id']}: rmse={r['rmse_g']:.4f} g")
    print("provenance: measured/synthetic per shot; q fitted-empirical "
          "(train only); simulator constants illustrative")
    print("NOTE: " + comp["limitations"])
    return 0


def cmd_experiment_reproduce(args: argparse.Namespace) -> int:
    from pathlib import Path as _Path
    from .experiment import reproduce_experiment
    from .shots import ShotError
    try:
        report = reproduce_experiment(args.dir)
    except (ShotError, FileNotFoundError, OSError, ValueError, KeyError) as exc:
        print(f"INVALID BUNDLE: {exc}", file=sys.stderr)
        return 2
    out = _Path(args.dir) / "reproduce_report.json"
    _write_json_strict(out, report)
    status = "PASS" if report["passed"] else "FAIL"
    print(f"reproduce {status}: {args.dir}")
    for name, ok in report["checks"].items():
        print(f"  {name}: {'ok' if ok else 'MISMATCH'}")
    if report["hash_notes"]:
        for note in report["hash_notes"]:
            print(f"  {note}")
    print("NOTE: " + report["limitations"])
    return 0 if report["passed"] else 1


def cmd_shot_calibrate(args: argparse.Namespace) -> int:
    from .calibration import calibrate
    from .shots import ShotError, load_shot_file
    try:
        series_list = args.series or []
        if series_list and len(series_list) != len(args.shots):
            print("error: --series length must match --shots length",
                  file=sys.stderr)
            return 2
        shots = []
        for idx, spath in enumerate(args.shots):
            cpath = series_list[idx] if series_list else None
            shots.append(load_shot_file(spath, cpath, args.aggregate))
        hold = set(args.holdout or [])
        train = [s for s in shots if s.metadata.shot_id not in hold]
        held = [s for s in shots if s.metadata.shot_id in hold]
        if hold and not held:
            print(f"error: holdout ids {sorted(hold)} match no loaded shots",
                  file=sys.stderr)
            return 2
        if not hold:
            train, held = shots[:-1], shots[-1:]
            if len(shots) < 2:
                print("error: need >= 2 shots (or pass --holdout)",
                      file=sys.stderr)
                return 2
        result = calibrate(train, held, args.already_wet_t0, args.sim_q)
    except (ShotError, FileNotFoundError, OSError, ValueError, KeyError) as exc:
        print(f"INVALID CALIBRATION: {exc}", file=sys.stderr)
        return 2
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    _write_json_strict(out, result)
    print(f"fit q={result['fit_q_g_s']:.4f} g/s on {result['n_train']} train shot(s) "
          f"[fitted-empirical, train only; k NOT refitted]")
    if result["ci95_mean_g_s"] is not None:
        lo, hi = result["ci95_mean_g_s"]
        print(f"95% CI mean: [{lo:.4f}, {hi:.4f}] g/s; "
              f"95% PI new shot: [{result['pi95_new_shot_g_s'][0]:.4f}, "
              f"{result['pi95_new_shot_g_s'][1]:.4f}] g/s")
    for v in result["verdicts"]:
        extra = ""
        if v["verdict"] == "FAIL":
            extra = " -- PREDICTION FAILURE DISCLOSED"
        print(f"held-out {v['shot_id']}: slope={v['slope_g_s']:.4f} g/s "
              f"verdict={v['verdict']}{extra}")
    if result["simulator_comparison"] is not None:
        sc = result["simulator_comparison"]
        print(f"simulator Q={sc['sim_q_g_s']:.4f} g/s bias={sc['bias_empirical_minus_sim_g_s']:+.4f} "
              f"g/s inside-CI={sc['sim_inside_95ci']}")
    print("NOTE: " + result["real_data_note"])
    if result["status"] == "FAIL":
        return 1
    return 0


def cmd_shot_compare(args: argparse.Namespace) -> int:
    from .shots import ShotError, compare_shots, derive_flow, load_shot_file, write_svg_plot
    try:
        series_list = args.series or []
        if series_list and len(series_list) != len(args.shots):
            print("error: --series length must match --shots length",
                  file=sys.stderr)
            return 2
        shots = []
        for idx, spath in enumerate(args.shots):
            cpath = series_list[idx] if series_list else None
            shots.append(load_shot_file(spath, cpath, args.aggregate))
        hold = set(args.holdout or [])
        train = [s for s in shots if s.metadata.shot_id not in hold]
        held = [s for s in shots if s.metadata.shot_id in hold]
        if hold and not held:
            print(f"error: holdout ids {sorted(hold)} match no loaded shots",
                  file=sys.stderr)
            return 2
        if not hold:
            # No explicit holdout: hold out last whole shot for honesty.
            train, held = shots[:-1], shots[-1:]
            if len(shots) < 2:
                print("error: need >= 2 shots to compare (or pass --holdout)",
                      file=sys.stderr)
                return 2
        result = compare_shots(train, held, args.already_wet_t0)
    except (ShotError, FileNotFoundError, OSError, ValueError, KeyError) as exc:
        print(f"INVALID COMPARISON: {exc}", file=sys.stderr)
        return 2
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    _write_json_strict(out_dir / "comparison.json", result)
    # Mass overlay: observed traces + simulated constant-flow baseline.
    mass_series: list[tuple[str, list[float], list[float], str]] = []
    flow_series: list[tuple[str, list[float], list[float], str]] = []
    q, t0 = result["fit"]["q_g_s"], result["fit"]["t0_s"]
    from .shots import constant_flow_baseline as _base
    for s in shots:
        times = [p.time_s for p in s.samples]
        masses = [p.beverage_mass_g for p in s.samples]
        mass_series.append((f"{s.metadata.shot_id} mass (observed)",
                            times, masses, "observed"))
        after = [(p.time_s, p.beverage_mass_g) for p in s.samples if p.time_s >= t0]
        m0 = after[0][1] if after else masses[0]
        mass_series.append((f"{s.metadata.shot_id} baseline (simulated)",
                            times, _base(times, q, t0, m0), "simulated"))
        try:
            flows, _ = derive_flow(s.samples, args.flow_window)
            flow_series.append((f"{s.metadata.shot_id} derived flow (observed)",
                                times, flows, "observed"))
        except Exception:
            pass
    write_svg_plot(str(out_dir / "mass.svg"), mass_series,
                   xlabel="time_s (s)", ylabel="beverage_mass_g (g)",
                   title="repeated shots mass vs time")
    if flow_series:
        write_svg_plot(str(out_dir / "flow.svg"), flow_series,
                       xlabel="time_s (s)", ylabel="derived_flow_g_s (g/s)",
                       title="repeated shots derived flow")
    print(f"constant-flow fit on {len(train)} train shot(s): "
          f"q={q:.4f} g/s from t0={t0:g} s [simulated baseline; "
          f"permeability NOT varied to fit trace]")
    print(f"train repeatability: n={result['train_repeatability']['n_shots']} "
          f"final_mass mean={result['train_repeatability']['final_mass_g']['mean']:.3f} g "
          f"sd={result['train_repeatability']['final_mass_g']['stdev']:.3f} g")
    for r in result["heldout_residuals"]:
        print(f"held-out {r['shot_id']}: rmse={r['rmse_g']:.4f} g "
              f"max|res|={r['max_abs_resid_g']:.4f} g")
    print(f"wrote {out_dir / 'comparison.json'}, {out_dir / 'mass.svg'}"
          + (f", {out_dir / 'flow.svg'}" if flow_series else ""))
    print("NOTE: " + result["limitations"])
    return 0


def cmd_experiment_site(args: argparse.Namespace) -> int:
    from .shots import ShotError
    from .sitegen import build_site
    try:
        info = build_site(args.bundles, args.out, args.title)
    except (ShotError, FileNotFoundError, OSError, ValueError, KeyError) as exc:
        print(f"INVALID SITE: {exc}", file=sys.stderr)
        return 2
    print(f"wrote static viewer ({info['n_experiments']} experiment(s)) -> {info['out']}")
    print("open it directly in a browser; no backend, no live fitting")
    return 0


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="espresso-m1",
        description="Milestone-1 espresso packed-bed Darcy solver",
    )
    sub = ap.add_subparsers(dest="cmd", required=True)

    ap_run = sub.add_parser("run", help="run a single case")
    ap_run.add_argument("--config", required=True, help="JSON parameter file")
    ap_run.add_argument("--set", nargs="*", default=None,
                        help="overrides like Nr=40 p_in=9e5")
    ap_run.add_argument("--p-in-bar", type=float, default=None)
    ap_run.add_argument("--out-dir", default="outputs/case")
    ap_run.set_defaults(func=cmd_run)

    ap_sw = sub.add_parser("sweep", help="run a pressure sweep to CSV")
    ap_sw.add_argument("--config", required=True, help="JSON parameter file")
    ap_sw.add_argument("--set", nargs="*", default=None)
    ap_sw.add_argument("--p-in-bar", required=True,
                       help="'1,2,..,12' list or 'start:stop:step' in bar")
    ap_sw.add_argument("--out", required=True, help="output CSV path")
    ap_sw.set_defaults(func=cmd_sweep)

    ap_tr = sub.add_parser("transport", help="run transient transport")
    ap_tr.add_argument("--config", required=True, help="hydraulic JSON file")
    ap_tr.add_argument("--tconfig", required=True, help="transport JSON file")
    ap_tr.add_argument("--set", nargs="*", default=None,
                        help="hydraulic overrides like Nr=40")
    ap_tr.add_argument("--tset", nargs="*", default=None,
                        help="transport overrides like t_end=30")
    ap_tr.add_argument("--out-dir", default="outputs/transport")
    ap_tr.set_defaults(func=cmd_transport)

    ap_sv = sub.add_parser("shot-validate", help="validate a GS3 shot file")
    ap_sv.add_argument("--shot", required=True, help="shot JSON file")
    ap_sv.add_argument("--series", default=None, help="optional external CSV series")
    ap_sv.add_argument("--aggregate", default=None,
                       choices=["mean", "first", "last", "median"],
                       help="duplicate-timestamp rule (required if duplicates present)")
    ap_sv.set_defaults(func=cmd_shot_validate)

    ap_sm = sub.add_parser("shot-metrics", help="brew ratio + measured cup EY + flows")
    ap_sm.add_argument("--shot", required=True)
    ap_sm.add_argument("--series", default=None)
    ap_sm.add_argument("--aggregate", default=None,
                       choices=["mean", "first", "last", "median"])
    ap_sm.add_argument("--flow-window", type=float, default=1.0,
                       help="derived-flow window in seconds")
    ap_sm.add_argument("--out", default=None, help="optional strict-JSON summary path")
    ap_sm.set_defaults(func=cmd_shot_metrics)

    ap_sp = sub.add_parser("shot-plot", help="plot mass vs time and derived flow (SVG)")
    ap_sp.add_argument("--shot", required=True)
    ap_sp.add_argument("--series", default=None)
    ap_sp.add_argument("--aggregate", default=None,
                       choices=["mean", "first", "last", "median"])
    ap_sp.add_argument("--flow-window", type=float, default=1.0)
    ap_sp.add_argument("--out", required=True, help="mass-vs-time SVG path")
    ap_sp.add_argument("--flow-out", default=None, help="derived-flow SVG path")
    ap_sp.add_argument("--baseline-q", type=float, default=None,
                       help="optional constant-flow baseline q in g/s (simulated)")
    ap_sp.add_argument("--baseline-t0", type=float, default=None)
    ap_sp.add_argument("--baseline-m0", type=float, default=None)
    ap_sp.set_defaults(func=cmd_shot_plot)

    ap_sc = sub.add_parser("shot-compare", help="compare repeated shots (train/holdout)")
    ap_sc.add_argument("--shots", nargs="+", required=True, help="shot JSON files")
    ap_sc.add_argument("--series", nargs="*", default=None,
                       help="optional parallel CSV series (same order as --shots)")
    ap_sc.add_argument("--aggregate", default=None,
                       choices=["mean", "first", "last", "median"])
    ap_sc.add_argument("--holdout", nargs="*", default=[],
                       help="shot IDs held out whole (not used for fitting)")
    ap_sc.add_argument("--already-wet-t0", type=float, required=True,
                       help="already-wet interval start in seconds")
    ap_sc.add_argument("--flow-window", type=float, default=1.0)
    ap_sc.add_argument("--out-dir", required=True)
    ap_sc.set_defaults(func=cmd_shot_compare)

    ap_cal = sub.add_parser("shot-calibrate", help="M4 harness: fit q with uncertainty + holdout verdicts")
    ap_cal.add_argument("--shots", nargs="+", required=True)
    ap_cal.add_argument("--series", nargs="*", default=None)
    ap_cal.add_argument("--aggregate", default=None,
                       choices=["mean", "first", "last", "median"])
    ap_cal.add_argument("--holdout", nargs="*", default=[])
    ap_cal.add_argument("--already-wet-t0", type=float, required=True)
    ap_cal.add_argument("--sim-q", type=float, default=None,
                        help="optional simulator constant Q in g/s for bias report (not fitted)")
    ap_cal.add_argument("--out", required=True, help="strict-JSON calibration report path")
    ap_cal.set_defaults(func=cmd_shot_calibrate)

    ap_es = sub.add_parser("experiment-save", help="save reproducible experiment bundle")
    ap_es.add_argument("--shots", nargs="+", required=True)
    ap_es.add_argument("--series", nargs="*", default=None)
    ap_es.add_argument("--already-wet-t0", type=float, required=True)
    ap_es.add_argument("--holdout", nargs="*", default=[])
    ap_es.add_argument("--flow-window", type=float, default=1.0)
    ap_es.add_argument("--label", default="")
    ap_es.add_argument("--out-dir", required=True)
    ap_es.set_defaults(func=cmd_experiment_save)

    ap_er = sub.add_parser("experiment-reproduce", help="reproduce bundle from its inputs/")
    ap_er.add_argument("--dir", required=True)
    ap_er.set_defaults(func=cmd_experiment_reproduce)

    ap_site = sub.add_parser("experiment-site", help="build static HTML viewer (no backend)")
    ap_site.add_argument("--bundles", nargs="+", required=True)
    ap_site.add_argument("--out", required=True, help="output index.html path")
    ap_site.add_argument("--title", default="Espresso experiments")
    ap_site.set_defaults(func=cmd_experiment_site)
    return ap


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return int(args.func(args))
    except (ValueError, KeyError, TypeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
