# GS3 prompt 03 actual-results note (2026-09-23)

Equations and illustrative coefficients unchanged. Code verification only;
no measured coffee validation.

## A. Hardening (actual runs)

- Full suite: **109 passed** (`PYTHONPATH=src ./.venv/bin/python -m pytest tests/ -q`).
  Prior 91 retained without weakened tolerances; 18 new regressions added
  (`test_hardening_gs3.py`: 10, `test_shots.py`: 8).
- Independent hydraulic probe (no wheel build):
  **17/17** (`review/probe_submission.py --out <new path>`).
  Wheel packaging is covered separately by
  `test_wheel_contains_catalog_and_imports_off_tree` (passes as part of the 109).
- Invalid-input lifecycle (verified success → invalid → success):
  - `run`: malformed JSON (`JSONDecodeError`), physical validation
    (`R_perf=-0.01`, under-resolved opening, boolean `Nr`), solver
    nonconvergence (`max_iter=1`). Each invalid rerun removes owned stale
    `summary.json`/`fields.npz`, writes strict-JSON `failure.json` with
    `error_type` + current `config`/`set` context (`params: null` when no
    valid params; current params for convergence/under-resolved), preserves
    unrelated files (`notes.txt`/`keep.txt`), returns 2 for invalid vs 1 for
    nonconvergence. Resuccess removes `failure.json`.
  - `transport`: `t_end=-1`, missing `tconfig` (`FileNotFoundError`),
    hydraulic `max_iter=1`. Same contract with
    `transport_summary.json` + `error_type`, stale
    `transport_timeseries.csv`/`transport_fields.npz` removed, unrelated kept.
- Small horizons: `t_end=1e-13, dt=0.1` → 1 step, final `t==1e-13`;
  `t_end=1.0, dt=0.3` → 4 steps `[0.3,0.3,0.3,0.1]`, final `t==1.0`.
  Schedule `[0.7, 0.9]` with `dt=0.5` hits `0.5,0.7,0.9,1.0` exactly;
  no pulse skipped. Nonadvancing `t+dt==t` fails explicitly via
  `TransportError`.
- Boundaries: deliberately mixed-sign top (`-1e-7` vs `+` aggregate positive)
  rejected face-by-face despite positive `Q_in`; `assemble_transport` with
  `Fz[0,0]<0` raises `TransportError` (faulty negative-top branch removed).
  Sealed bed still freezes to `Q_in==0, Q_out==0`; forward case passes.
  Interior signed flows preserved (two-cell ±F tests pass).
- Normal demos unchanged: extraction 934 steps TDS 4.292341 %, EY 3.684558 %;
  tracer 516 steps TDS 0.25 %; final transport `t==60.0` exactly.

## B. Shot pipeline (synthetic only)

- `shot-validate` / `shot-metrics` / `shot-plot` / `shot-compare` run on
  `examples/shots/synthetic_*.json` (all `source: synthetic`).
- Linear `m=1.2*t` → derived flow 1.2 g/s at every sample (window 2 s);
  irregular sampling recovers 1.2 within 5 % interior; explicit missing TDS
  yields `EY null` with no inference; EY fixture `36 g * 8.5 % / 18 g`
  → `17.0 %` matches hand calculation.
- `shot-compare` (train linear, hold out whole irregular, `t0=0`):
  `q=1.2000 g/s`, held-out RMSE `0.0000 g`. Calibration remains DISABLED;
  comparison discloses constant-Q limitation and whole-shot holdout.

Missing measurements (not fitted): per-shot mass trace + `time_zero`, dose,
basket, grinder ordinal setting, bean batch/roast date, temperature setting +
location, confirmed pre-brew only, puck-located pressure only with
uncertainty (nominal 9 bar and AV meter ≠ measurements), TDS method +
uncertainty. See `SHOT_DATA.md`.
