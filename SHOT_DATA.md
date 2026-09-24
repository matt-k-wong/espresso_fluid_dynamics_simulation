# GS3 shot-data pipeline (M3 baseline, uncalibrated)

Versioned schema `gs3_shot_v1`. No measured GS3 data is included; every public
example has `"source": "synthetic"`. No code connects to or controls the machine.
All transport/exchange coefficients remain illustrative assumptions.

## Schema

Shot JSON (metadata + inline `samples`, or metadata alone + external CSV):

- `schema_version`: must be `"gs3_shot_v1"`.
- `shot_id`: non-empty string.
- `source`: `measured` | `synthetic`.
- `machine.profile_id`: string; `machine.feed`: string (e.g. `"AV reservoir"`).
  No exact serial in shared files.
- `beans.batch` / `beans.roast_date`: strings or null if unknown.
- `grinder.name`: string; `grinder.setting`: ordinal number/string
  (never particle diameter in µm).
- `dose_g`: grams, (0, 40] plausible; negative rejected.
- `basket_id`, `beverage_target_g` (>0 when present).
- `temperature.setting_C` + `temperature.location` (location required when
  set; a boiler/display setting is not a puck measurement).
- `prebrew`: object or null, only when confirmed on this early unit.
- `time_zero`: `pump_activation` | `first_liquid` | `already_wet`.
- `notes`: free text.
- `cup.TDS_pct` + `cup.method` + `cup.uncertainty_pct`: TDS in 0–30 %;
  missing stays null. EY is never inferred without TDS.

Time-series CSV header:

```
time_s,beverage_mass_g,pressure_Pa,pressure_location,temperature_K,temperature_location
```

Empty optional fields stay missing. Units are strict: `time_s` in seconds
0–600, `beverage_mass_g` in grams 0–1000, `pressure_Pa` in Pa 0–2e6
(use Pa, not bar), `temperature_K` in K 250–373 (use K, not °C).
Out-of-range values are rejected as inconsistent units.

Validation rejects: nonfinite values, decreasing times, duplicate timestamps
unless `--aggregate mean|first|last|median` is given, negative dose, and
inconsistent units. Raw mass samples are retained verbatim.

## Importer and CLI

```bash
PYTHONPATH=src ./.venv/bin/python -m espresso_m1.cli shot-validate \
  --shot examples/shots/synthetic_linear.json

PYTHONPATH=src ./.venv/bin/python -m espresso_m1.cli shot-validate \
  --shot examples/shots/synthetic_linear_meta.json \
  --series examples/shots/synthetic_linear_series.csv

PYTHONPATH=src ./.venv/bin/python -m espresso_m1.cli shot-metrics \
  --shot examples/shots/synthetic_ey.json --out /tmp/metrics.json

PYTHONPATH=src ./.venv/bin/python -m espresso_m1.cli shot-plot \
  --shot examples/shots/synthetic_linear.json \
  --out /tmp/mass.svg --flow-out /tmp/flow.svg \
  --baseline-q 1.2 --baseline-t0 0.0 --baseline-m0 0.0

PYTHONPATH=src ./.venv/bin/python -m espresso_m1.cli shot-compare \
  --shots examples/shots/synthetic_linear.json examples/shots/synthetic_irregular.json \
  --already-wet-t0 0.0 --out-dir /tmp/shot_cmp
```

- `shot-validate`: schema + ordering checks only.
- `shot-metrics`: brew ratio `beverage/dose`, measured cup EY
  `beverage*TDS/dose` only when TDS present (else `null` + `missing_no_TDS_no_inference`),
  derived flow by local linear regression in `--flow-window` seconds
  (method + window + one-sided endpoint treatment reported; noisy derivative
  is labeled derived, never measured). Negative raw increments are flagged,
  never removed.
- `shot-plot`: SVG mass-vs-time (`time_s` vs `beverage_mass_g`) and derived
  flow (`time_s` vs `derived_flow_g_s`) with explicit axes/units.
  Observed traces are solid; simulated constant-flow baselines are dashed.
- `shot-compare`: fits a constant-flow `q` on training shots only over the
  declared already-wet interval `[t0, end]`, holds out whole shots
  (`--holdout id...`; default holds out the last whole shot), reports
  repeatability (mean/sd of final mass, brew ratio, EY) and held-out
  residuals (RMSE, max|res|, mean). Writes strict-JSON `comparison.json`
  plus `mass.svg`/`flow.svg`.
- `experiment-save` / `experiment-reproduce` (M6 explorer): save a
  self-contained bundle from `inputs/` copies plus `comparison.json`,
  `mass.svg`/`flow.svg`, `export/combined_series.csv`,
  `export/metrics.csv`, `manifest.json` (file hashes, train/holdout IDs,
  per-shot measured/synthetic provenance, fitted-empirical `q` vs
  illustrative constants), and `README.txt` (limitations + reproduce
  command). `experiment-reproduce --dir` re-hashes inputs, recomputes from
  bundle copies, and PASS/FAILs on `q`, held-out RMSE, train mean, and
  hashes. Example:

```bash
PYTHONPATH=src ./.venv/bin/python -m espresso_m1.cli experiment-save \
  --shots examples/shots/synthetic_linear.json examples/shots/synthetic_irregular.json \
  --already-wet-t0 0.0 --label demo --out-dir /tmp/exp_demo
PYTHONPATH=src ./.venv/bin/python -m espresso_m1.cli experiment-reproduce \
  --dir /tmp/exp_demo
```

- `experiment-site` (static viewer, no backend): builds one self-contained
  `index.html` from saved bundles — experiment dropdown, inline mass/flow
  SVGs, metrics and held-out-residual tables, per-input provenance, and the
  limitations/calibration text. Open the file directly in a browser; no
  server, no live fitting, no machine connection:

```bash
PYTHONPATH=src ./.venv/bin/python -m espresso_m1.cli experiment-site \
  --bundles /tmp/exp_demo --out /tmp/exp_site/index.html --title "Demo"
open /tmp/exp_site/index.html
```

## Baseline and non-identifiability

The simulator predicts constant `Q` per transport run (frozen hydraulics,
already-wet bed). Comparison overlays `m(t)=m0+q*(t-t0)` and reports
residuals; permeability is never secretly varied to fit a trace.

- Without puck-pressure measurement, effective conductance may be estimated
  but `k` and pressure losses cannot both be claimed identified.
- Cup TDS alone does not determine `f_sol, K, lambda, D, k0`.
- `CALIBRATION_STATUS`: calibration against real coffee is DISABLED until
  real measured shot data is supplied. Synthetic fixtures test
  importer/metric correctness only.
- `shot-calibrate` (M4 harness, synthetic-verified): fits the empirical
  flow `q` on train whole-shots with a 95 % t confidence interval and a
  95 % prediction interval, verdicts each held-out whole shot
  (`pass` / `FAIL` / `unassessed_single_train`), and optionally reports
  empirical-vs-simulator bias for a `--sim-q` that is overlaid, never
  refitted. Prediction failures set exit code 1 and are disclosed in the
  strict-JSON report. Verified on synthetic only; the first run with
  `source: measured` data requires scientific review:

```bash
PYTHONPATH=src ./.venv/bin/python -m espresso_m1.cli shot-calibrate \
  --shots examples/shots/synthetic_linear.json examples/shots/synthetic_irregular.json \
  --already-wet-t0 0.0 --sim-q 1.2 --out /tmp/cal.json
```

## Normal brewing workflow (nine-shot pilot, spread across normal use)

1. Keep beans/age, dose, basket, water, temperature setting, puck prep fixed.
   Record grinder name + ordinal setting (not µm), dose_g, basket ID,
   beverage target_g, temperature setting + location, pre-brew only if
   confirmed, and `time_zero` convention. Use `templates/shot_template.json`
   + `templates/series_template.csv`.
2. Control: repeat the current normal recipe on three occasions (mass trace
   + time + taste notes). Log beverage mass vs elapsed time from the scale
     (do not edit raw dips).
3. Two nearby grind settings at fixed dose/target, three repetitions each,
   interleaved. Record elapsed time, mass trace, taste notes.
4. Add TDS only with a consistent instrument/protocol (method + uncertainty);
   never turn missing TDS into guessed EY.
5. Run `shot-validate`, then `shot-metrics` (brew ratio, EY when TDS present),
   then `shot-compare` with whole-shot holdout. If residuals show systematic
   curvature beyond constant flow, that is evidence for a missing mechanism
   to discuss in scientific review — not a license to secretly vary `k`.

## Missing measurements (needed before any fitting)

- Beverage mass vs elapsed time with documented `time_zero` and scale resolution.
- Dose_g per shot (measured, not nominal) and basket ID.
- Grinder + ordinal setting per shot; bean batch/roast date where known.
- Temperature setting + measurement location; pre-brew settings only if
  confirmed present on this early AV reservoir unit.
- Pressure only with sensor/gauge location + uncertainty (boiler/pump gauge
  is not puck pressure; nominal 9 bar is not a measurement).
- AV volume-meter setting is not interchangeable with final cup mass.
- Cup TDS with method + uncertainty where available; time-resolved TDS if ever.
- Exact serials stay in local private data, never in shared fixtures.

No new fluid physics, machine-control API, or modern GS3 capabilities are
added here without evidence the early unit supports them.
