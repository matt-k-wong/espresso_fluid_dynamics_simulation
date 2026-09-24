# Architecture note — Milestone 1

The numerical core (`src/espresso_m1/`) is independent of plotting and file
output so it can later be ported (e.g. to Rust). All I/O lives in `cli.py`.
No Rust implementation is included in this milestone. Data flow:

```
params.py ──▶ grid.py ──▶ closures.py ──▶ assembly.py ──▶ solver.py ──▶ diagnostics.py
   │                            │               │              │               ▲
   │                        driver.py ──────────┴──────────────┴───────────────┘
   │                            ▲
cli.py ─────────────────────────┴── (JSON configs, summary.json, fields.npz, sweep CSV)
```

## Modules

- **State (`solver.NonlinearResult`)**: cell fields `p`, `k`, `phi`,
  `sigma_eff` on `(Nr, Nz)` cell centers; face fluxes `Fr` on `(Nr+1, Nz)`
  radial faces (+r) and `Fz` on `(Nr, Nz+1)` axial faces (+z, downward);
  solve count plus dimensional indicators (`dp_abs_Pa`, `dq_abs_m3s`,
  `local_abs_m3s`, `mismatch_abs_m3s`) and their documented scales
  (`dp_limit_Pa`, `q_scale_m3s`, `balance_limit_m3s`).
- **Grid (`grid.Grid`, frozen dataclass)**: `dr`, `dz`, face positions,
  annulus column areas `π(re²−rw²)`, radial shell areas `2πr·dz`, cell volumes,
  `open_mask` (`r_center < R_perf`), discrete vs ideal open area. Mixed-bottom
  rule and solid-rim approximation documented in the module docstring; a
  positive `R_perf` with no open columns raises an actionable resolution
  error (only `R_perf = 0` is a sealed bed).
- **Coefficients (`closures.py`)**: `effective_stress`, `porosity_from_stress`
  (bounded, clipped), `kozeny_ratio`, `cell_permeability(p, params)`.
  Uniform branch returns `k0`/`phi0`; compaction branch guarantees `k > 0`.
- **Assembly (`assembly.py`)**: face transmissibility helpers
  (`radial_face_T`, `axial_face_T`, `top_face_T`, `bottom_face_T`), `assemble`
  → `(A_csr, b)` with `idx = j*Nr + i`, `compute_face_fluxes` → `(Fr, Fz)`
  using the same transmissibilities. Harmonic averaging rationale and the
  `spsolve` choice are documented in the module docstring.
- **Nonlinear solve (`solver.py`)**: `initial_pressure_guess` (linear in z),
  `solve_linear` (SuperLU, finite-check), `solve_pressure` (Picard under the
  contract in the module docstring: shifted `u = p − p_out` solve with the
  current k iterate, unrelaxed `k(p)` candidate fluxes, acceptance on
  pressure/flow updates AND physical balances, returned state is exactly the
  checked state). `ConvergenceError` on budget exhaustion carries the last
  indicators. `FLOW_ATOL_M3S = 1e-18` is a numerical floor, not physics.
- **Diagnostics (`diagnostics.py`)**: `evaluate` → `Diagnostics`: `Q_in`,
  `Q_out` (m³/s, mL/s), mismatch abs/rel, max cell residual abs/rel from
  explicit face-flux balances, `Δp`, `R_hyd = Δp/Q` (None if `Q ≤ 0`),
  min/max of `k`, `phi`, `p`.
- **Outputs (`driver.py`, `cli.py`)**: `run_case` (failure → `converged=False`
  record, never an exception to the CLI; only `ConvergenceError` counts as an
  ordinary failure), `run_sweep` (empty/non-negative-finite pressure lists
  rejected; per-point failure records with NaN fields, nonzero exit when any
  point fails), `analytical_cylindrical_Q`. `cli run` stages `fields.npz`
  then publishes `summary.json` last (atomic replace), removes stale owned
  artifacts on both paths, and touches nothing else in the directory;
  `cli sweep` writes the CSV with `p_in, Q, converged, iters`, open areas,
  and all diagnostics. `R_hyd` is JSON null with an explicit status when
  undefined; all JSON is written strictly (`allow_nan=False`). Pressure
  units: Pa in configs/code, bar accepted at the CLI (1 bar = 1e5 Pa).
- **Transport (`transport.py`, `tparams.py`)**: `TransportParams` + validation
  (extraction requires zero liquid solute; tracer requires `f_sol = 0`,
  `exchange_rate = 0`); `freeze_snapshot` (dry mass by bulk volume summing to
  the dose, forward/zero flow gate) and test-only `zero_flow_snapshot`;
  `assemble_transport` → `(L, top_inflow, outflow)` with signed upwind faces;
  `transport_step` (backward Euler, B elimination, linear residual out);
  `run_transport` (breakpoint-aligned dt under an advective accuracy limit,
  per-step finite/linear/positivity/budget checks, cup collection with the
  same outlet flux, `cup_metrics` with invalid-state detection);
  `TransportCase`/`TransportResult` carry series, final fields, and ledger
  maxima. CLI `transport` owns all file output under the same staged
  lifecycle (`transport_summary.json`, `transport_timeseries.csv`,
  `transport_fields.npz`).
- **Parameters (`params.py`, `tparams.py`, `catalog.json`, shipped in the
  wheel via `package-data`)**: `admissible_range` states enforced
  constraints; `recommended_range` is unenforced illustrative guidance;
  finite/boolean/fractional inputs are rejected before coercion. All
  transport entries carry illustrative-assumption provenance.
- **Shots (`shots.py`)**: versioned `gs3_shot_v1` importer/metrics/baseline
  (no machine I/O, raw retained, EY only with TDS, whole-shot holdout).
- **Experiment (`experiment.py`, M6)**: `save_experiment` copies inputs
  verbatim into `inputs/`, writes `comparison.json`, SVGs (observed solid vs
  simulated dashed, explicit units), `export/` CSVs, and `manifest.json`
  with sha256, train/holdout IDs, measured/synthetic/illustrative/
  fitted-empirical labels, limitations, and reproduce command;
  `reproduce_experiment` re-hashes and recomputes from bundle copies only.
  CLI owns all file output; core stays free of I/O assumptions.
- **Viewer (`sitegen.py`)**: `build_site` inlines bundle SVGs/tables into one
  static `index.html` (dropdown, escaped text, no network/JS backend).
- **Calibration (`calibration.py`, M4 harness)**: `calibrate` fits empirical
  `q` on train whole-shots with t-intervals, verdicts held-out whole shots,
  optionally reports bias vs an overlaid (never refitted) simulator `Q`;
  synthetic-verified, review-gated for measured data.

## Key invariants

- Permeability reaching the assembler is finite and strictly positive
  (asserted in `assemble` and `harmonic_face_k`).
- Axis/sidewall/rim fluxes are structurally zero (no matrix entries, zeroed
  flux arrays), verified by tests.
- Porosity stays in `[phi_min, phi0]` by construction (clip + `σ ≥ 0` floor).
- Gauge pressure used consistently; `p_out = 0` = atmosphere.
