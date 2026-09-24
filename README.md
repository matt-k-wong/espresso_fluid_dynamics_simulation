# Espresso packed-bed simulation engine — Milestone 1

2D axisymmetric, single-phase, saturated, isothermal Darcy-flow solver with a
configurable pressure-dependent permeability closure.

Scope: **milestone 1 only** — steady liquid flow in a fixed bed. No wetting /
gas dynamics, no bed deformation, no fines transport, no heat or extraction
(those are later milestones requiring separate approval and validation, and
must not be smuggled in as placeholder terms).

## Model equations

Coordinates `r ∈ [0, R]`, `z ∈ [0, L]`, `z` increasing downward. Superficial
Darcy velocity `u`, gauge pressure `p` (Pa, `p = 0` = atmosphere):

```
u = -(k/μ) ∇p        (Darcy, gravity neglected)
∇·u = 0              (incompressible, quasi-steady, saturated)
```

Boundary conditions:

- Top (`z = 0`): `p = p_in` (Dirichlet) over the full radius.
- Bottom (`z = L`): mixed — `p = p_out` (Dirichlet) on the perforated part
  (`r < R_perf`), zero normal flux on the solid rim (`r ≥ R_perf`).
- Axis (`r = 0`), sidewall (`r = R`): zero normal flux.

Permeability modes (`permeability_mode`):

1. `uniform`: `k = k0` (verification).
2. `compaction`: bounded porosity–effective-stress law + Carman–Kozeny,
   `σ_eff = σ_eff_top + p_in − p` (prescribed top load, constant total axial
   stress — a reduced hydraulic closure, not a deformation model),
   `φ = φ_min + (φ0 − φ_min)·exp(−σ_eff/σ_c)`,
   `k = k0·[φ³/(1−φ)²]/[φ0³/(1−φ0)²]`, with `φ ∈ [φ_min, φ0]`, `k > 0`.

See `src/espresso_m1/closures.py` for the load assumption and bounds, and
`REPORT.md` for what the closure can and cannot produce.

## Assumptions and known limits

- Saturated, single-phase, isothermal, quasi-steady; no gravity head.
- Fixed bed geometry (no deformation coupling in milestone 1).
- Poroelastic coupling is a scalar `k(p)` closure; it does not resolve
  channels, fines migration, or early multiphase wetting (cf. Waszkiewicz et
  al., Phys. Fluids 38, 063113 (2026), who leave these unresolved; Cameron et
  al., Matter 2, 631–648 (2020), is used as context for flow inhomogeneity,
  not as proof of any universal pressure threshold).
- The model does **not** predict a universal pressure optimum. With the
  milestone-1 closure no flow maximum was found (see `REPORT.md`).
- All compaction coefficients are **illustrative assumptions**, not measured
  coffee properties (see `src/espresso_m1/catalog.json` provenance labels).

## Units

SI everywhere: m, Pa (gauge), Pa·s, m², m³/s (plus mL/s as a convenience
column in outputs). The CLI accepts bar for inlet-pressure specifications and
converts with 1 bar = 1e5 Pa.

## Solver strategy

Cell-centered finite volumes with axisymmetric areas/volumes (2πr factor);
interior faces use harmonic averaging of `k/μ`; Dirichlet boundaries use
half-face transmissibilities; the mixed bottom is a stair-step split by cell
center (`r_center < R_perf` open, solid rim otherwise). A positive `R_perf`
that resolves to zero open columns is rejected as under-resolved (only
`R_perf = 0` selects the sealed bed). The SPD sparse system is solved with
`scipy.sparse.linalg.spsolve` (direct SuperLU — robust and deterministic at
milestone-1 sizes; for larger problems CG would need an SPD-compatible
preconditioner such as AMG or incomplete Cholesky, not generic ILU). The
linear system is solved for `u = p − p_out` and shifted back, so a common
pressure offset cannot affect acceptance. Nonlinear `k(p)` coupling uses
Picard iteration: each pressure candidate is checked with the **unrelaxed**
`k(p)` (relaxation only shapes the next iterate), and acceptance requires
the pressure update (`tol_p_atol_Pa + tol_p_rel·|p_in − p_out|`), the flow
update (absolute floor `1e-18 m³/s` + `tol_q_rel`), AND the physical
local/global balances (`1e-18 m³/s + tol_res_rel·max(|Q_in|, |Q_out|)`).
Failure raises `ConvergenceError` and is recorded, never silently returned.
Details: `ARCHITECTURE.md`.

## Run commands

Setup (once): `python3 -m venv .venv && ./.venv/bin/pip install numpy scipy pytest`

- One case: `PYTHONPATH=src ./.venv/bin/python -m espresso_m1.cli run --config configs/default.json --out-dir outputs/case`
- Override params: `... run --config configs/default.json --set R_perf=0.0145 Nr=40 Nz=80 --out-dir outputs/x`
- Pressure sweep to CSV: `PYTHONPATH=src ./.venv/bin/python -m espresso_m1.cli sweep --config configs/default.json --p-in-bar "1:12:1" --out outputs/sweep.csv`
  (ranges never overshoot: `1:2:0.6` means 1, 1.6 bar; `3:1:-1` descends)
- Tests: `PYTHONPATH=src ./.venv/bin/python -m pytest tests/ -q`
- Full reproduction (suite + examples + sweeps + independent probe):
  `sh scripts/reproduce.sh`

`R_hyd` is JSON `null` with status `undefined_zero_flow` when flow is zero;
all JSON output is strict (no NaN/Infinity).

## Milestone 2 — transient transport and lumped extraction

Time-resolved solute on one frozen hydraulic snapshot per run (grid,
porosity, permeability, face fluxes fixed; t = 0 is an already wetted bed):

```
W_i dC_i/dt + Σ_outward J_if = R_i,   dB_i/dt = −R_i
R_i = λ (B_i − K m_i C_i)             [kg/s, reversible partition]
J = F·C_up + H·(C_i − C_j),  H = D·harm(φ)·A/d   [m³/s, F contains area]
```

Backward Euler + first-order upwind + centered dispersion, one sparse direct
solve per step after exact elimination of the solid pool; per-step checks of
finite values, linear residual, positivity (never clipped), and the global
ledger `ΣB + Σ(WC) + M_cup`. Cup: flux-weighted outlet concentration,
`TDS = 100·M_cup/(ρ·V_cup)`, `EY = 100·M_cup/M_dose` (extraction only — no
coffee EY for injected tracer; TDS null when the cup is empty). Omitted and
kept visible: grain diffusion, swelling/moving solids, CO₂, fines,
temperature, composition-dependent density. All transport coefficients are
illustrative and unvalidated — see `REPORT_TRANSPORT.md`.

- Extraction demo: `... cli transport --config configs/default.json --tconfig configs/transport_extraction.json --set R_perf=0.0145 --out-dir outputs/transport_extraction` → TDS 4.29 %, EY 3.68 % at 60 s.
- Tracer demo: `... cli transport --config configs/uniform.json --tconfig configs/transport_tracer.json --out-dir outputs/transport_tracer`
- Full reproduction: `sh scripts/reproduce.sh`

## Repository layout

```
src/espresso_m1/   numerical core (no plotting, no file I/O)
  params.py, catalog.json   parameter handling + machine-readable catalog
  grid.py, closures.py      grid/geometry, permeability closures
  assembly.py               FV transmissibilities, sparse assembly, face fluxes
  solver.py                 linear + Picard nonlinear solver
  diagnostics.py, driver.py flow/balance diagnostics, case + sweep drivers
  cli.py                    CLI (owns all file output)
configs/           default.json (compaction example), uniform.json
tests/             analytical, conservation, boundary, linearity, catalog tests
outputs/           generated summaries, fields, sweep CSVs (not tracked)
README.md, ARCHITECTURE.md, REPORT.md
```

## Verified results (actual runs, 20×40 unless noted)

- Analytical check (uniform, fully open): relative error 3.8e-15 at 20×40;
  3.0e-15 … 1.3e-14 across 10×20 … 40×80 (error flat at solver noise).
- Nonlinear fully-open flow vs independent quadrature: 0.468642 mL/s against
  0.46877637 mL/s continuum; normalized errors shrink ≈4× per refinement.
- Nontrivial 2D manufactured solution: L2 error/amplitude 9.87e-4, 2.45e-4,
  6.10e-5 on 8×16, 16×32, 32×64 (second order).
- Conservation: every accepted state passes independently recomputed local
  and global balances below `1e-18 + tol_res_rel·max(|Q_in|, |Q_out|)`; tiny
  relaxations (1e-8, 1e-12) report explicit nonconvergence instead.
- Boundaries: radial flux zero at axis/sidewall; zero bottom flux on the rim;
  equal-pressure and sealed equilibria converge with exactly zero flow.
- Linearity (uniform): Q(6 bar)/Q(3 bar) = 2, Q(9)/Q(3) = 3 to 1e-9.
- Compaction sweep (illustrative params): Q(p) monotonic and concave, no peak
  at any tested pressure (1–16 bar), perforation, or compaction strength;
  shape unchanged under grid refinement. See `REPORT.md`.
