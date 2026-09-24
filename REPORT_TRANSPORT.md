# Milestone-2 report: transient transport and lumped extraction

A conservative transient solute-transport module on top of the verified
milestone-1 hydraulics. Each run freezes one converged saturated isothermal
snapshot (grid, porosity, permeability, face fluxes); time zero is an
already wetted bed. The extraction law is a phenomenological,
dimensionally defined, inventory-conserving exchange between one lumped
solid pool and the liquid. **Code verification below is separate from
physical validation: no measured shot data is used, so all transport and
exchange parameters and every cup prediction remain unvalidated.**

## Omitted physics (kept visible, no placeholders inserted)

Particle diffusion inside grains, swelling/moving solids, CO2/gas,
fines migration, temperature dependence, composition-dependent density.
The passive-storage approximation ignores mass removal effects on geometry
and liquid properties. The pressure-dependent hydraulic law cannot produce
a pressure-flow downturn under its present assumptions; this phase does not
alter that fact. No 6 bar result was fitted or targeted.

## Demos (actual runs; `scripts/reproduce.sh`)

Extraction (`transport --config configs/default.json
--tconfig configs/transport_extraction.json --set R_perf=0.0145`):
compaction hydraulics, half-radius perforation, illustrative
M_dose=0.018 kg, f_sol=0.25, λ=0.1 1/s, K=0.005 m³/kg, D=1e-9 m²/s,
t_end=60 s. Result: 934 steps, cup M=6.632205e-04 kg in V=1.545125e-05 m³,
**TDS=4.2923 %, EY=3.6846 %** (bound 25 %), e_mean=0.0667, e_var=0.00303.
Budget worst |res| = 4.2e-17 kg (tol 5.5e-12); linear residual max 8.2e-16;
C, B nonnegative with no clipping.

Tracer (`--config configs/uniform.json --tconfig configs/transport_tracer.json`):
fully open uniform hydraulics, D=0, 5 kg/m³ inlet pulse on [0,10)∪[20,40) s.
516 steps, cup M=2.972339e-04 kg, TDS=0.2500 %, budget 4.3e-19 kg.
No coffee EY is reported for injected tracer (`EY_pct` null).

## Verification evidence

| Check | Result |
|---|---|
| Zero tracer stays zero; uniform C preserved | exact 0; 5e-9 rel. (hydraulic residual scale) |
| Two-cell ±F with D>0 | hand-computed F+H entries; donor loss = receiver gain to solver roundoff |
| Closed batch vs analytic C(t), B(t) | dt ratios 2.7, 2.3 (first-order); inventory < 2e-15 rel. |
| K=0 → B0·e^{−λt}; λ=0 → B frozen | ratios ~2; exact |
| Reabsorption batch (C0-rich, B0=0) | B grows above initial, ledger exact, analytic match |
| Pulse: input = outlet + stored | to 1e-9 rel. |
| Outlet centroid vs τ+Tp/2 (τ=9.3333 s) | dt errs 2.2e-5, 1.0e-6, 1.0e-7; grid errs 3.3e-3, 4.2e-4, 1.7e-4 |
| J0 diffusion eigenmode L2/amplitude | 3.75e-4, 9.68e-5, 2.71e-5 on 8×16→32×64 (ratios 3.9, 3.6) |
| Extraction ledger every step | ≤ 1e-12 + 1e-9·max(initial, input); reported maxima above |
| Metrics recomputed from saved arrays | TDS/EY/e_mean/e_var match; EY ≤ 100·f_sol; empty cup → null TDS strict JSON |

## Resolution verdict

dt refinement at fixed 20×40 (t_end=20): M_cup = 2.14389692e-04,
2.14395231e-04, 2.14398003e-04 kg (successive ratios 2.0 — first-order in
time, as expected for backward Euler). Grid refinement at dt=0.0625:
M_cup = 2.0802e-04, 2.1440e-04, 2.1772e-04 kg (ratio 1.9 — first-order
upwind space); cup mass changes 1.5 % on the last refinement (≈1–2 %
remaining-error estimate at 40×80: adequate for cup metrics, not for
fine spatial detail). e_var = 5.16e-4, 6.56e-4, 7.49e-4 is still moving
(~14 % on the last step): the extraction-variance field metric is **not**
fully resolved. K and λ were identical across grids.

## Distinguishing future mechanisms (data still needed)

Synchronized bed pressure, beverage mass flow, temperature, time-resolved
TDS, initial dose/PSD, and an independently measured
compression/permeability relationship where available. Before any
pressure-downturn claim, specify whether curves compare equal elapsed time,
equal collected mass, long-time equilibrium, or pressure history.
