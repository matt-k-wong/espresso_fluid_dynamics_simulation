# Milestone-1 report: pressure–flow behavior of the Darcy + compaction closure

This is a **verified hydraulic prototype with illustrative constitutive
parameters**. It is not an experimentally validated extraction, deformation,
or channel-formation simulator: it resolves steady saturated Darcy flow with
a scalar pressure-dependent permeability, and nothing else. All compaction
coefficients are illustrative assumptions (`catalog.json`); no pressure
optimum is predicted or claimed.

Every number below is produced by `scripts/reproduce.sh` (maintained suite,
documented CLI examples and sweeps, independent probe). Full run
configurations are stated; no table requires guessing its generating command.

## 1. Verification recap (repaired solver)

| Check | Result | Source |
|---|---|---|
| Analytical Q, uniform, fully open (20×40) | rel. error 3.8e-15 | `tests/test_analytical.py` |
| Refinement 10×20 → 40×80 | 3.0e-15 … 1.3e-14 (flat, solver noise) | same |
| Nonlinear fully-open vs independent quadrature (9 bar) | 0.468642 mL/s vs 0.46877637 continuum; normalized errors 1.14e-3, 2.86e-4, 7.15e-5 (ratios ≈ 4.0) | `tests/test_nonlinear_reference.py` |
| Nontrivial 2D manufactured solution, L2/amplitude | 9.87e-4, 2.45e-4, 6.10e-5 on 8×16, 16×32, 32×64 (ratios > 4.0) | `tests/test_mms.py` |
| Inlet/outlet mismatch, max cell residual (compaction ± perforation) | below `1e-18 + tol_res_rel·max(\|Q_in\|,\|Q_out\|)` by independent recompute | `tests/test_convergence_truthful.py` |
| Radial flux at axis/sidewall; bottom flux on solid rim | identically zero | `tests/test_boundaries.py` |
| Q(6)/Q(3), Q(9)/Q(3) with compaction off | 2 and 3 to 1e-9 (linear) | `tests/test_physics.py` |
| Equal pressures 0/0 and 1e5/1e5; sealed bed (R_perf=0) | converge, exactly zero flow | `tests/test_convergence_truthful.py` |
| +5e5 Pa common offset | Q and k invariant to ~1e-14 relative | same |
| Relaxation 1 / 0.5 / 0.1 (8×16 perforated, tol_res_rel=1e-8) | same Q = 0.24676146 mL/s within 1e-6, valid balances | same |
| Relaxation 1e-8 / 1e-12 | explicit nonconvergence, no silent success | same |

Compaction 9 bar, fully open: Q = 0.468 mL/s vs 1.982 mL/s uniform
(`cli run --config configs/default.json`, 18 Picard iterations).

## 2. Kirchhoff transform: why no flow maximum is possible here

For the milestone-1 closure, hold geometry, viscosity, top effective stress,
and the constitutive function fixed across a sweep and define the Kirchhoff
variable `H(s) = ∫ k(a) da` with `s = σ_top + p_in − p`. Then `u = ∇H/μ` and
continuity gives the axisymmetric Laplace equation for `H`, with `H = 0` at
the top, `H = I(Δp)` on the open outlet, and zero normal gradient on
insulated boundaries. Writing `H = I·θ` with the geometry-only field `θ`
(top 0, outlet 1) gives `Q = (G/μ)·∫_{σ_top}^{σ_top+Δp} k(a) da` with a
positive geometric factor `G` (fully open cylinder: `G = πR²/L`). Hence Q is
**strictly increasing** for positive k, concave for decreasing k, and
re-linearizes (never plateaus at finite slope loss) once the porosity floor
binds. Normalized flow paths are pressure-independent in this continuum
model. The finite-grid harmonic-average scheme does not preserve the
transform exactly, so the maintained tests demonstrate *approach* to the
identity under refinement (perforated normalized-ratio errors 3.09e-3,
1.17e-3, 3.74e-4 at 10×20 → 40×80) rather than asserting exact equality.

## 3. Pressure sweeps: no flow maximum found

Command: `cli sweep --config configs/default.json --p-in-bar "1:12:1"`
(plus `--set R_perf=0.0145` for the perforated sweep). With compaction on
(illustrative `σ_c = 2e5 Pa`, `σ_top = 1e4 Pa`, `φ0 = 0.35`, `φ_min = 0.20`),
Q(p_in) is monotonic increasing and concave: fully open 0.142 → 0.551 mL/s;
perforated 0.078 → 0.303 mL/s. Stronger compaction (`σ_c = 1e5 Pa`) and wider
probes (to 16 bar, `outputs/sweep_compaction_perf_strong.csv`) stay
monotonic. No parameters were tuned to force a peak and none produced one.

Earlier draft probes at very small perforations are **withdrawn**: a positive
`R_perf` that resolves to zero open columns is an under-resolved request,
not a sealed bed, and the repaired solver rejects it explicitly instead of
reporting a misleading zero flow. Only resolved configurations
(`open_area_m2`/`open_area_ideal_m2` are now recorded in every summary and
sweep row) support physical conclusions.

## 4. Which parameters control the rollover (no turnover)

Single-parameter sensitivities at 9 bar, perforated (`R_perf = 0.0145 m`):

- `σ_c`: 5e4 → 1e6 Pa gives Q = 0.159 → 0.601 mL/s (sets where concavity
  engages).
- `σ_eff_top`: 0 → 3e5 Pa gives Q = 0.268 → 0.153 mL/s (overall level via
  top-boundary mobility).
- `φ_min`: 0.10 → 0.30 gives Q = 0.099 → 0.681 mL/s (high-pressure asymptote
  and solver stiffness; Carman–Kozeny amplifies small changes).
- `k0`, `μ`, geometry enter as overall conductance scales.

## 5. Grid refinement, separated by error kind

Compaction, 9 bar: fully open Q = 0.46824 → 0.46864 → 0.46874 mL/s
(10×20 → 20×40 → 40×80); perforated (`R_perf = 0.0145 m`) Q = 0.2503 →
0.2575 → 0.2613 mL/s with the same monotonic concave sweep shape at 40×80.
Geometric and discretization errors are reported separately: for
`R_perf = R/2` with even `Nr` the junction lies exactly on a grid face and
the discrete open area is exact (relative area error 0 to 1e-16), so the
slower perforated convergence is attributed to the reduced regularity at
the Dirichlet/Neumann junction and the boundary flux treatment — not to
stair-step area error. The no-peak finding survives refinement.

## 6. Measurements needed to calibrate the compaction closure (before milestone 2)

1. **Controlled pressure–flow data**: Q vs `p_in` (1–12+ bar) on the same
   grind/dose/tamp, with shot-to-shot statistics — to test for any real
   turnover and fit the rollover scale (`σ_c`).
2. **Bed compression curve**: confined (oedometric) stress vs porosity for
   wet coffee grounds — to replace the exponential law and fix `φ0`,
   `φ_min`, `σ_c` (currently all illustrative).
3. **Permeability–porosity pairs**: measured `k` at several compacted states
   — to validate or replace Carman–Kozeny and anchor `k0`.
4. **In-situ pressure taps or outlet-flow distribution**: axial pressure
   profile and perforation-resolved outflow — to test the constant-total-stress
   assumption `σ_eff = σ_top + p_in − p` against a real stress field.
5. **Basket geometry + beverage viscosity**: actual `R`, `L(dose)`, `R_perf`,
   and liquid viscosity at brew temperature — to remove the remaining
   illustrative geometry/fluid values.

Until (1)–(3) exist, `σ_c`, `φ_min`, `σ_top`, and `k0` must stay labeled
illustrative and no pressure optimum may be inferred from this model.
