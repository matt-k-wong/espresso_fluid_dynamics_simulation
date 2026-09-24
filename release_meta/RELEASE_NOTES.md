# Release notes (candidate 0.1.0, uncalibrated, decision pending)

## What this is

2D axisymmetric single-phase saturated isothermal Darcy solver with a
pressure-dependent permeability closure; frozen-hydraulic transient
transport with lumped reversible exchange; GS3 shot-data baseline
(constant-flow comparison, whole-shot holdout); local experiment explorer
(save/reproduce bundle with provenance and exports).

## Actual equations, scope, and limits

- Hydraulics: `u = -(k/μ)∇p`, `∇·u = 0`; top Dirichlet `p_in`, mixed bottom
  (`p_out` on `r < R_perf`, no-flux rim), no-flux axis/sidewall; gravity off.
- Transport per frozen snapshot; `t = 0` is an already-wetted bed.
- Fixed bed; no wetting/gas, deformation, fines, heat, particles, or
  dissolution feedback. Closure cannot produce a pressure-flow downturn.
- Extraction variance less resolved than bulk cup metrics; all EY/TDS demo
  values illustrative, not predictions for a real machine.
- No optimal-pressure, flavor-prediction, full-digital-twin, or 3D-channel
  claims. Calibration against real coffee DISABLED (no measured data).

## Verification (not validation)

- `112 passed` full suite; independent hydraulic probe `17/17`
  (wheel check covered by packaging test).
- Demos: extraction 934 steps TDS 4.29 % / EY 3.68 %; tracer 516 steps.
- Fresh wheel-install smoke from this candidate: catalog load (34 params),
  `run`, `transport`, `shot-validate`, experiment save/reproduce PASS
  (see REPRODUCE_LOG.md).

## Known failures / unresolved

- Picard nonconvergence at extreme relaxations is reported, not hidden.
- Under-resolved perforations rejected as errors (use `R_perf = 0` sealed).
- M4 calibration and M5 extensions await real measurements; not included.
