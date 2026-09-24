# Espresso roadmap and definition of done

Updated 2026-09-24. Goal: a defensible, useful espresso experiment tool, with
an optional public research/education release. The original full multiphysics
vision remains a research direction, not a prerequisite for a useful v0.1.

## Current state

Milestone 1 hydraulics and its repair are implemented. Milestone 2 conservative
transport/lumped exchange is implemented. M2.1 status/time-boundary hardening,
the M3 GS3 shot-data pipeline, a minimal M6 experiment explorer
(save/reproduce bundle + static viewer), and the M4 calibration harness
(synthetic-verified; first measured-data run needs review) are implemented. On review, the full suite returned
**112 passed in ~10 s** and the unchanged hydraulic review probe returned
**17/17 without its optional wheel build**. The prior repaired evidence file
includes the wheel check. An M7 reviewable candidate directory is prepared
locally (`release_candidate/`, publication a separate owner decision).
These are code-verification results, not measured coffee validation.
See `review/M2_REVIEW_AND_NEXT.md` for remaining edge cases.

Known scientific limits: fixed saturated isothermal hydraulics per shot;
phenomenological exchange; no dissolution feedback on permeability or volume;
no resolved particles, startup wetting, gas, heat, fines, or deformation.
The current stress/permeability closure cannot predict a pressure-flow downturn.
The extraction-variance metric is less resolved than bulk cup metrics in the
reported study. Illustrative EY/TDS values are not predictions for a real machine.

## Milestones and gates

| Stage | Work | Definition of complete | Suggested agent |
| --- | --- | --- | --- |
| M1 | Conservative axisymmetric Darcy hydraulics and truthful nonlinear convergence | Existing analytical, geometric, nonlinear and artifact checks pass | Implemented; regression suite retained |
| M2 | Frozen-hydraulic conservative transport and finite soluble inventory | Analytical limits, solute ledger, positivity, refinement evidence | Implemented; edge-case hardening queued |
| M2.1 | Remaining status/time-boundary hardening | Invalid reruns cannot expose stale success; every successful positive horizon reaches its declared endpoint | Implemented; regressions retained |
| M3 | GS3 shot-data pipeline and honest comparison baseline | Versioned schemas, unit validation, synchronized data, plots, repeatability estimate, no synthetic data labeled measured | Implemented; prompt 03 |
| M4 | Small identifiable calibration and held-out validation | Separate fit/validation shots, uncertainty and baseline comparisons; failure to predict is explicitly reported | Strong scientific review; implementation may use Sol/Muse |
| M5 | Evidence-directed physical extension, if needed | Select one missing mechanism from residual patterns/data; conserve mass/energy and validate its isolated limit before coupling | Astra for model design/review; bounded coding tasks afterward |
| M6 | Usable local experiment explorer | Reproduce a saved shot/comparison, distinguish measured/fitted/illustrative inputs, export data and explain limitations | Implemented (minimal save/reproduce bundle) |
| M7 | Optional open-source release preparation | Fresh-install reproduction, provenance/license/privacy audit, release docs and reviewable candidate; publication is a separate owner decision | Candidate prepared locally; prompt 90; owner license/approval pending |

M6 can run after M3/M4 without waiting for every research extension. M7 can
follow a modest, honest v0.1; it need not wait for a universal espresso model.
Begin maintaining citations and provenance now, even if packaging waits.

## The user's machine

Confirmed by the user: **early-production La Marzocco GS3 AV, reservoir-fed**.
The user reports a first-year unit. Exact year, firmware, retrofit history,
group hardware, pump condition, basket, and grinder are not yet established.
Keep exact serial numbers and personal shot metadata in local private data;
they are unnecessary for public fixtures.

Start with a measured boundary model. Record beverage mass versus elapsed
time and known temperature/pre-brew settings. Record pressure only with its
sensor/gauge location and uncertainty. A coffee-boiler or pump gauge is not
automatically a measurement at the puck inlet. Do not invent a profile from
a nominal 9 bar setting or infer pump curves from unrestricted flow alone.
The AV volume-meter setting is not interchangeable with final cup mass.

Modern [La Marzocco guidance](https://home.lamarzoccousa.com/pre-brew-pre-infusion-and-pressure-manipulation-explained/)
distinguishes AV pre-brew controls from MP paddle pressure manipulation. It
does not establish which functions/retrofits are present on this early unit.
Reservoir operation supplies no mains-line pressure. Verify the exact machine
manual/settings before proposing pressure or pre-brew experiments.

## First useful home experiment

Use the current normal recipe as the control. Keep beans/age, dose, basket,
water, temperature setting, and puck preparation recorded. Repeat the control
on three occasions to estimate noise. Compare two nearby grind settings at a
fixed measured dose and target beverage mass, with three repetitions each,
interleaving order. Record elapsed time, mass trace where possible, and taste
notes. This nine-shot pilot can be spread across normal use; it is an example
design, not a statistical power guarantee. Add TDS only if a suitable instrument
and a consistent measurement protocol are available. Do not turn missing TDS
into a guessed extraction yield. Cup EY = beverage mass * TDS fraction / dose.

Do not optimize parameters against taste as if bulk EY uniquely predicted
flavor. Use the experiment to decide what to measure next and whether the
current simulator adds value beyond a simple empirical shot model.

## Completion levels

- **Learning/illustration prototype:** already substantially achieved, with
  status hardening pending. Honest plots and numerical references are useful.
- **Useful home decision tool:** M3, M4 and M6; must predict held-out shots or
  useful changes within measured uncertainty, and disclose failures.
- **Public research/education v0.1:** M2.1 plus M7 and reproducible examples;
  may explicitly remain uncalibrated.
- **Validated multiphase deforming-bed research engine:** several additional
  constitutive and experimental problems. There is no defensible percentage
  complete or calendar estimate from the current evidence.

The next scarce input is relevant measurements, not another large batch of
generated equations. Do not add all missing physics merely to lengthen the
feature list.
