# Next assignment: harden M2, then prepare GS3 data comparison

Read ROADMAP.md and review/M2_REVIEW_AND_NEXT.md. Complete this assignment in
the espresso workspace. Keep the equations and illustrative coefficients
unchanged. The purpose is reliable evidence from the user's early GS3 AV on
reservoir water. Work can be done by a bounded coding agent; unresolved
identifiability or constitutive decisions belong in scientific review.

## A. Bounded correctness repairs

1. Treat a recognized run/transport invocation with invalid input as a failed
   new run in its requested output directory. Replace its status and remove
   only owned stale success artifacts; preserve unrelated files. Test
   success -> invalid parameter -> success, including malformed config,
   physical validation failure, and solver nonconvergence. Keep strict JSON.
   Store error type and the available current input/config context, not old
   parameters. Do not broadly catch programming bugs as ordinary failures.
2. A successful positive t_end must advance to the requested final time.
   Remove the absolute-epsilon zero-step escape. Either support representable
   small horizons or reject them explicitly as unsupported before solving.
   Test t_end=1e-13 and normal non-integer step counts. Step exactly to inlet
   schedule boundaries or use a documented representability rule; never skip
   a pulse silently. Detect a nonadvancing t+dt and fail explicitly.
3. Enforce forward-or-zero boundary flow face by face when freezing a snapshot.
   Recompute rate sums from the actual arrays. A deliberately mixed-sign
   boundary with positive aggregate flow must be rejected. Remove the faulty
   negative-top-flow branch or replace it with an explicit unsupported-state
   error. Preserve valid signed interior radial/axial flows. Check tolerance
   treatment does not introduce unaccounted mass changes.
4. Rerun the full suite and the independent hydraulic review probe. Add the
   targeted regressions without weakening existing tolerances. Write a short
   actual-results note.

## B. Versioned measurements and comparison tools

Create a small CSV/JSON schema and importer for shot observations, a synthetic
example labeled `synthetic`, and a documented capture template. Do not invent
measured GS3 data. Do not connect to or control the machine.

Shot metadata: shot ID, source type (measured/synthetic), machine profile ID,
AV/reservoir, bean batch/roast date if known, grinder and setting (an ordinal
setting, not particle diameter), dose_g, basket ID, beverage target_g, user
temperature setting with measurement location, pre-brew settings if confirmed,
and notes. No exact machine serial is needed in public examples.

Time series: time_s, beverage_mass_g, optional pressure_Pa and sensor location,
optional temperature_K and sensor location. Keep missing data missing. Record
whether time zero is pump activation, first liquid, or the modeled already-wet
state. Optional cup TDS_pct must carry measurement method and uncertainty where
known. Never infer EY without TDS or another actual solute measurement.

Reject nonfinite values, decreasing times, duplicate timestamps unless an
explicit aggregation rule is selected, negative dose, and inconsistent units.
Retain raw mass samples. A balance can have small negative mass increments due
to noise; flag them and use a documented optional smoothing method to derive
flow. Never rewrite the raw data to force monotonicity. Report the differentiation
window and endpoint treatment; don't call a noisy derivative measured flow.

Provide CLI commands to validate a shot, compute brew ratio and measured cup EY
when possible, plot mass versus time and derived flow, and compare repeated
shots. Plot explicit axes/units and distinguish observations from simulations.
For an initial baseline, compare a constant-flow mass trace over a declared
already-wet interval. Fit on training shots only; hold out whole shots, not
adjacent time rows from the same shot. Report repeatability and residuals.

The current simulator predicts constant Q in a transport run. Expose that
limitation in any comparison; do not secretly vary permeability to fit the
trace. If measured pressure at the puck is unavailable, effective conductance
may be estimated, but k and pressure losses cannot both be claimed identified.
Similarly, cup TDS alone does not uniquely determine f_sol, K, lambda, D and k0.

Accept synthetic fixtures only for importer/metric correctness: a known linear
mass trace with known derivative, an irregularly sampled trace, explicit
missing TDS, and a cup EY example with independently calculated units. Keep
calibration against real coffee disabled until real data is supplied.

## Completion

Deliver the repairs, schema, templates, importer, comparison command, tests,
and short instructions for the user's normal brewing workflow. Give a list of
missing measurements rather than fabricated fitted coefficients. Do not add
new fluid physics, a machine-control API, or modern GS3 capabilities without
evidence that the early machine supports them.
