# Task: repair and verify milestone 1 before adding physics

You are working in an existing Python/NumPy/SciPy espresso solver repository.
Implement the repairs below. Read `review/ASTRA_REVIEW.md` and
`review/submission_results.json` first. Preserve the reviewed snapshot and the
independent reviewer probe; write new run evidence to a different file. Existing
passing tests are a starting point, not evidence that the reported defects are
absent. Make focused changes; retain the numerical model and public run/sweep
entry points. This invocation authorizes completing all tasks below without
intermediate permission requests.

If a requirement conflicts with a conservation law or an executable example,
document the contradiction and continue independent work. Do not hide it by
changing an acceptance threshold or adding an artificial flow maximum.

## Task A — make convergence truthful

Relevant files: `solver.py`, `diagnostics.py`, `driver.py`, parameter/catalog
definitions as necessary.

Defect: the present residual evaluates relaxed k, not the physical k(p).
Acceptance is followed by another pressure solve and k update without a new
check. With Nr=8, Nz=16, compaction, R_perf=0.0145, relaxation=1e-8, it returns
success with about 20.6% global mass mismatch.

Implement this iteration contract:

1. Solve pressure using the current iterate of k. Iteration relaxation may
   influence the next iterate, but never the definition of the physical
   residual.
2. For the pressure candidate, evaluate **unrelaxed** k(p), phi(p), and stress.
3. Compute the candidate's actual face fluxes using that k(p). Compute local
   net outward flux in m^3/s and the global inlet-minus-outlet mismatch.
4. Accept only when pressure/flow update criteria AND physical local/global
   balance criteria pass. Return exactly this checked state. A later polish
   solve, if retained, must pass the same checks again.
5. Populate convergence diagnostics from the returned state and the actual
   consecutive candidates used in its acceptance. Never reuse stale metrics.

Use combined absolute and relative tolerances with documented units. Define
the local/global balance tolerance as

```
Q_scale = max(abs(Q_in), abs(Q_out))
balance_limit = flow_atol_m3s + tol_res_rel * Q_scale
```

Require each cell's absolute net flux and the absolute global mismatch to be
below this limit. Report dimensional residuals and dimensionless diagnostics;
the latter must have a documented scale. `flow_atol_m3s=1e-18` is an explicit
numerical default, not a physical parameter. Pressure-update relative scale
should depend on the imposed pressure difference, not the arbitrary pressure
offset; add a pressure absolute tolerance, e.g. 1e-6 Pa. Use physical Q for
flow-update comparisons, with the same absolute floor. Solving in pressure
relative to the outlet is a simple way to reduce cancellation.

Handle these exact equilibria explicitly or robustly: equal inlet/outlet
pressures give constant pressure and zero flux; R_perf=0 gives a sealed bed
with pressure equal to p_in and zero flux. Do not confuse a positive but
unresolved opening with an intentionally sealed outlet.

Acceptance:

- Run relaxation 1, 0.5, and 0.1 on the 8x16 perforated compaction case, with
  adequate iteration budgets. Tighten `tol_res_rel` to 1e-8 for these tests.
  Accepted solutions must have independently recomputed local/global balances
  below the configured limits and Q agreeing within 1e-6 relative.
- For relaxation 1e-8 and 1e-12, either converge to the same solution with
  valid balances or return explicit nonconvergence. Tiny iterate movement
  alone must never count as success. The default 100-iteration budget need
  not converge for these deliberately slow relaxations.
- Equal pressure 0/0 and 1e5/1e5 Pa must succeed with zero flow. A valid sealed
  outlet must also succeed. Keep nonconvergence reporting for a deliberately
  insufficient nonlinear iteration budget.
- Add a common pressure offset of 5e5 Pa to both boundaries of a forward-flow
  case; Q and k must remain invariant to numerical tolerance.
- Independent tests must recompute residuals from returned fields. Do not
  merely assert the solver's `converged` flag or its self-reported residual.

The old test asserting exactly three solves for uniform flow encodes an
implementation choice. It may be replaced with correctness and finite-cost
checks if the repaired solver uses a different number of solves; document why.

## Task B — validate inputs and resolve geometric intent

Relevant files: `params.py`, `catalog.json`, `grid.py`, CLI parsing.

- Reject NaN/infinity in every numeric input. Reject booleans where numerical
  scalars or integers are required. Reject fractional grid counts/iteration
  counts before coercion. Retain convenient CLI integer parsing.
- Align catalog constraints and enforced constraints. Distinguish mathematical
  admissibility from an illustrative recommended range; do not call an
  unenforced example range a hard valid range. Add actual source provenance
  for any literature value or label it an illustrative assumption.
- Retain the current cell-center mask for this repair; do not introduce a new
  mesh scheme. If R_perf>0 but open_mask has no true entries, reject with an
  actionable resolution error. Include requested and represented open areas
  in single-case and sweep metadata. Describe the solid-rim approximation.
- For zero flow, serialize undefined resistance as JSON null with a clear
  status or documented meaning, not NaN/Infinity. Require strict JSON writes.
- Parse ranges without overshooting their requested endpoint. `1:2:0.6` means
  [1, 1.6] bar. Allow a consistent descending range if supported; reject zero
  steps, inconsistent direction, nonfinite values, and empty sweeps. An empty
  sweep must not return successful completion. Ensure sensible negative-input
  rejection under the existing nonnegative gauge-pressure contract.

## Task C — keep artifacts and packaging trustworthy

Relevant files: `cli.py`, `pyproject.toml`.

- In a reused run directory, after a failed new run no old success summary or
  fields may masquerade as its output; after success no stale failure record
  may remain. Clean only artifacts owned by this command, not arbitrary user
  files. Stage file writes and publish the summary/status last. Carry the
  current configuration in the status even on failure where possible.
- Exercise success -> forced failure -> success in the same temporary path.
- Ensure `catalog.json` is present in the built wheel. Install a wheel into a
  temporary target and load its catalog while running outside the source tree
  with source PYTHONPATH removed. Editable/PYTHONPATH imports do not prove
  packaging works.
- Preserve per-point sweep failure records and a nonzero exit code when a
  sweep has failed points. Handle expected numerical solver failures as
  explicit failures; do not broadly catch programming errors as if they were
  ordinary nonconvergence.

## Task D — add meaningful verification and correct the report

Keep the original useful analytical/boundary tests. Add regression tests for
A–C. The review probe contains independent examples; turn the relevant
mathematical checks into maintained tests without changing the probe.

1. Nonlinear fully open reference: independently integrate the stated
   k(sigma), without calling the production closure to define the answer.
   At submitted defaults and 9 bar, the continuum Q is
   0.4687763667883215 mL/s. Mesh errors for 10x20, 20x40, 40x80 should decrease
   at approximately second order (observed error ratios about 4; require >3.5
   in this smooth case). A solver tolerance is not a discretization error.
2. Nontrivial 2D manufactured solution: use the field and exact integrated
   source in `review/probe_submission.py`. On 8x16, 16x32, 32x64, require
   volume-weighted L2 error reductions >3.8. Generate the source from the
   continuum differential equation, never from A @ sampled_exact_pressure.
   A test-only right-hand-side source is sufficient; no production source
   API is required.
3. Explain the Kirchhoff transform in `review/ASTRA_REVIEW.md`: for the current
   scalar, spatially uniform k(sigma), constant viscosity, fixed geometry,
   fixed top stress, and no sources, Q=(G/mu)*integral(k ds) in both fully open
   and perforated cases. Q is strictly increasing for positive k; with the
   current decreasing k it is concave. Normalized flow paths are independent
   of inlet pressure in the continuum. Do not assert an exact identity for
   the finite-grid harmonic-average nonlinear scheme; demonstrate approach
   to the identity with refinement.
4. Correct claims that staircase opening error explains the even-grid
   R_perf=R/2 refinement: those openings are exactly aligned. Report geometric
   error separately from discretization error and nonlinear residuals.
5. Withdraw unsupported tiny-opening results. Regenerate reports with full
   configurations and actual observed outcomes. Add a short reproducible
   validation script/command for the reported grids and sweeps; do not leave
   numerical tables that require guessing the command that generated them.
6. Describe this as a verified hydraulic prototype with illustrative
   constitutive parameters. It is not yet an experimentally validated
   extraction, deformation, or channel-formation simulator. A CG recommendation
   must specify an SPD-compatible preconditioner; generic ILU is not always
   suitable for CG.

## Completion

Run the maintained suite, the documented CLI example, and

```
PYTHONPATH=src ./.venv/bin/python review/probe_submission.py --wheel --out review/repaired_results.json
```

Complete the additional wheel-import and relaxation tests above as well. The
review probe is supplemental, not the entire acceptance specification. Compare
the new evidence with the original snapshot. Preserve its SHA-256 manifest
so reviewers can identify the submitted version. Do not overwrite original
review findings to make them describe the repaired implementation.

Deliver implemented fixes, regression tests, regenerated numerical evidence,
and updated docs. Report actual commands/results and any remaining limitations.
Finish this milestone before executing the separate transport prompt.
