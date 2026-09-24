# Independent review of the submitted milestone 1

Reviewed on 2026-09-23. This evaluates the repository and completion report
attributed to Muse by the user. Model identity and execution history are not
independently established. No submitted source code, configuration, or test
has been changed in this review. Repairs are assigned in
`prompts/01_REPAIR_M1.md`; the subsequent build is specified in
`prompts/02_TRANSPORT_AND_EXTRACTION.md`.

## Assessment

**Keep the implementation, repair the convergence and output contracts, and
retain independent review.** This submission demonstrates competent numerical
implementation and a useful willingness to report a negative scientific
result. It does not establish that the model can reliably validate its own
work. One sample is insufficient to rank Muse against frontier models in
general; no inference about latency, cost, or repeated-run reliability is
possible from the supplied material.

Positive evidence:

- All 22 original tests pass on this machine.
- The finite-volume geometry, indexing, harmonic interfaces, and pressure
  boundary signs appear consistent. An independent, genuinely two-dimensional
  manufactured solution converges at second order.
- Independent integration of the nonlinear constitutive law confirms the
  fully open solution converges toward Q = 0.4687763667883215 mL/s at 9 bar
  with the submitted default parameters.
- The implementation is modular and clearly labels uncalibrated compaction
  parameters. It did not manufacture a flow maximum to satisfy the original
  narrative.
- Its 1D monotonicity proof is correct. Identifying the impossibility of the
  requested turnover was a stronger scientific contribution than adding more
  simulation features would have been.

## Confirmed defects and limitations

| Priority | Finding | Reproduced evidence / consequence |
| --- | --- | --- |
| Blocker | False nonlinear convergence under relaxation | With compaction, Nr=8, Nz=16, R_perf=0.0145 m, relaxation=1e-8, the solver returns success in 3 solves while inlet/outlet mismatch is 0.20572265 and maximum cell imbalance / Q is 0.07979817. Correct default-relaxation Q is about 0.24676146 mL/s; the falsely accepted result is about 0.16253412 mL/s. |
| Blocker | Returned state is not the state whose convergence was checked | `solver.py:119` claims to evaluate k(p), but assembles with relaxed `k_new`. After acceptance, it solves once more, changes k again, and returns without checking the new residual. Reported convergence metrics belong to an earlier state. The final solve does not guarantee mass balance with the returned k. |
| High | An allowed equilibrium case falsely fails | The same 8x16 compaction case with p_in=p_out=1e5 Pa exhausts 100 iterations. The exact solution is constant pressure and zero flow; noise divided by a tiny flux floor prevents acceptance. |
| High | Unresolved outlet silently becomes sealed | R_perf=0.007 R has positive ideal area but no open columns on the default grid or the 8x16 probe grid. It is reported as a successful zero-flow problem. The report's tiny-opening probes cannot establish physical behavior without reporting resolution. R_perf=0 is a valid, distinct sealed-boundary case. |
| High | Failure leaves stale success artifacts | A successful CLI run followed by max_iter=1 in the same directory leaves `summary.json`, `fields.npz`, and `failure.json` together. A subsequent success leaves the old failure file. |
| Medium | Weak input validation | Nr=3.9 is silently truncated to 3; positive infinity is accepted for viscosity and pressure; boolean max_iter is coerced to 1. Catalog valid ranges and actual accepted values differ. |
| Medium | Sweep range lies about requested points | `1:2:0.6` produces 1, 1.6, 2.2 bar. `3:1:1` produces an empty CSV and a successful 0/0-points completion. |
| Medium | Installed wheel omits its parameter catalog | A wheel built from a temporary copy contains no `espresso_m1/catalog.json`; the installed `load_catalog()` contract is broken. |
| Medium | Zero-flow JSON is nonstandard | At p_in=p_out=0, the summary writes NaN for resistance. Strict JSON rejects it. Undefined physical quantities should have an explicit null/status representation. |
| Scientific reporting | Original verification was too narrow | Machine-precision agreement for a linear 1D profile does not establish 2D accuracy or nonlinear solver correctness. The independent tests below provide better evidence. |
| Scientific reporting | Explanation of mixed-boundary refinement is unsupported | For R_perf=R/2 and every even Nr in the report, the opening lies on a grid face and its area is exact. The observed first-order-like Q convergence cannot be attributed to stair-step area error in those runs. The Dirichlet/Neumann junction has reduced solution regularity; boundary flux treatment also warrants scrutiny. |

The convergence defect is not merely an extreme-value nuisance: relaxation is
an explicitly supported solver parameter. A small update does not prove a
small equation residual. This failure must be fixed before transient solute
transport relies on the hydraulic fluxes.

## Correction to Sol's prompt and a stronger mathematical result

Sol correctly recommended a staged implementation and refused to assume a
universal pressure optimum. However, it left the worker to discover that its
suggested class of closures cannot produce the intended downturn. This
feasibility check belonged in the planning stage.

The restriction applies to the current **2D perforated geometry as well as
the 1D open cylinder**, provided geometry, viscosity, top effective stress,
and the spatially uniform constitutive function are held fixed across the
pressure sweep. Let

```
s = sigma_top + p_in - p
H(s) = integral from sigma_top to s of k(a) da
I(DeltaP) = H(sigma_top + DeltaP)
u = grad(H)/mu
```

Continuity implies the axisymmetric Laplace equation for H. The top has H=0,
the open outlet has H=I, and the insulated boundaries have zero normal H
gradient. Write H=I*theta, where theta solves the geometry-only problem with
top 0 and open outlet 1. Then

```
G = integral over open outlet of grad(theta) dot n dA > 0
Q = (G/mu) * integral from sigma_top to sigma_top+DeltaP of k(a) da
dQ/dDeltaP = (G/mu) * k(sigma_top+DeltaP) > 0
```

G has units of length; for the fully open cylinder G=pi*R^2/L. If k decreases
with stress, Q is concave. With the implemented positive porosity floor, Q
eventually approaches a line with positive slope, not a finite plateau.
Normalized flow paths are pressure-independent in this continuum model;
pressure cannot create new preferred paths under these assumptions. Discrete
harmonic averaging does not preserve this transform exactly at finite mesh
spacing, so compare convergence toward the reference rather than asserting
machine-precision equality of two discretizations.

This is an independent mathematical deduction from the submitted equations.
It does not rule out downturns in experiments or models with different state
dependence. Adding mechanics, gas, or fines by name does not guarantee one.
The revised model must specify which assumption changes and demonstrate that
its equations can produce the proposed response without an imposed target.

The reviewed [Waszkiewicz et al. preprint, v2](https://arxiv.org/html/2512.21528v2)
also distinguishes transient brewing behavior from its long-time saturated
pressure response. Equations 7 and 16 give a nondecreasing response within
their constitutive validity range; extrapolating the polynomial beyond that
range can produce a misleading downturn. This is not grounds to dismiss the
paper's experimental observations. The journal DOI was inaccessible during
this review, so this statement is specifically about the accessible v2 text.

## Independent numerical evidence

Commands, run from the repository root:

```
PYTHONPATH=src ./.venv/bin/python -m pytest tests/ -q
PYTHONPATH=src ./.venv/bin/python review/probe_submission.py --wheel --out review/recheck_results.json
```

The original suite returned `22 passed`. The second command is intentionally
nonzero while the listed defects remain. `review/submission_results.json`
preserves the original review results, environment versions, and SHA-256
hashes of submitted files. Use a different output path for subsequent runs.
The wheel check builds in a temporary copy and may download build tooling.

For the nontrivial manufactured pressure field,

```
p(r,z) = p_in*(1-z/L) + amplitude*(1-(r/R)^2)^2*sin(pi*z/L)
```

the exact source is obtained by integrating the continuum cylindrical
Laplacian over each cell. It is not generated from the code's discrete
matrix. Volume-weighted L2 pressure error divided by amplitude:

| Grid | Error |
| --- | --- |
| 8x16 | 9.8671514e-4 |
| 16x32 | 2.4457782e-4 |
| 32x64 | 6.1041041e-5 |
| 64x128 | 1.5254301e-5 |

Independent nonlinear full-opening integral error falls by about four on
each refinement: 1.1423e-3, 2.8580e-4, 7.1464e-5, 1.7868e-5 for
10x20 through 80x160. For the half-radius outlet, the normalized nonlinear /
uniform flow ratio approaches the same constitutive integral, while absolute
Q converges more slowly. At 80x160 the ratio error is 1.1323e-4.

## Recommended next work and operating method

1. Give Muse `prompts/01_REPAIR_M1.md`. Keep the original evidence snapshot
   and use the independent probe to judge the repaired result. This is a
   bounded engineering assignment suitable for a competent smaller coding
   model, split into its listed tasks if necessary.
2. After those checks pass, give it `prompts/02_TRANSPORT_AND_EXTRACTION.md`.
   This adds conservative transient transport and a clearly labeled lumped
   extraction model on fixed hydraulics. It provides a testable path toward
   cup EY, TDS, and spatial extraction differences. It does not claim to solve
   pressure reversal, startup wetting, or grain diffusion.
3. Revisit coupled mechanics/thermal/gas/fines using observations that
   distinguish candidate mechanisms. Separate same-time, same-cup-mass,
   long-time equilibrium, and pressure-history comparisons; they are different
   experiments. Ask what existing assumption needs to change before selecting
   another subsystem.

The user's method is sound when "small-model task" means complete definitions,
bounded scope, objective checks, and explicit failure behavior. It becomes
counterproductive when the prompt fixes an unverified scientific conclusion
or too much implementation detail. Preserve an escape clause: the implementer
must identify a contradiction with a derivation or failing example, rather
than obeying a mistaken specification or silently redesigning it.

Evaluate the planner and implementer separately. Judge workers on previously
fixed acceptance checks, held-out perturbations, truthful failure reporting,
and total review/rework cost. A strong worker may use less effort on an explicit
task, but its own tests and completion prose are not independent verification.
For this submission, the worker discovered a mathematical limitation that the
planner missed, while the reviewer discovered an implementation defect the
worker's tests missed. That is useful evidence for retaining separate roles,
not for assigning an overall model league ranking.

## Addendum — Kirchhoff-transform note (repair implementation, 2026-09-23)

Appended during repairs; no finding above is altered. The repaired solver and
report use the following statement of the transform. With geometry, viscosity,
top effective stress, and the spatially uniform constitutive function fixed
across a sweep, set s = sigma_top + p_in - p and H(s) = integral of k. Then
u = grad(H)/mu, continuity gives the axisymmetric Laplace equation for H
(H = 0 top, H = I(DeltaP) open outlet, zero normal gradient insulated), and
with the geometry-only field theta (top 0, outlet 1),

    Q = (G/mu) * integral from sigma_top to sigma_top+DeltaP of k(a) da,

G > 0 a geometric factor (fully open cylinder: G = pi*R^2/L). Q is strictly
increasing for positive k and concave for decreasing k; the porosity floor
re-linearizes it at high pressure. Normalized flow paths are
pressure-independent in this continuum model. The finite-grid
harmonic-average nonlinear scheme does not preserve this identity exactly,
so the maintained tests (`tests/test_nonlinear_reference.py`) demonstrate
approach under refinement (fully-open normalized errors 1.14e-3, 2.86e-4,
7.15e-5) rather than exact equality of discretizations.
