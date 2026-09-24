# Milestone 2 review, 2026-09-23

The new full suite reproduced: 91 tests passed in 9.55 seconds. The original
review script, run without its optional wheel build, passed 17/17. Results are
in `review/m2_hydraulics_recheck.json`. That run records source hashes for the
current submission; it does not replace the original M1 evidence.

Code inspection confirms the requested signed upwind/dispersion operator and
algebraic elimination of the solid pool. The cup accumulation uses the same
outlet flux as the transport equation. Independent checks reproduced two
additional edge cases:

1. **Stale success after invalid transport rerun.** Run a valid transport case
   into a temporary folder, then use the same folder with `--tset t_end=-1`.
   The CLI exits 2, but the old summary still says `converged: true`, and its
   fields and time series remain. The existing lifecycle test covers a
   hydraulic convergence failure, not validation failure before `run_transport`.
2. **Positive horizon skipped.**
   `run_transport(SimParams(Nr=4,Nz=8), TransportParams(t_end=1e-13,dt=.1))`
   returns success with zero steps. The while-loop tolerance is larger than
   the entire requested run. Very close schedule breakpoints use a separate
   absolute epsilon and deserve the same targeted check. This is a numerical
   contract defect, not a material source of error in the 60-second demo.

Also inspected: the top-boundary backflow branch in `assemble_transport` adds
F<0 to the diagonal even though outward advection there has rate -F. Production
transport is intended to reject boundary reversal, but `freeze_snapshot`
checks total rates rather than every boundary face. Normal supplied hydraulic
cases have forward flow. Remove the contradictory branch and enforce the
declared face-level domain, or correct and explicitly test the boundary
convention. Do not claim reversed-boundary transport is supported by accident.

These are bounded hardening tasks in prompt 03. They do not negate the normal
demo or analytical-reference tests. The reported resolution verdict is for
specific synthetic parameters/grids and should not become a universal accuracy
guarantee for the engine.

Comparison with the [Waszkiewicz/Lisicki paper](https://arxiv.org/html/2512.21528v2):
the repository resolves more spatial detail and numerical transport. The paper
has experimental evidence and an evolving-porosity hydraulic description that
this frozen-flow milestone lacks. More code and a 2D grid do not establish
greater physical accuracy. This review refers to the accessible v2 preprint;
it is not a claim that the final journal version was independently inspected.
