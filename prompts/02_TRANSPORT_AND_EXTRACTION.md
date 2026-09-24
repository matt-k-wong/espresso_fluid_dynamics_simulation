# Task: conservative transient transport and lumped extraction on fixed hydraulics

Execute this after `prompts/01_REPAIR_M1.md` is complete and its acceptance
checks pass. This invocation authorizes the build below. If the hydraulic
repair checks still fail, finish that already specified repair first and
report its results separately before starting this phase. Do not treat the
original 22 passing tests as sufficient.

## Purpose and scientific limits

Add time-resolved liquid solute concentration, remaining soluble inventory,
cup solute mass, TDS, cup extraction yield, and spatial extraction differences
to the verified 2D axisymmetric hydraulic engine.

Use one hydraulically converged, saturated, isothermal snapshot per run. Freeze
its grid, porosity, permeability, and face volume fluxes for the entire
transport run. Time zero means an already wetted bed. The pressure-dependent
hydraulic law cannot produce a pressure-flow downturn under its present
assumptions; this phase does not alter that fact.

The extraction law below is a phenomenological exchange between one lumped
solid pool and the liquid. It is dimensionally defined and inventory
conserving. It does not resolve particle diffusion, swelling, moving solids,
CO2, fines, or temperature, and its parameters must remain illustrative until
calibrated. Those omissions must remain visible in the report. Do not insert
placeholder equations for those processes or fit a desired 6 bar result.

Preserve the existing steady run/sweep interfaces. Add the smallest coherent
transport module(s), tests, configuration, and a separate CLI subcommand such
as `transport`. Keep file output outside the numerical update functions.

If a requirement is mathematically contradictory, produce a derivation or
minimal counterexample and continue unaffected tasks. Do not silently replace
the physics or weaken the checks below.

## Task A — state, units, and conserved quantities

For every axisymmetric cell i define:

| Symbol | Meaning | Units |
| --- | --- | --- |
| V_i | Bulk cell volume from the existing grid | m^3 |
| phi_i | Frozen hydraulic snapshot porosity | dimensionless |
| W_i=phi_i V_i | Mobile liquid storage volume | m^3 |
| C_i | Solute mass / liquid solution volume | kg/m^3 |
| B_i | Remaining solute in the lumped solid pool | kg |
| m_i | Initial dry coffee mass assigned to the cell | kg |
| F_ij | Signed liquid volume rate from cell i toward j | m^3/s |
| D | Isotropic pore-liquid dispersion coefficient | m^2/s |
| lambda | Solid/liquid exchange rate coefficient | 1/s |
| K | Equilibrium solid solute mass fraction / liquid concentration | m^3/kg dry coffee |

Use `m_i = M_dose * V_i / sum(V)` as the explicit initial uniform dry-mass
distribution. Sum m_i to the requested dose exactly to floating-point
precision. Do not infer a new dry dose from the pressure-dependent phi field;
the current fixed-geometry closure is not a mass-conserving deformation model.
For extraction runs, initialize B_i=f_sol*m_i and C_i=0. The accessible soluble
fraction f_sol is an input, not a predicted coffee property. W_i remains fixed
as material dissolves: this is a passive transport approximation and ignores
the effect of mass removal on geometry and liquid properties.

For passive tracer runs use B_i=0, lambda=0, configurable initial C, and a
prescribed nonnegative inlet concentration history. Extraction runs use zero
inlet solute and zero initial dissolved solute so cup solute can be attributed
to the coffee. Do not report a coffee EY for externally injected tracer.

## Task B — spatial fluxes and boundary conditions

The semidiscrete equations are

```
W_i dC_i/dt + sum_outward(J_if) = R_i
dB_i/dt = -R_i
R_i = lambda * (B_i - K*m_i*C_i)         [kg/s]
```

This is reversible partition exchange. At local equilibrium B_i=K*m_i*C_i.
Liquid from another cell can return solute to a solid pool; a local B_i can
therefore exceed its initial value while global mass stays conserved. Do not
clip B to its initial value or call such redistribution a conservation bug.
K=0 is the irreversible first-order release limit. Default K>0 lets local
liquid concentration influence extraction; merely making release proportional
to B everywhere would impose identical solid depletion in a uniform bed even
when residence times differ.

Reuse the existing volumetric face fluxes. **They already contain face area.**
Do not multiply them by area again and do not divide them by phi when forming
solute transport. Pore velocity u/phi is relevant to a travel-time reference,
not to the advective flux F*C.

For an interior face oriented i -> j:

```
C_up = C_i if F_ij >= 0 else C_j
J_ij = F_ij*C_up + H_ij*(C_i-C_j)
H_ij = D * harmonic_mean(phi_i, phi_j) * A_face/distance
```

F and H both have units m^3/s. This convention discretizes
`d(phi*C)/dt + div(u*C - phi*D*grad(C)) = R/V`.
Use the same J with opposite signs in the two neighboring cells. Interior
radial F may be negative, so test both signs.

At the top prescribe **total inward solute flux** `F_in*C_in(t)` per face
(a Danckwerts-type inlet). This is one total-flux condition; do not also add
an independent Dirichlet diffusion flux there. At the perforated outlet use
outward advective flux `F_out*C_cell` and zero normal dispersive gradient.
Axis, sidewall, and solid bottom rim have zero total solute flux. For this
phase the hydraulic boundary flow must be forward or zero; reject an actual
reversed inlet/outlet flow with a clear unsupported-boundary message. Handle
roundoff consistently with the repaired hydraulic balance tolerance.

Keep no-flow synthetic states available to core tests for closed-batch and
diffusion checks. These are test inputs, not a new pressure boundary feature.

## Task C — use this conservative implicit update

Use backward Euler in time, first-order upwind advection, and centered
dispersion. Begin with a direct sparse solve. Construct L so `L*C - b_in`
equals the vector of net outward solute fluxes, including the outlet terms;
`b_in` contains top inward `F_in*C_in` in kg/s. L has units m^3/s.

The solid pool can be eliminated algebraically, leaving one linear system
for C per step. For dt>0 define

```
beta = dt*lambda / (1 + dt*lambda)

[diag(W + beta*K*m) + dt*L] C_new
    = W*C_old + beta*B_old + dt*b_in

B_new = (B_old + dt*lambda*K*m*C_new) / (1 + dt*lambda)
```

Array products on the right and in the diagonal are elementwise. These
equations follow by eliminating B_new from the two implicit balances; do
not separately add a second release source after this solve. With nonnegative
inputs, the intended discretization preserves nonnegative concentrations and
inventories. Check finite values, linear residual, positivity, and the global
mass budget after each step. Do not silently clamp negative values or change
inventories to make the budget balance. Report failed states explicitly.

Use the actual accepted dt for all updates. Shorten steps at input-history
breakpoints and the final time. A configurable conservative advective accuracy
limit such as `dt <= 0.5*min_i(W_i/sum_positive_outflow_i)` is useful even
though backward Euler is stable for larger dt; skip zero-outflow cells in
that estimate. Permit explicit dt choices in verification tests. Assess time
resolution by refinement rather than presenting implicit stability as
accuracy. No nonlinear transport solver or TVD limiter is needed in this phase.

## Task D — cup collection, extraction metrics, and artifacts

At each accepted step collect the **same numerical outlet flux** used in the
transport equation:

```
delta_M_solute_cup = dt * sum_outlet(F_out*C_new_adjacent)
delta_V_cup = dt * Q_out
M_solution_cup = rho_solution * V_cup
TDS_cup_percent = 100*M_solute_cup/M_solution_cup
EY_cup_percent = 100*M_solute_cup/M_dose
```

rho_solution is a constant illustrative solution density. Because Q is
solution volume flow in this approximation, do not add solute mass a second
time to rho_solution*V_cup. This is an approximate physical density closure,
not a prediction of composition-dependent solution density. For an empty cup,
TDS is undefined: serialize null/status. Detect physically invalid outputs
such as solute mass exceeding solution mass rather than clipping TDS.

Flux-weight the instantaneous outlet concentration by F_out. Do not average
outlet cell concentrations arithmetically. Report cup EY separately from
solid-pool depletion: solute already released but still inside W*C has not
reached the cup.

The global solute ledger is

```
sum(B) + sum(W*C) + M_solute_cup
    = sum(B_initial) + sum(W*C_initial) + M_solute_inlet_cumulative
```

Report its signed and absolute residuals. For the extraction example the
right-hand inlet contribution and initial liquid contribution are zero.

For spatial extraction, define local dry-dose extraction fraction
`e_i=(B_initial_i-B_i)/m_i`. Report its dry-mass-weighted mean and variance:
`mean=sum(m_i*e_i)/M_dose`,
`variance=sum(m_i*(e_i-mean)^2)/M_dose`.
Local e_i may be negative due to reabsorption in this reversible model;
explain this and preserve the value. The global cup EY remains bounded by
100*f_sol for extraction runs with zero initial/input liquid solute.

Export strict JSON metadata/status with full hydraulic and transport inputs,
a time-series CSV, and final fields including C, B, m, W, and local extraction.
Use the repaired output lifecycle so failed runs do not expose stale success
data. Record grid, dt policy, actual step count, budget residuals, and model
assumptions. No UI or plotting framework is required; a small ordinary plot
script is acceptable if it materially helps interpret the example.

Illustrative demo values (SI; label all as illustrative):

```
M_dose=0.018 kg; f_sol=0.25; lambda=0.1 1/s
K=0.005 m^3/kg; D=1e-9 m^2/s; rho_solution=1000 kg/m^3
t_end=60 s; initial C=0; inlet C=0
```

Use the repaired default compaction hydraulics and a half-radius perforated
outlet for one demo. Validate finite, nonnegative coefficients; require positive
dose/density/storage and positive time steps. K, D, and lambda may be zero for
the limiting tests below. State that no experimental calibration is supplied.

## Task E — independent acceptance checks

Tests must establish conservation and limiting behavior, not just compare a
function against another function that shares its implementation.

1. **Zero source/tracer:** zero initial and inlet solute with B=0 stays zero.
   For passive tracer with a uniform initial and inlet C, a divergence-free
   flow preserves that uniform C (up to hydraulic residual tolerance).
2. **Two-cell flux direction:** test positive and negative interior F and a
   concentration gradient with D>0. The receiving cell gains exactly the mass
   the donating cell loses when boundaries and sources are absent. Detect
   double multiplication by area using dimensional scaling.
3. **Closed-batch exchange:** with zero flow/dispersion, one cell has total
   M=B+W*C. The continuum reference is
   `C_eq=M/(W+K*m)`,
   `C(t)=C_eq+(C(0)-C_eq)*exp(-lambda*(1+K*m/W)*t)`,
   `B(t)=M-W*C(t)`.
   Demonstrate first-order time convergence (asymptotic error reduction near
   two for dt halving) and exact inventory conservation to roundoff. K=0 must
   approach B(t)=B0*exp(-lambda*t). lambda=0 must leave B unchanged.
4. **Passive residence time:** uniform fully open hydraulic case, D=0,
   lambda=0. Advective pore speed is u/phi, and transit time is
   `tau=sum(W)/Q`. Inject a known finite nonnegative pulse, align time steps
   with its boundaries, and verify input equals outlet-plus-stored tracer.
   Compare breakthrough/travel time with the continuum reference under grid
   and dt refinement. Upwind numerical spreading is expected; do not claim
   second-order accuracy for a sharp front.
5. **Independent dispersion check:** sealed synthetic flow, constant phi,
   lambda=0. A no-flux cylindrical diffusion eigenmode is
   `C=Cbar+A*J0(a*r/R)*cos(pi*z/L)*exp(-D*((a/R)^2+(pi/L)^2)*t)`,
   where a is the first positive zero of J1 (about 3.83170597). Choose Cbar>A>0.
   Compare at successively refined grids with time error controlled. Both
   radial and axial diffusion must contribute correctly. The mode is regular
   at the axis and has zero derivative at every boundary.
6. **Full extraction ledger and positivity:** run the perforated compaction
   example with exchange and dispersion. Require every step's global budget
   error <= `1e-12 kg + 1e-9*max(initial_total_solute, cumulative_input_solute)`.
   Report actual maxima and tolerances; budget acceptance does not establish
   spatial accuracy. Concentration/inventory must remain nonnegative to
   floating-point solver tolerance without corrective clipping.
7. **Metrics:** independently recompute flux-weighted outlet concentration,
   cup mass, cup TDS/EY, and dry-mass-weighted extraction variance from saved
   arrays/time series. Check cup EY <=100*f_sol within roundoff for extraction
   mode. Zero cup volume must produce valid strict JSON with undefined TDS.
8. **Accuracy evidence:** separate dt refinement at fixed grid from spatial
   refinement at a controlled smaller dt. Use at least three resolutions.
   Report changes in cup solute mass/EY and a field norm, and explicitly say
   whether the examples are adequately resolved. Do not tune K or lambda
   when changing the grid to make agreement appear better.

Build in this order: validated state and passive flux assembly; implicit
passive transport; exchange elimination and batch reference; cup ledger and
metrics; CLI/example and refinement evidence. Each subtask has an executable
check before the next one adds complexity. A stronger model may complete the
sequence in one session; a smaller model may be assigned one subtask at a time.

## Completion

Deliver working code, regression and analytical tests, one reproducible
extraction example, a passive tracer example, parameter provenance, and a
short numerical report. Show actual commands and results. Preserve all
milestone-1 acceptance checks. Report code verification separately from
physical validation: no measured shot data means the extraction parameters
and predictions remain unvalidated.

Finish with the data needed to distinguish later mechanisms: synchronized
pressure at the bed, beverage mass flow, temperature, time-resolved TDS,
initial dose/PSD, and an independently measured compression/permeability
relationship where available. Before proposing a pressure-downturn mechanism,
specify whether curves compare equal elapsed time, equal collected mass,
long-time equilibrium, or pressure history. Those experiments need not have
the same response.
