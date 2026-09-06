# Methodology and Development History

This document explains the reasoning behind the solver development, including
unsuccessful approaches. It is a record of why the current KT-D1V5,
finite-volume, generalized-minmod, SSP-RK3 combination was selected.

## 1. Define a pointed application

The broad topic was compressible flow with the discrete Boltzmann method. That
was narrowed to a staged objective: validate a compressible FVDBM with shock
tubes, then apply it to supersonic wedge flow and eventually examine
thermodynamic nonequilibrium near oblique shocks.

A wedge was selected before an airfoil because its oblique-shock angle has a
strong analytical benchmark and its geometry is simpler. The airfoil remains
a possible later application, not a current claim.

## 2. Start in one dimension

The first implementation used the Sod shock tube because it combines a
rarefaction, contact discontinuity, and shock in one known solution. Beginning
in 1D isolates the kinetic equilibrium, flux, collision, time integration, and
boundary treatment before geometry and y-direction fluxes are introduced.

This does not prove that a 2D wedge solver works. It proves that the core
compressible transport and thermodynamic recovery are credible enough to
justify the next dimension.

## 3. Separate velocity discretization from spatial discretization

The physical x-domain was divided into finite-volume cells. Each cell stores
five population values associated with the D1V5 molecular velocities. `V5`
does not require lattice streaming: the populations move through finite-volume
faces using numerical fluxes.

This distinction matters for later body-fitted or cut-cell geometries because
the spatial mesh and the discrete velocity model serve different purposes.

## 4. Early moment-matching attempt

An early D1V5 version constructed equilibrium populations by solving a
five-moment polynomial system. Although the code could run, it produced poor
wave speeds, excessive velocities, incorrect temperature behavior, and at
times nonphysical states. Increasing to D1V7 did not fix the underlying
closure problem.

The lesson was that adding velocities does not repair an equilibrium model
whose moment constraints and internal-energy treatment are not physically
consistent with the intended gas.

## 5. Adopt the published KT equilibrium

The implementation was changed to the Kataoka–Tsutahara D1V5 equilibrium with
the published coefficient structure. The rest population carries the internal
energy term `eta0^2`, allowing the energy moment and specific-heat ratio to be
represented consistently.

This was the major model correction. It separated problems caused by the
equilibrium closure from problems caused by spatial transport.

## 6. Build finite-volume population transport

For every population and cell, a flux enters through the west face and leaves
through the east face. The sign of the discrete molecular velocity decides
which reconstructed face state is upstream. The transport contribution is the
negative face-flux difference divided by cell width.

First-order upwind transport was useful as a stable baseline, but it smeared
the rarefaction, contact, and shock. Higher spatial resolution alone could not
remove this built-in numerical diffusion efficiently.

## 7. Add MUSCL reconstruction and minmod

MUSCL reconstructed a linear variation inside each cell instead of treating
the entire cell as constant. Standard minmod limited the slope to prevent new
oscillations near discontinuities.

The result was stable, but the smallest safe slope was often selected. This
flattened gradients and left visible diffusion error.

## 8. Improve time integration

Forward Euler was replaced by SSP-RK3. Three right-hand-side evaluations are
blended during every accepted step, improving temporal accuracy and making the
MUSCL update more robust.

The time step is limited by both transport and collision:

```text
dt_adv = CFL*dx/max(|c_i|)
dt_col = 0.25*tau
dt     = min(dt_adv, dt_col).
```

Physical-state checks are performed after every RK stage. A failed stage is
rejected instead of silently accepted.

## 9. Test the MC limiter

The monotonized-central limiter allowed larger slopes and reduced average
smearing. Its density, velocity, pressure, and temperature L1 errors were
smaller than standard minmod in the matched comparison.

However, the sharper population reconstruction created a temperature
overshoot. The overshoot became an explicit metric instead of being assessed
only by visual inspection.

## 10. Try hybrid limiters

Two hybrid approaches were explored:

1. **Population hybrid:** every population independently chose MC or minmod.
   This could treat the five populations inconsistently inside one physical
   cell, disturbing moments such as temperature.
2. **Macroscopic hybrid:** density, pressure, and temperature created one
   troubled-cell mask shared by all populations. This was more physically
   consistent, but the binary switch between full MC and full minmod remained
   abrupt and did not provide the desired accuracy/overshoot balance.

These versions are retained in `matlab/archive/` because unsuccessful methods
explain the final design.

## 11. Replace binary switching with generalized minmod

Generalized minmod uses

```text
sigma = minmod(theta*dL, (dL+dR)/2, theta*dR),  1 <= theta <= 2.
```

`theta=1` is the minmod end of the family and `theta=2` is the MC end. The
current `theta=1.20` allows a controlled increase in slope without jumping
directly to full MC.

At `Nx=6250`, this reduced all four global L1 errors by about 17–20% relative
to minmod and reduced MC's temperature peak excess by about 75%. It remains a
compromise: MC still had smaller average errors, and minmod still had smaller
peak excess.

## 12. Quantify error

The solver reports:

- global L1 error for an overall average difference;
- wave-region L1 error so undisturbed domain length does not hide wave error;
- global and wave-region L2 error for stronger weighting of larger local
  differences;
- temperature L-infinity error for the worst local temperature difference;
- signed temperature peak excess and percentage; and
- retry and fallback counts as stability diagnostics.

Using several metrics prevents one favorable number from hiding a different
failure mode.

## 13. Isolate boundaries with a longer domain

The working domain was increased to `Lx=20`, placing the initial diaphragm at
`x0=10`. This keeps the waves away from the transmissive boundaries during the
reported interval.

There is an important control issue: if `Nx` stays fixed while `Lx` increases,
then `dx` also increases and the simulation becomes spatially coarser. A clean
boundary study should increase `Nx` with `Lx` so `dx` remains fixed.

## 14. Perform the grid study

The generalized-minmod solver was run at 6,250, 12,500, 25,000, and 50,000
cells. Most L2 errors decreased with refinement, but the observed order was
only about 0.26–0.37 and finest-grid velocity error increased slightly.

The study therefore supports grid improvement, not formal second-order
convergence. Discontinuities, fixed relaxation time, time-step behavior, and
exact-wave/grid alignment all contribute to the measured trend.

## 15. Current decision point

The 1D solver is accurate enough to support planning a controlled move toward
2D, but one more 1D study should separate spatial and finite-relaxation effects
by scaling `tau` and `dt` with `dx`. The two-dimensional discrete-velocity
model should be selected from published work and its required moments derived
before implementation begins.

## 16. Validation required before a wedge claim

The next solver should pass, in order:

1. published two-dimensional equilibrium-moment checks;
2. planar shock-tube y-uniformity and D1 reference agreement;
3. a genuinely two-dimensional shock-tube benchmark;
4. conservation and grid studies;
5. stable wall and inflow/outflow boundary tests; and
6. oblique-shock angle comparison for a wedge.

Only after those checks should thermodynamic nonequilibrium moments be used to
draw physical conclusions about a wedge shock.
