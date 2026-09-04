# MATLAB Solver Guide

The MATLAB code is intentionally comment-heavy. It is both a research solver
and a record of how each numerical operation works.

## Recommended order

| Order | File | Purpose |
|---:|---|---|
| 1 | [`d1v5/d1v5_sodshock_gminmod_rk3.m`](d1v5/d1v5_sodshock_gminmod_rk3.m) | Run one KT-D1V5 Sod case and compare it with the exact solution |
| 2 | [`d1v5/d1v5_sodshock_limiter_history_comparison.m`](d1v5/d1v5_sodshock_limiter_history_comparison.m) | Compare every limiter developed during the project under matched settings |
| 3 | [`d1v5/run_d1v5_gminmod_L2_convergence.m`](d1v5/run_d1v5_gminmod_L2_convergence.m) | Run the four-grid L2 study |
| 4 | [`d2v9/d2v9_kt_muscl_sod_teaching.m`](d2v9/d2v9_kt_muscl_sod_teaching.m) | Inspect the current two-dimensional extension |

## Main D1V5 solver

Run from this repository with:

```matlab
cd matlab/d1v5
d1v5_sodshock_gminmod_rk3
```

The script performs the following sequence:

1. Defines the gas, grid, time, and limiter parameters.
2. Creates the left and right Sod states.
3. Builds the KT-D1V5 equilibrium populations.
4. Recovers macroscopic moments from the populations.
5. Reconstructs population values at cell faces with MUSCL.
6. Applies velocity-sign upwinding at every face.
7. Adds BGK collision to the negative flux divergence.
8. Advances the populations with three SSP-RK3 stages.
9. Rejects a time step and reduces `dt` if a stage becomes nonphysical.
10. Compares the final fields with the exact Euler solution.
11. Prints error and temperature-peak diagnostics.
12. Saves a `.mat` file for later post-processing.

### Default parameters

| Parameter | Default | Meaning |
|---|---:|---|
| `gamma` | `1.4` | Ratio of specific heats |
| `b` | `5` | Internal degrees-of-freedom parameter computed from `gamma` |
| `v1` | `1.0` | Slow KT molecular speed |
| `v2` | `3.0` | Fast KT molecular speed |
| `eta0` | `1.75` | Internal-energy value carried by the rest population |
| `Nx` | `50000` | Number of finite-volume cells |
| `Lx` | `20.0` | Physical domain length |
| `tau` | `5e-5` | BGK relaxation time |
| `CFL` | `0.10` | Advective time-step safety factor |
| `tEnd` | `0.15` | Final physical time |
| `limiterTheta` | `1.20` | Generalized-minmod sharpness parameter |

### Controlled overrides

The reference solver accepts environment-variable overrides. MATLAB's
`setenv` is useful because the script begins with `clear`.

```matlab
setenv('DBM_NX','6250')
setenv('DBM_THETA','1.20')
setenv('DBM_TAU','5e-5')
d1v5_sodshock_gminmod_rk3

% Clear the overrides when the run is finished.
setenv('DBM_NX','')
setenv('DBM_THETA','')
setenv('DBM_TAU','')
```

Only change one controlled variable at a time when measuring its effect.

## How the finite-volume transport works

Each physical cell stores five population values. For population `q`, MUSCL
uses the neighboring cell values to estimate a limited slope. That slope
creates a left and right reconstructed value at each face.

- If `c(q) > 0`, the face receives the state reconstructed from the cell on
  its left because that population moves right.
- If `c(q) < 0`, the face receives the state reconstructed from the cell on
  its right because that population moves left.
- The difference between the east and west face fluxes, divided by `dx`, is
  the finite-volume flux divergence for that cell.

The transport and collision pieces are then combined:

```text
RHS = -flux divergence + (equilibrium - population)/tau.
```

## Generalized minmod

For left and right differences `dL` and `dR`, the implemented slope is

```text
sigma = minmod(theta*dL, (dL+dR)/2, theta*dR).
```

- `theta = 1` behaves like standard minmod: safer but more diffusive.
- `theta = 2` reaches standard MC: sharper but more prone to overshoot here.
- `theta = 1.20` is the current measured compromise.
- If the three candidates do not share a sign, the slope is zero. This avoids
  constructing a new local maximum or minimum.

## SSP-RK3 update

The right-hand side is evaluated three times:

```text
f1    = f + dt*RHS(f)
f2    = 3/4*f + 1/4*(f1 + dt*RHS(f1))
f_new = 1/3*f + 2/3*(f2 + dt*RHS(f2))
```

Boundary conditions and physical-state checks are applied at every stage.
The solver reports how often a failed stage required a smaller time step or a
more diffusive fallback limiter. The reproduced `Nx=6250` run required none.

## Limiter-history comparison

```matlab
setenv('DBM_NX','6250')  % optional faster documented comparison
d1v5_sodshock_limiter_history_comparison
setenv('DBM_NX','')
```

This runs minmod, MC, population-by-population hybrid limiting, shared
macroscopic hybrid limiting, and generalized minmod with identical settings.
It prints a table of global L1 error, wave-region L1 error, and temperature
peak excess, then saves the complete workspace data.

## L2 grid study

```matlab
convergenceTable = run_d1v5_gminmod_L2_convergence;
```

This is an expensive calculation because it runs the reference solver at
`Nx = 6250, 12500, 25000, 50000`. It writes CSV, MAT, and PNG outputs to the
top-level `results/` directory. The committed CSV and figure record the
completed study, so rerunning is only necessary after changing the solver.

## D2V9 status

The D2V9 program expands storage, moments, and fluxes into x and y. It uses a
KT-style nine-velocity equilibrium, MUSCL-minmod reconstruction in both
directions, BGK collision, and a planar Sod initial condition.

It is not yet equivalent to a validated 2D research solver. Before reporting
its results, check:

1. every required equilibrium moment against the published KT equations;
2. conservation of mass, x/y momentum, and energy;
3. uniformity across y for a planar x-directed Sod problem;
4. agreement of the centerline with the validated D1V5/Euler result;
5. grid convergence in both x and y; and
6. sensitivity to `tau`, CFL, and boundary placement.

## Archive

The [`archive/`](archive/) directory preserves limiter versions that were
important to the reasoning but were superseded. They are not recommended as
starting points for new simulations.
