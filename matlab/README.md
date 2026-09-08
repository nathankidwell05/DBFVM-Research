# MATLAB Solver Guide

I intentionally left a lot of comments in the MATLAB code. I wanted the files
to work as solvers, but I also wanted to be able to come back later and
understand what every major step is doing.

## Recommended order

| Order | File | Purpose |
|---:|---|---|
| 1 | [`d1v5/d1v5_sodshock_gminmod_rk3.m`](d1v5/d1v5_sodshock_gminmod_rk3.m) | Run one KT-D1V5 Sod case and compare it with the exact solution |
| 2 | [`d1v5/d1v5_sodshock_limiter_history_comparison.m`](d1v5/d1v5_sodshock_limiter_history_comparison.m) | Compare every limiter developed during the project under matched settings |
| 3 | [`d1v5/generate_documentation_figures.m`](d1v5/generate_documentation_figures.m) | Generate extra error and limiter-comparison figures from saved MATLAB runs |

## Main D1V5 solver

Run from this repository with:

```matlab
cd matlab/d1v5
d1v5_sodshock_gminmod_rk3
```

The main script works through the problem in this order:

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

The reference solver lets me override a few parameters with environment
variables. I use `setenv` because the script starts with `clear`.

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

When I am measuring what a parameter does, I only change one at a time.

## How the finite-volume transport works

Each physical cell stores five population values. For one population `q`,
MUSCL looks at neighboring cells and estimates a limited slope. I use that
slope to reconstruct a value on the left and right side of the cell.

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

I wrote out how every limiter works, how the hybrid versions decide to switch,
and where each method falls short in the
[limiter and stability guide](../LIMITERS.md).

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

I apply the boundary conditions and check the physical state after every
stage. If a trial stage gives a bad state, the solver reduces the time step
and can fall back to a safer limiter. The reproduced `Nx=6250` run did not need
any retries or fallbacks.

## Limiter-history comparison

```matlab
setenv('DBM_NX','6250')  % optional faster documented comparison
d1v5_sodshock_limiter_history_comparison
setenv('DBM_NX','')
```

This runs minmod, MC, both hybrid attempts, and generalized minmod with the
same setup. It prints global L1 error, wave-region L1 error, and temperature
peak excess so I can compare more than just the appearance of the curves.

## Documentation figures

After running both the main solver and limiter-history comparison, generate
the additional spatial-error and tradeoff figures with:

```matlab
generate_documentation_figures
```

This script reads the two saved `.mat` files and exports three PNG figures to
the top-level `results/` folder. It only makes figures; it does not rerun the
solver or change the results.

## Archive

The [`archive/`](archive/) folder keeps the limiter versions that helped me
figure out what was going wrong. They are there to show the development path,
not because they are the best files to use for a new run.
