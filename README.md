# DBFVM Research

This repository contains finite-volume discrete Boltzmann method (FVDBM)
research code for compressible-flow validation and shock-wave applications.
The current reference problem is the one-dimensional Sod shock tube; the next
stage is controlled two-dimensional validation before wedge-flow studies.

## Recommended starting point

Run the current 1D reference solver in MATLAB:

```matlab
cd matlab/d1v5
d1v5_sodshock_gminmod_rk3
```

It solves the classic Sod problem
`(rho,u,p)_L = (1,0,1)` and `(rho,u,p)_R = (0.125,0,0.1)` with:

- KT-D1V5 equilibrium and internal-energy placement;
- finite-volume population transport;
- MUSCL generalized-minmod reconstruction with `theta = 1.20`;
- velocity-sign upwinding;
- BGK collision;
- SSP-RK3 time integration; and
- transmissive boundaries on a long domain.

The default validation settings are `Nx = 50000`, `Lx = 20`, `tau = 5e-5`,
`CFL = 0.10`, and `tEnd = 0.15`.

## What is verified

The four-grid MATLAB study completed at `Nx = 6250, 12500, 25000, 50000`
without retries, fallback limiters, or positivity repairs. Global L2 error
decreased overall for density, velocity, pressure, and temperature. The
effective orders were only about `0.26–0.37`, so the present results do **not**
support a formal second-order convergence claim.

See [the convergence report](reports/d1v5_gminmod_L2_convergence_report.md) and
[the CSV results](results/d1v5_gminmod_L2_convergence_results.csv).
The cleaned reference copy also passed a MATLAB `Nx = 400` smoke test with
12,000 accepted steps, positive density and temperature, and no retries or
fallbacks.

## Directory guide

- `matlab/d1v5/`: current solver, limiter-history comparison, and L2 driver.
- `matlab/d2v9/`: 2D planar-Sod teaching prototype; implemented but not fully
  validated.
- `matlab/archive/`: earlier MC and hybrid-limiter development snapshots.
- `python/`: earlier JAX D2V9 finite-volume prototype and exact-Riemann test.
- `reports/`: concise descriptions of the limiter development and validation.
- `results/`: compact tables that are useful to version-control.
- `CHANGELOG.md`: concise record of what entered each curated snapshot.

The archived files document the development path; they are not recommended as
the starting point for new simulations.

## Next controlled study

Hold the generalized-minmod limiter fixed and vary `tau` with grid spacing.
Record whole-domain and wave-region L2/Linf errors, temperature peak excess,
and wave-location errors. This separates spatial error from finite-relaxation
and collision-time-step effects before moving to 2D accuracy work.

Model, numerical-method, and benchmark citations are listed in
[SOURCES.md](SOURCES.md). The repository intentionally distinguishes verified
results from working prototypes; in particular, the D2V9 code is not yet a
validated two-dimensional solver.
