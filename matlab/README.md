# MATLAB solvers

## D1V5 reference workflow

From `project/matlab/d1v5`:

```matlab
d1v5_sodshock_gminmod_rk3
convergenceTable = run_d1v5_gminmod_L2_convergence;
d1v5_sodshock_limiter_history_comparison
```

- `d1v5_sodshock_gminmod_rk3.m` is the current single-case reference solver.
- `run_d1v5_gminmod_L2_convergence.m` reruns that solver on four grids and
  writes convergence data under `results/`.
- `d1v5_sodshock_limiter_history_comparison.m` compares minmod, MC, two failed
  hybrid ideas, and generalized minmod under matched settings.

The files are intentionally comment-heavy because they also serve as teaching
material. The exact solution is the inviscid Euler benchmark; the DBM solver
has finite relaxation, so comparison errors contain both numerical and model
effects.

## D2V9 status

`d2v9/d2v9_kt_muscl_sod_teaching.m` extends storage and fluxes into two physical
dimensions for a planar Sod case. It should be treated as an unvalidated
prototype until its equilibrium moments, y-uniformity, conservation, and grid
convergence have been checked in MATLAB.

## Archive

`archive/` preserves the sequence that exposed the temperature overshoot:
full MC, population-specific hybrid limiting, and a shared macroscopic sensor.
Those files are evidence of the development process, not production solvers.
