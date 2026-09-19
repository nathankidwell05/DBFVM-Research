# KT-D1V5 solver guide

[`d1v5_kt_solver_guide.pdf`](d1v5_kt_solver_guide.pdf) is a 24-page report on the
KT-D1V5 finite-volume discrete Boltzmann solver: the physical model and numerical
method, grid convergence and order of accuracy on the Sod shock tube, a limiter study,
and an accuracy-versus-cost comparison with first-order and MUSCL versions of the DBM
and with a conventional finite-volume Euler solver (HLLC flux). Every figure has a
discussion of what it shows and what it means.

All figures are generated in MATLAB by `scripts/make_guide_figures.m` and are stored in
[`../figures/`](../figures/) as PDF and PNG (MATLAB `.fig` copies are produced locally
but are not tracked by Git).

## Build

```bash
tectonic d1v5_kt_solver_guide.tex
```

## Reproduce

From `scripts/` in MATLAB:

```matlab
run_d1v5_limiter_grid('firstorder',6250)   % one DBM comparison run (repeat for each limiter and grid)
run_cost_benchmark                          % Euler solver runs and clean DBM timing (run alone)
make_guide_figures                          % all figures and data/tab_methods.tex
```

The generalized-minmod DBM runs come from `matlab/d1v5/run_d1v5_gminmod_L2_convergence.m`.
`euler_sod_solver.m` is the conventional Euler solver used for comparison.
