# Project snapshot history

## 2026-09-02

- Added the current KT-D1V5 generalized-minmod SSP-RK3 solver with clean file
  naming and validated defaults.
- Added the matched five-limiter history comparison.
- Added the four-grid L2 convergence driver, report, and CSV results.
- Preserved the full-MC and hybrid-limiter experiments under `matlab/archive/`.
- Added the D2V9 planar-Sod teaching prototype with an explicit unvalidated
  status.
- Added the earlier JAX D2V9 prototype and its exact-Riemann unit test.
- Removed dead positivity-repair and scalar-limiter code from the reference
  solver; nonphysical states are rejected rather than silently repaired.
- Standardized the reference plot as exact black solid and generalized minmod
  green dashed.
- Added concise README files and ignored generated MATLAB/Python artifacts.
