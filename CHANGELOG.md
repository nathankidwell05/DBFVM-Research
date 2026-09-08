# Project snapshot history

## 2026-09-08

- Added a limiter and stability guide covering first-order upwind, minmod, MC,
  both hybrid attempts, generalized minmod, physical-state detection, retry
  switching, and the shortcomings of every method.
- Explained how the fixed Sod wave window was selected and why wave-region
  errors are reported separately from whole-domain errors.

## 2026-09-05

- Removed unrelated prototype code so the repository contains only the MATLAB
  research workflow.
- Added reproducible MATLAB post-processing for spatial error, limiter error
  bars, and the accuracy-versus-overshoot tradeoff.
- Added three more figures to the main README and results guide.

## 2026-09-03

- Expanded the main and MATLAB READMEs into a complete project guide.
- Added a chronological methodology and decision record.
- Added a controlled `DBM_NX` override to the limiter-history script.
- Regenerated the MATLAB D1V5 reference and five-limiter comparison at
  `Nx=6250`.
- Added three result figures and two machine-readable metric tables.

## 2026-09-02

- Added the current KT-D1V5 generalized-minmod SSP-RK3 solver with clean file
  naming and validated defaults.
- Added the matched five-limiter history comparison.
- Preserved the full-MC and hybrid-limiter experiments under `matlab/archive/`.
- Removed dead positivity-repair and scalar-limiter code from the reference
  solver; nonphysical states are rejected rather than silently repaired.
- Standardized the reference plot as exact black solid and generalized minmod
  green dashed.
- Added concise README files and ignored generated MATLAB artifacts.
