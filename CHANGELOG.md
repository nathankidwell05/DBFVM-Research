# Project snapshot history

## 2026-09-18

- Added `reports/d1v5_mathematics/`: a LaTeX report (with PDF) deriving the
  KT-D1V5 equilibrium, a Chapman–Enskog analysis of the discrete model
  (including the contact diffusivity `D_c = tau (v1^2-U^2)(v2^2-U^2)/((b+2)T)`),
  the exact Sod solution, MUSCL/generalized-minmod TVD theory, SSP-RK3
  stability, and the error norms, plus the scripts that verify each result and
  build every figure and table.
- Changed the solver defaults to `tau = 5e-6`, `CFL = 0.025`, and added the
  `DBM_CFL` and `DBM_DISABLE_LIVE_PLOT` overrides.
- Added the convergence driver (`Nx = 3125, 6250, 12500, 25000, 50000`) with
  per-grid caching, L1 errors, effective-CFL records, settings-stamped output
  files, and `onCleanup` removal of the `DBM_*` overrides.
- Set the limiter-history comparison to `tau = 5e-6`, `CFL = 0.025`,
  `Nx = 6250`, added L2 errors, and gave its outputs settings-stamped names.

## 2026-09-08

- Added a limiter and stability guide covering first-order upwind, minmod, MC,
  both hybrid attempts, generalized minmod, physical-state detection, retry
  switching, and the shortcomings of every method.
- Explained how the fixed Sod wave window was selected and why wave-region
  errors are reported separately from whole-domain errors.
- Expanded the limiter guide with the complete code path, numerical examples,
  MATLAB decision logic, and plain-language explanations of each method's
  weaknesses.

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
