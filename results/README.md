# Versioned Numerical Results

This directory contains the compact numerical evidence referenced by the main
README. Generated `.mat` workspaces remain ignored because they are large and
can be recreated from the MATLAB scripts.

## Files

| File | What it contains |
|---|---|
| `d1v5_gminmod_reference_Nx6250.png` | Four-field wave-region comparison between generalized minmod and exact Sod |
| `d1v5_gminmod_reference_Nx6250_metrics.csv` | Parameters, runtime, L1 errors, temperature diagnostics, retries, and fallbacks from that run |
| `d1v5_limiter_history_Nx6250.png` | Density, velocity, pressure, and temperature for all five limiter versions |
| `d1v5_temperature_limiter_comparison_Nx6250.png` | Enlarged view of limiter behavior near the temperature plateau |
| `d1v5_limiter_history_Nx6250_metrics.csv` | Matched global/wave L1 errors and peak excess for each limiter |
| `d1v5_limiter_error_bars_Nx6250.png` | Global and wave-region L1 errors for all fields and limiters |
| `d1v5_limiter_tradeoff_Nx6250.png` | Temperature accuracy-versus-overshoot tradeoff and normalized global errors |
| `d1v5_gminmod_absolute_error_Nx6250.png` | Location and magnitude of generalized-minmod absolute error |
| `d1v5_gminmod_L2_convergence_plot.png` | Whole-domain and wave-region L2 errors over four grids |
| `d1v5_gminmod_L2_convergence_results.csv` | Full convergence table, observed orders, runtime, and peak excess |

## Meaning of the measurements

- **Global L1:** mean absolute numerical-versus-exact difference over every
  cell. Long undisturbed regions can make this number look very small.
- **Wave-region L1:** the same mean absolute difference measured only on
  `x0-0.30 <= x <= x0+0.40`, where the waves exist.
- **Global L2:** root-mean-square difference over the entire domain. Larger
  local errors receive more weight than in L1.
- **Wave-region L2:** root-mean-square difference restricted to the wave zone.
- **Temperature peak excess:** maximum numerical temperature minus the exact
  maximum temperature in the wave zone. Positive means overshoot.
- **Observed order:** change in error between two grid spacings, computed as
  `log(E_coarse/E_fine)/log(dx_coarse/dx_fine)`.

## Reproducibility note

The `Nx=6250` reference and limiter figures were regenerated on September 3,
2026 from the committed MATLAB code. The four-grid convergence files record
the previously completed `Nx=6250–50000` study. Runtime values depend on the
computer; error values should be compared only when the physical and numerical
parameters are matched.
