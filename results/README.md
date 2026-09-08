# Versioned Numerical Results

This folder holds the figures and numerical values I use in the main README.
I do not commit the generated `.mat` workspaces because they are large and I
can recreate them from the MATLAB scripts.

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

## What the error measurements mean

- **Global L1:** the average absolute difference between my result and the
  exact result over every cell. Long, unchanged parts of the tube can make
  this number look smaller than the error around the waves.
- **Wave-region L1:** the same calculation, but only on
  `x0-0.30 <= x <= x0+0.40`, where the waves are located.
- **Global L2:** the root-mean-square difference over the whole domain. It
  gives more weight to large local errors than L1 does.
- **Wave-region L2:** root-mean-square difference restricted to the wave zone.
- **Temperature peak excess:** maximum numerical temperature minus the exact
  maximum temperature in the wave zone. Positive means overshoot.
- **Observed order:** change in error between two grid spacings, computed as
  `log(E_coarse/E_fine)/log(dx_coarse/dx_fine)`.

### How I chose the wave region

The easiest way I think about this is that global L1 grades the entire graph,
while wave-region L1 zooms in on the part where something is actually
happening. Most of the long shock tube is still at its original constant
state. If I average over all of those easy cells, they can make the total error
look better than the fit around the rarefaction, contact, and shock really is.

For the current case, the diaphragm starts at `x0 = 10` and the solution is
measured at `tEnd = 0.15`. The exact Sod solution tells me how far the three
waves have traveled by that time. I use the slightly wider interval

```text
x0 - 0.30 <= x <= x0 + 0.40
```

or `9.70 <= x <= 10.40`. That interval contains the rarefaction, contact, and
shock, plus a little extra space for numerical smearing around their edges.
The wave-region L1 error is just the normal average absolute error calculated
inside that window.

The solver is **not automatically detecting the wave region**. This fixed
window was selected for the current initial conditions and final time. If I
change `x0`, `tEnd`, or the left and right gas states, I need to calculate the
new exact wave locations and update the window.

## Reproducibility note

I regenerated the `Nx=6250` reference and limiter figures on September 3,
2026 from the committed MATLAB code. The four-grid files contain the finished
`Nx=6250–50000` study. Runtime depends on the computer. I only compare the
error values when the physical and numerical parameters match.
