# Limiter-development archive

These MATLAB snapshots preserve unsuccessful or superseded stages of the 1D
solver. They are kept so the numerical decisions remain reproducible.

- `d1v5_sodshock_muscl_MC_RK3_before_error_metrics_patch.m`: full-MC stage
  before quantitative error and peak diagnostics.
- `d1v5_sodshock_before_shared_hybrid_sensor.m`: population-by-population
  hybrid limiting, which could reconstruct the five populations inconsistently.
- `d1v5_sodshock_before_generalized_minmod.m`: shared macroscopic troubled-cell
  sensor, still using an abrupt MC/minmod switch.
- `test_shared_hybrid_sensor.m`: focused experiment for the shared sensor.

The current method is generalized minmod with `theta = 1.20`, located in
`../d1v5/`.
