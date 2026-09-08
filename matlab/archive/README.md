# Limiter-development archive

These are older versions of my one-dimensional solver. Some were unsuccessful
and some were simply replaced by a better approach. I kept them because they
show what I tried and why I changed the limiter.

- `d1v5_sodshock_muscl_MC_RK3_before_error_metrics_patch.m`: full-MC stage
  before quantitative error and peak diagnostics.
- `d1v5_sodshock_before_shared_hybrid_sensor.m`: population-by-population
  hybrid limiting, which could reconstruct the five populations inconsistently.
- `d1v5_sodshock_before_generalized_minmod.m`: shared macroscopic troubled-cell
  sensor, still using an abrupt MC/minmod switch.
- `test_shared_hybrid_sensor.m`: focused experiment for the shared sensor.

The version I currently use is generalized minmod with `theta = 1.20`. It is
in `../d1v5/`.
