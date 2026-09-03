# KT-D1V5 Generalized-Minmod L2 Convergence Test

## Test setup

- Model: KT-D1V5 finite-volume discrete Boltzmann method
- Transport: MUSCL reconstruction with generalized minmod
- Limiter parameter: `theta = 1.20`
- Time integration: SSP-RK3
- Relaxation time: `tau = 5.0e-5`
- CFL number: `0.10`
- Domain length: `Lx = 20`
- Final time: `tEnd = 0.15`
- Resolutions: `Nx = 6250, 12500, 25000, 50000`
- Benchmark: exact inviscid Euler solution for the Sod shock tube
- No time-step retries, limiter fallbacks, or positivity repairs occurred during any run.

Only the number of cells and grid spacing changed between runs.

## Global L2 results

| Nx | dx | L2(rho) | Order | L2(u) | Order | L2(p) | Order | L2(T) | Order |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 6,250 | 0.0032 | 1.960774e-3 | — | 8.419527e-3 | — | 1.883587e-3 | — | 4.674423e-3 | — |
| 12,500 | 0.0016 | 1.593604e-3 | 0.299 | 6.593089e-3 | 0.353 | 1.443274e-3 | 0.384 | 3.640293e-3 | 0.361 |
| 25,000 | 0.0008 | 1.257972e-3 | 0.341 | 3.583238e-3 | 0.880 | 1.114715e-3 | 0.373 | 2.696154e-3 | 0.433 |
| 50,000 | 0.0004 | 1.131536e-3 | 0.153 | 3.871277e-3 | -0.112 | 1.023180e-3 | 0.124 | 2.428375e-3 | 0.151 |

## Wave-region L2 results

The wave region was fixed at `x0 - 0.30 <= x <= x0 + 0.40` for every grid.

| Nx | L2(rho) | Order | L2(u) | Order | L2(p) | Order | L2(T) | Order |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 6,250 | 1.047479e-2 | — | 4.497857e-2 | — | 1.006245e-2 | — | 2.497158e-2 | — |
| 12,500 | 8.513308e-3 | 0.299 | 3.522142e-2 | 0.353 | 7.710222e-3 | 0.384 | 1.944708e-2 | 0.361 |
| 25,000 | 6.724143e-3 | 0.340 | 1.915321e-2 | 0.879 | 5.958404e-3 | 0.372 | 1.441155e-2 | 0.432 |
| 50,000 | 6.048317e-3 | 0.153 | 2.069285e-2 | -0.112 | 5.469127e-3 | 0.124 | 1.298021e-2 | 0.151 |

## Temperature peak excess

| Nx | Peak excess | Percentage of exact peak |
|---:|---:|---:|
| 6,250 | 2.048193e-3 | 0.179443% |
| 12,500 | 2.128815e-4 | 0.018651% |
| 25,000 | 8.115870e-6 | 0.000711% |
| 50,000 | 3.290399e-4 | 0.028827% |

Peak excess improved strongly through 25,000 cells but was not monotonic at 50,000 cells. This metric depends on a single maximum cell and is sensitive to exact discontinuity placement, so it should be considered together with L2 rather than used alone.

## Interpretation

From 6,250 to 50,000 cells, the global L2 errors decreased by approximately:

- Density: 42.3%
- Velocity: 54.0%
- Pressure: 45.7%
- Temperature: 48.0%

The overall effective orders across all three grid doublings were approximately `0.264` for density, `0.374` for velocity, `0.293` for pressure, and `0.315` for temperature.

The test demonstrates grid improvement, but it does not show clean second-order asymptotic convergence. This is not automatically a failure: the Sod solution contains discontinuities, where a MUSCL method cannot retain its second-order smooth-region accuracy. In addition, `tau` and the collision-limited time step remained fixed while the grid was refined. At high resolution, spatial error becomes smaller while finite-tau model error, time-integration error, and grid alignment with exact discontinuities become more visible. This likely explains the reduced orders and the slight increase in velocity L2 on the finest grid.

## Recommended follow-up

For a cleaner measurement of the spatial discretization order, compare each grid with a much finer DBM reference solution at the same `tau`, interpolate that reference onto each coarser grid, and reduce the time step with grid spacing. The present test is still the appropriate engineering comparison against the exact Euler Sod benchmark.
