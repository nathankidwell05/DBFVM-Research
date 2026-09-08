# Finite-Volume Discrete Boltzmann Research

This repository is a record of my work on a finite-volume discrete Boltzmann
method (FVDBM) for compressible flow. I am starting with a problem I can check
carefully in one dimension, then building toward two-dimensional shock tubes
and supersonic flow over a wedge.

## Research path

```text
1D Sod shock tube
        ↓
2D planar and genuinely 2D shock-tube validation
        ↓
2D supersonic wedge flow
        ↓
3D validation and wedge flow
        ↓
Possible airfoil applications
```

Right now I am finishing the one-dimensional stage. The D1V5 solver runs and
can be compared directly with the exact Sod solution. I have not selected or
implemented the two-dimensional velocity model yet.

| Part of the project | Current status |
|---|---|
| KT-D1V5 Sod solver | Working MATLAB reference implementation |
| Exact Euler comparison | Implemented for density, velocity, pressure, and temperature |
| Limiter comparison | Minmod, MC, two hybrid attempts, and generalized minmod compared |
| Error measurement | Global and wave-region L1 errors plus temperature peak excess |
| Wedge solver | Not yet implemented |

## Why begin with the Sod shock tube?

The Sod problem begins with two stationary gas states separated by a
diaphragm:

| State | Density, `rho` | Velocity, `u` | Pressure, `p` |
|---|---:|---:|---:|
| Left | 1.000 | 0.000 | 1.000 |
| Right | 0.125 | 0.000 | 0.100 |

Removing the diaphragm produces three different wave features:

- a smooth rarefaction wave moving left;
- a contact discontinuity moving right; and
- a shock moving farther to the right.

This gives me several useful tests in one problem: smooth transport, contact
and shock capture, wave speed, pressure and temperature recovery, stability,
and boundary treatment. There is also an exact inviscid Euler solution. That
means I can calculate the error instead of only deciding whether the graph
looks right.

## Current numerical model

The code advances five discrete particle populations with

```text
df_i/dt + c_i df_i/dx = (f_i^eq - f_i)/tau.
```

| Component | Choice in the reference solver | Purpose |
|---|---|---|
| Kinetic model | Kataoka–Tsutahara D1V5 | Represents compressible mass, momentum, and energy moments |
| Discrete velocities | `[-3, -1, 0, 1, 3]` | Moves five populations in one physical dimension |
| Internal-energy parameter | `eta0 = 1.75` on the rest population | Allows the desired specific-heat ratio |
| Collision | Single-relaxation-time BGK | Relaxes populations toward local equilibrium |
| Spatial discretization | Cell-centered finite volume | Computes population transport through cell faces |
| Face reconstruction | MUSCL | Reconstructs left and right states at each face |
| Limiter | Generalized minmod, `theta = 1.20` | Balances minmod stability against MC sharpness |
| Upwinding | Based on the sign of each `c_i` | Selects the population arriving at a face |
| Time integration | SSP-RK3 | Reduces temporal error and preserves strong-stability properties |
| Boundaries | Transmissive copying | Lets waves leave with minimal reflection when they reach a boundary |

### D1V5 does not mean lattice streaming

`D1V5` only describes the **velocity model**: one physical dimension and five
discrete molecular velocities. It does not automatically make the code a
lattice Boltzmann solver. I still divide the domain into finite-volume cells
and calculate fluxes at their faces. The populations do not have to move
exactly from one lattice node to the next in one time step.

## Run the MATLAB code

The main implementation is
[`d1v5_sodshock_gminmod_rk3.m`](matlab/d1v5/d1v5_sodshock_gminmod_rk3.m).

From MATLAB:

```matlab
cd matlab/d1v5
d1v5_sodshock_gminmod_rk3
```

The saved setup uses `Nx = 50000`, `Lx = 20`, `tau = 5e-5`, `CFL = 0.10`,
`tEnd = 0.15`, and `theta = 1.20`. To do a faster test without editing the
file, I can temporarily override the grid:

```matlab
setenv('DBM_NX','6250')
d1v5_sodshock_gminmod_rk3
setenv('DBM_NX','')
```

Additional workflows and output definitions are explained in the
[MATLAB guide](matlab/README.md).

## Reproduced one-dimensional result

I regenerated the result below from the saved MATLAB code using `Nx = 6250`.
It completed 12,000 accepted steps in 39.96 seconds on my computer, with no
retries and no fallback limiters. The runtime will change from one computer to
another, so the error values are the more useful comparison.

| Global L1 error | Value |
|---|---:|
| Density | `1.534213e-4` |
| Velocity | `3.150135e-4` |
| Pressure | `1.384924e-4` |
| Temperature | `2.222085e-4` |

The temperature peak exceeded the exact peak by `2.048193e-3`, or `0.179443%`.

![Generalized-minmod FVDBM compared with the exact Sod solution](results/d1v5_gminmod_reference_Nx6250.png)

Full machine-readable values are in
[`d1v5_gminmod_reference_Nx6250_metrics.csv`](results/d1v5_gminmod_reference_Nx6250_metrics.csv).

## What the limiter study showed

All five cases below used the same `Nx = 6250`, physical model, initial
condition, domain, relaxation time, CFL number, and SSP-RK3 integrator. Only
the reconstruction limiter changed.

| Limiter | Global L1(T) | Wave L1(T) | Temperature peak excess |
|---|---:|---:|---:|
| Minmod | `2.791803e-4` | `7.967473e-3` | `7.964756e-4` |
| MC | `1.851267e-4` | `5.283297e-3` | `8.087784e-3` |
| Population hybrid | `2.355076e-4` | `6.721108e-3` | `4.650723e-3` |
| Macroscopic hybrid | `2.787301e-4` | `7.954627e-3` | `1.005029e-3` |
| Generalized minmod, `theta=1.20` | `2.222085e-4` | `6.341566e-3` | `2.048193e-3` |

There was no limiter that won every comparison:

- MC had the smallest average errors, but its peak excess was almost four
  times the generalized-minmod value.
- Minmod produced the smallest peak excess, but it had the largest global and
  wave-region errors because it smeared the solution more strongly.
- Generalized minmod reduced MC's peak excess by about 75% while reducing all
  four global L1 errors by about 17–20% relative to minmod.

That is why I currently use `theta = 1.20`. It is a compromise for this test,
not proof that it is the best setting for every flow problem.

![Limiter development comparison](results/d1v5_limiter_history_Nx6250.png)

![Temperature overshoot comparison](results/d1v5_temperature_limiter_comparison_Nx6250.png)

The next figure separates the error over the whole domain from the error only
around the waves. The wave-region values are larger because the long, flat
parts of the tube are no longer making the average look smaller.

![Limiter L1 error comparison](results/d1v5_limiter_error_bars_Nx6250.png)

The tradeoff plot shows the decision more clearly. Moving left means less
temperature error around the waves, while moving down means less overshoot.
MC is farthest left, but it also moves much higher because of its overshoot.
Generalized minmod stays between the minmod and MC results.

![Limiter accuracy and overshoot tradeoff](results/d1v5_limiter_tradeoff_Nx6250.png)

The spatial error plots show where my result differs from the exact solution.
Most of the error is around the wave edges, not in the constant parts of the
tube.

![Spatial absolute-error profiles](results/d1v5_gminmod_absolute_error_Nx6250.png)

## Repository map

| Path | Contents |
|---|---|
| [`matlab/d1v5/`](matlab/d1v5/) | Current 1D solver, limiter comparison, and figure-generation code |
| [`matlab/archive/`](matlab/archive/) | Earlier limiter implementations retained to document decisions |
| [`results/`](results/) | Versioned figures and machine-readable numerical results |
| [`reports/`](reports/) | Detailed limiter-development explanations |
| [`LIMITERS.md`](LIMITERS.md) | How each limiter detects risky gradients, switches behavior, and can fail |
| [`SOURCES.md`](SOURCES.md) | Papers and books used for model and numerical-method decisions |
| [`METHODOLOGY.md`](METHODOLOGY.md) | Chronological explanation of what was attempted and why |

## What this work does not prove yet

- The result I have validated so far is one-dimensional.
- The exact curve is an inviscid Euler solution, while the DBM calculation has
  finite relaxation time. The reported difference is therefore not purely
  spatial discretization error.
- The long domain keeps the waves away from the boundaries. Changing `Lx`
  without changing `Nx` also changes `dx`, so domain and resolution studies
  must be separated carefully.
- I have not selected and validated the two-dimensional discrete-velocity
  model yet.
- I do not have a validated wedge or airfoil result yet.

## What I plan to do next

1. Test how the choices of `tau` and time step affect the 1D result.
2. Measure individual rarefaction, contact, and shock position/width errors.
3. Select a published two-dimensional discrete-velocity model and verify all
   required equilibrium moments before implementing it.
4. Perform a genuine two-dimensional shock-tube test.
5. Add wall and inflow/outflow boundary conditions for a 2D wedge.
6. Compare the computed oblique-shock angle with oblique-shock theory before
   studying thermodynamic nonequilibrium quantities.

## References

I collected the model, limiter, finite-volume, and benchmark sources in
[`SOURCES.md`](SOURCES.md). The comments help explain the code, but the cited
papers are the source for the actual equations.
