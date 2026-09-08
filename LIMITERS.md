# Limiter and Stability Guide

This is my short guide to every limiter I tested in the KT-D1V5 solver. The
main questions are simple: what does the limiter look at, when does it become
more cautious, and what does it do badly?

## What a limiter actually does

MUSCL needs a slope inside every finite-volume cell so it can estimate the
population value on both faces. For one population, I first calculate

```text
dL = f(i)   - f(i-1)   change on the left side of cell i
dR = f(i+1) - f(i)     change on the right side of cell i
```

If `dL` and `dR` have the same sign, the population is still moving in the
same direction across the cell, so I may be able to keep a nonzero slope. If
the signs disagree, I am near a peak, dip, or sharp change. The limiters set
the slope to zero there so the reconstruction does not invent a new peak or
dip.

The important distinction is that a limiter only controls the local shape of
the reconstruction. It does **not** prove that the entire state is stable or
physical. I use separate checks after every SSP-RK3 stage to make sure density
and temperature are still valid.

## Methods tested

| Method | How it reacts | When it becomes conservative | Main advantage | Main shortcoming |
|---|---|---|---|---|
| First-order upwind | Sets every MUSCL slope to zero | Always | Most robust transport option | Strong numerical diffusion smears shocks, contacts, and smooth gradients |
| Minmod | Uses the smaller one-sided slope when `dL` and `dR` agree | Sets the slope to zero when their signs disagree | Stable and produces very little temperature overshoot in this test | Most diffusive MUSCL option tested; larger global and wave-region errors |
| MC | Applies `minmod(2*dL, (dL+dR)/2, 2*dR)` | Sets the slope to zero unless all three candidates have the same sign | Sharpest tested reconstruction and smallest average errors | Produced the largest temperature overshoot near a discontinuity |
| Population hybrid | Calculates MC and minmod, then chooses separately for each population | Uses minmod if curvature ratio is above `0.35` or the one-sided slope-size ratio is above `3.0` | Attempts to retain MC in smooth areas and add safety near abrupt population changes | Different populations in one physical cell can make different choices, disturbing the combined energy moment and recovered temperature |
| Macroscopic hybrid | Builds one warning mask from density, pressure, and temperature and shares it across all five populations | Uses minmod in warned cells and their immediate neighbors; otherwise uses MC | Keeps the five population reconstructions physically consistent | Binary MC/minmod switching is abrupt, sensor thresholds are problem-dependent, and the result remained close to diffusive minmod |
| Generalized minmod | Applies `minmod(theta*dL, (dL+dR)/2, theta*dR)` | Sign disagreement gives zero slope; `theta` continuously controls permitted slope size | Provides a simple continuous compromise between minmod and MC | A fixed `theta` cannot be optimal everywhere; increasing it reduces diffusion but can restore overshoot |

## How each one makes its decision

### First-order upwind

First-order upwind does not try to detect anything. It sets every slope to
zero, so the face just gets the value from the upstream cell. I only keep it
as an emergency fallback because it is safe but very smeary.

### Minmod

I think of minmod as asking two questions:

1. Do the left and right differences point in the same direction?
2. If so, which has the smaller magnitude?

If both answers are acceptable, minmod keeps the smaller slope. If the signs
disagree, it returns zero. That makes it hard to destabilize, but it also
flattens real gradients and spreads the shock or contact over more cells.

### Monotonized central (MC)

MC compares three candidates:

```text
2*dL,   (dL+dR)/2,   2*dR
```

It keeps the smallest candidate only when all three have the same sign. The
centered choice lets MC keep a steeper slope than minmod. That made the Sod
result sharper and lowered the average errors, but it also caused the largest
temperature overshoot.

### Population-by-population hybrid

My first hybrid attempt judged each population separately:

```text
curvature ratio = abs(dR-dL) / (abs(dR)+abs(dL)+small)
slope-size ratio = max(abs(dL),abs(dR)) / (min(abs(dL),abs(dR))+small)
```

It selected minmod when the curvature ratio exceeded `0.35` or the slope-size
ratio exceeded `3.0`; otherwise it selected MC. The problem is that one
population can switch to minmod while another population in the same physical
cell stays on MC. Temperature is recovered from all five populations together,
so these mixed decisions can throw off the energy moment and temperature.

### Shared macroscopic hybrid

For the second hybrid attempt, I made one decision from density, pressure, and
temperature and shared it across all five populations. For each field I
calculate a relative jump and a roughness value:

```text
relative jump = max(abs(dL),abs(dR)) / (abs(field(i))+small)
roughness     = abs(dR-dL) / (abs(dR)+abs(dL)+small)
```

A cell is marked as troubled when

```text
(roughness > 0.20 and relative jump > 0.0025)
or relative jump > 0.08.
```

I also mark one neighboring cell on each side. All five populations then use
minmod in the marked cells and MC everywhere else. This fixed the inconsistent
population choices, but the cutoff values are hand-chosen and the switch from
MC to minmod is still abrupt.

### Generalized minmod

The method I currently use is

```text
slope = minmod(theta*dL, (dL+dR)/2, theta*dR).
```

`theta = 1` is the minmod end of this family and `theta = 2` is normal MC. I
currently use `theta = 1.20`. The code does not label a cell as either “MC” or
“minmod.” It checks all three possible slopes in every cell and keeps the
smallest safe one. This gives me a smoother compromise instead of a sudden
switch between two limiters.

## How the reference solver detects a failed time step

The solver normally uses generalized minmod with `theta = 1.20`. After each
SSP-RK3 trial stage, I recover the macroscopic variables and check that:

- density and temperature are finite;
- density is positive; and
- temperature is positive.

If a stage fails, I reject that trial and cut `dt` in half. The retry count
also tells the solver when to use something safer:

```text
retryCount < 4    generalized minmod
retryCount < 8    minmod
otherwise         first-order upwind
```

This switch happens because an **RK trial state became nonphysical**, not just
because the code found a shock. If smaller steps and safer limiters still do
not work, the solver stops and prints the problem instead of silently changing
the answer.

## What the matched limiter test showed

At `Nx = 6250`, all cases used the same physical model, time integrator,
domain, relaxation time, and CFL number:

| Limiter | Global L1 temperature error | Temperature peak excess |
|---|---:|---:|
| Minmod | `2.791803e-4` | `7.964756e-4` |
| MC | `1.851267e-4` | `8.087784e-3` |
| Population hybrid | `2.355076e-4` | `4.650723e-3` |
| Macroscopic hybrid | `2.787301e-4` | `1.005029e-3` |
| Generalized minmod, `theta=1.20` | `2.222085e-4` | `2.048193e-3` |

MC had the smallest average temperature error but the worst overshoot. Minmod
kept the overshoot smallest but smeared the result more. I picked generalized
minmod as the current middle ground. I am not claiming that it is automatically
the best limiter for every problem.

![Limiter accuracy and overshoot tradeoff](results/d1v5_limiter_tradeoff_Nx6250.png)

## Where the implementations are located

- Current solver and retry logic:
  [`matlab/d1v5/d1v5_sodshock_gminmod_rk3.m`](matlab/d1v5/d1v5_sodshock_gminmod_rk3.m)
- Matched five-method comparison:
  [`matlab/d1v5/d1v5_sodshock_limiter_history_comparison.m`](matlab/d1v5/d1v5_sodshock_limiter_history_comparison.m)
- Numerical comparison data:
  [`results/d1v5_limiter_history_Nx6250_metrics.csv`](results/d1v5_limiter_history_Nx6250_metrics.csv)
