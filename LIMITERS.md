# Limiter and Stability Guide

This is my plain-language guide to every limiter I tested in the KT-D1V5
solver. The main questions are simple: what does the limiter look at, when
does it become more cautious, and what does it do badly?

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

## Where the limiter fits into my code

The limiter is only one part of the finite-volume transport calculation. The
full path through `build_MUSCL_flux_divergence` is:

1. `fq = f(:,q)` pulls one of the five populations out across the whole tube.
2. `dL` and `dR` measure how that population changes on either side of each
   interior cell.
3. The selected limiter turns those two differences into one safe cell slope,
   `sigma`.
4. MUSCL uses `sigma` to estimate the population on both faces of the cell.
5. Upwinding chooses the face value that is actually moving into the face.
6. The east-face flux minus the west-face flux gives the net transport out of
   the cell.
7. The code combines that transport with the BGK collision term and advances
   the populations with SSP-RK3.

The important lines for face reconstruction are

```matlab
fLeftFace  = fq(leftCells)  + 0.5*sigma(leftCells);
fRightFace = fq(rightCells) - 0.5*sigma(rightCells);
```

If `sigma` is large, the face values are allowed to be farther from the cell
average. This gives a sharper result, but it also gives the reconstruction
more room to overshoot. If `sigma` is zero, both reconstructed sides fall back
to the original cell value. That is safer, but more diffusive.

The limiter works on the populations `f_i`, not directly on density, velocity,
pressure, or temperature. This is because the kinetic equation transports the
five populations. I recover the macroscopic variables afterward by taking
moments of all five populations together. This is also why inconsistent
limiting among the populations can show up as an error in recovered
temperature.

### Limiting and upwinding are not the same operation

The limiter decides **how steep the reconstructed line inside a cell may be**.
Upwinding then decides **which side of a face supplies the flux**:

```matlab
if c(q) > 0
    Fface = c(q)*fLeftFace;      % population moves left to right
elseif c(q) < 0
    Fface = c(q)*fRightFace;     % population moves right to left
end
```

Changing the limiter does not change the direction of transport. It only
changes the reconstructed value supplied to the same upwind flux.

## A small numerical example

Suppose one population has these values in three neighboring cells:

```text
f(i-1) = 0.20,   f(i) = 0.40,   f(i+1) = 0.50
```

Then `dL = 0.20` and `dR = 0.10`. Both are positive, so the population keeps
rising across the cell and a nonzero slope is allowed.

| Limiter | Candidates | Chosen slope |
|---|---|---:|
| First order | No candidates; slope is forced off | `0.00` |
| Minmod | `0.20`, `0.10` | `0.10` |
| MC | `0.40`, `0.15`, `0.20` | `0.15` |
| Generalized minmod, `theta=1.20` | `0.24`, `0.15`, `0.12` | `0.12` |

This shows the basic tradeoff. MC keeps the steepest slope, minmod keeps the
smallest, and generalized minmod sits between them. If `dL = 0.20` and
`dR = -0.10`, the signs disagree and all three MUSCL limiters return zero.
That prevents the reconstructed line from creating a new peak inside the cell.

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

The code completes that idea with

```matlab
if limiterType == "firstorder"
    sigma(:) = 0.0;
end
```

There is no smooth-cell or troubled-cell decision here. Every cell receives
the same treatment. Its main failure is excessive numerical diffusion: the
method is stable partly because it erases sharp information.

### Minmod

I think of minmod as asking two questions:

1. Do the left and right differences point in the same direction?
2. If so, which has the smaller magnitude?

If both answers are acceptable, minmod keeps the smaller slope. If the signs
disagree, it returns zero. That makes it hard to destabilize, but it also
flattens real gradients and spreads the shock or contact over more cells.

The vectorized MATLAB version first assumes every slope is zero:

```matlab
s = zeros(size(dL));
positive = (dL > 0) & (dR > 0);
negative = (dL < 0) & (dR < 0);
```

`positive` and `negative` are logical masks. They identify the cells where
both differences agree in sign. The code then keeps the smaller magnitude:

```matlab
s(positive) = min(dL(positive),dR(positive));
s(negative) = max(dL(negative),dR(negative));
```

For two negative numbers, `max` selects the one closer to zero, which is still
the smaller-magnitude slope. Cells not included in either mask stay at zero.

Minmod does not know whether it found a physical shock. It only knows that the
local differences either agree or disagree. A smooth local maximum will also
be flattened because its two differences have opposite signs.

### Monotonized central (MC)

MC compares three candidates:

```text
2*dL,   (dL+dR)/2,   2*dR
```

It keeps the smallest candidate only when all three have the same sign. The
centered choice lets MC keep a steeper slope than minmod. That made the Sod
result sharper and lowered the average errors, but it also caused the largest
temperature overshoot.

The code builds the three candidates with

```matlab
a = 2*dL;
b = 0.5*(dL+dR);
c = 2*dR;
```

It then creates masks requiring `a`, `b`, and `c` to all be positive or all be
negative. If they agree, the code takes the smallest magnitude. If they do not
agree, the slope stays zero. MC therefore still has a safety rule, but it
allows a larger safe slope than minmod in many cells.

MC's weakness in this solver was not a complete numerical crash. The run could
remain stable while the sharper population reconstruction produced a local
temperature value above the exact plateau. This is why checking only whether
the code runs is not enough; I also measure the temperature peak excess.

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

In plain terms, `curvatureRatio` becomes large when the slope entering a cell
looks very different from the slope leaving it. `slopeSizeRatio` becomes large
when one side is nearly flat but the other side changes sharply. The tiny
number called `small` only prevents division by zero in a flat region.

The actual switch is

```matlab
useMinmod = curvatureRatio > 0.35 | slopeSizeRatio > 3.0;
slope = slopeMC;
slope(useMinmod) = slopeMinmod(useMinmod);
```

This block runs inside the loop over `q`, so every population builds its own
`useMinmod` mask. That independence is the central shortcoming of this version.
The numbers `0.35` and `3.0` are also chosen sensor thresholds, not universal
physical constants.

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

The code loops over the three macroscopic fields and combines their warnings
with a logical OR:

```matlab
fields = {rho,p,T};
troubledCell(2:Nx-1) = troubledCell(2:Nx-1) | fieldTrouble;
```

This means density can flag a contact, pressure can flag a shock, and
temperature can flag a thermal jump. Once any field marks a physical cell,
all five populations use minmod there. The code also expands the mask by one
cell in both directions because a reconstructed face uses information from
neighboring cells.

This sensor is more physically consistent than the population sensor, but it
still has two weaknesses. First, its cutoff values were selected for this Sod
test and may need to change for another problem. Second, a cell jumps suddenly
from full MC to full minmod when it crosses a threshold. A tiny change in the
sensor can therefore cause a much larger change in the reconstructed slope.

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

The code performs this directly:

```matlab
leftCandidate   = theta*dL;
centerCandidate = 0.5*(dL+dR);
rightCandidate  = theta*dR;

allPositive = leftCandidate > 0 & centerCandidate > 0 & rightCandidate > 0;
allNegative = leftCandidate < 0 & centerCandidate < 0 & rightCandidate < 0;
```

It keeps the smallest positive candidate or the negative candidate closest to
zero. Everything else stays at zero. The middle candidate is chosen when the
centered slope is smaller than both scaled one-sided candidates. For example,
if `dL = dR = 0.10` and `theta = 1.20`, the candidates are `0.12`, `0.10`, and
`0.12`, so the middle value is selected.

Generalized minmod is not dynamically changing `theta` from cell to cell.
`theta` remains `1.20` for the whole normal run. The local result still adapts
because the smallest of the three candidates can be different in each cell.
Its weakness is that one fixed value of `theta` may be too cautious in a
smooth region but too aggressive near a different kind of discontinuity.

## What “switching limiters” means in this project

There are three different ideas that can sound like the same thing:

1. **A limiter returning zero:** minmod, MC, and generalized minmod all turn
   off the local slope when their candidates disagree in sign. The named
   limiter has not changed.
2. **A hybrid sensor switching reconstruction:** the two hybrid experiments
   explicitly choose MC or minmod based on a local warning rule.
3. **The retry system changing the active limiter:** the reference solver
   changes from generalized minmod to minmod and then first order only after
   an SSP-RK3 trial produces a nonphysical state.

The current normal calculation uses the first idea. The third idea is only a
safety fallback. It does not constantly switch between all the tested methods
as the waves move through the tube.

## How the reference solver detects a failed time step

The solver normally uses generalized minmod with `theta = 1.20`. After each
SSP-RK3 trial stage, I recover the macroscopic variables and check that:

- density and temperature are finite;
- density is positive; and
- temperature is positive.

The function `is_physical_D1V5` carries out that check:

```matlab
[rho,~,T,~] = recover_macros_D1V5(f,b,c,h2);
ok = all(isfinite(rho(:))) && all(isfinite(T(:))) ...
    && all(rho(:) > 0) && all(T(:) > 0);
```

The tildes mean that this check ignores the returned velocity and pressure.
The current pass/fail decision is specifically based on density and
temperature. This catches a state that the model cannot physically use, but
it does not catch every inaccurate answer. A small positive temperature
overshoot, for example, still passes this test. I have to catch that separately
by comparing with the exact solution and measuring peak excess.

If a stage fails, I reject that trial and cut `dt` in half. The retry count
also tells the solver when to use something safer:

```text
retryCount < 4    generalized minmod
retryCount < 8    minmod
otherwise         first-order upwind
```

The code starts a new time step with `retryCount = 0`. If a trial fails, it
runs

```matlab
retryCount = retryCount + 1;
dt = 0.5*dt;
continue;
```

`continue` sends the program back to the start of the retry loop without
accepting the bad trial state. Retry counts `0` through `3` use generalized
minmod, `4` through `7` use minmod, and `8` or higher use first-order upwind.
The documented `Nx=6250` result had zero retries, so that result used
generalized minmod throughout rather than quietly mixing in the fallback
methods.

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
