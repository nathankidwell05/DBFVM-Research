# Generalized-Minmod Improvement Report

## Purpose

The goal of the limiter work was to keep the Sod shock-tube solution sharper than standard minmod without recreating the temperature overshoot produced by the full monotonized-central (MC) limiter.

## Development path

### 1. Standard minmod baseline

Standard minmod selects the smallest safe one-sided slope. If the left and right changes disagree in sign, it returns zero.

- Benefit: very resistant to new oscillations near shocks and contacts.
- Drawback: frequently makes the reconstructed solution too flat, which adds numerical diffusion and smears the waves.
- In the generalized-minmod formula, this is the `theta = 1` limit.

### 2. Full MC reconstruction

The MC limiter permits larger slopes than standard minmod and includes the centered slope as a candidate.

- Benefit: sharper reconstruction and less wave smearing.
- Drawback found in this solver: the sharper population reconstruction produced an excessive temperature peak.
- In the generalized-minmod formula, this is the `theta = 2` limit.

### 3. Binary MC/minmod hybrid

The first hybrid attempted to use MC in smooth areas and switch to minmod near a detected discontinuity.

- Benefit: tried to combine MC sharpness with minmod protection.
- Drawback: an abrupt two-choice switch can change the reconstruction strongly from one cell to the next. When the decision is made separately for each DBM population, the five populations can also be reconstructed inconsistently even though their moments combine to produce one temperature.

### 4. Shared macroscopic troubled-cell sensor

Density, pressure, and temperature were used to make one troubled-cell decision shared by all five populations.

- Benefit: the populations switched limiter together, making the treatment more physically consistent.
- Drawback: the switch remained binary. A cell still received either full minmod or full MC, and the MC portions could still create a temperature overshoot.

### 5. Generalized minmod

The binary switch was replaced by the parameterized slope

```text
sigma = minmod(theta*dL, 0.5*(dL+dR), theta*dR),  1 <= theta <= 2.
```

The same formula and value of `theta` are used throughout the reconstruction.

- `theta = 1`: standard minmod; safest and most diffusive.
- `theta = 2`: full MC; sharpest end of this limiter family.
- `1 < theta < 2`: controlled compromise between the two.

For the current test, `theta = 1.2` allows somewhat steeper slopes than standard minmod while remaining much closer to its conservative behavior than to full MC.

## What the three candidates do

- `theta*dL`: limits the slope using the change on the left of the cell.
- `0.5*(dL+dR)`: estimates the centered slope using both neighbors. It is normally selected when the two one-sided slopes have the same sign and similar sizes.
- `theta*dR`: limits the slope using the change on the right of the cell.
- If the candidates do not all have the same sign, the slope becomes zero. This prevents the reconstruction from creating a new local peak or valley.

## Improvements provided by generalized minmod

1. **Less diffusion than standard minmod.** Values of `theta` above one permit larger safe slopes, so smooth gradients and wave transitions are not flattened as aggressively.
2. **More overshoot control than full MC.** Values well below two restrict the reconstructed face values before they reach full-MC sharpness.
3. **No abrupt MC/minmod choice during normal reconstruction.** Sharpness is controlled continuously by one parameter rather than by switching between two separate limiter formulas.
4. **Correct behavior in smooth linear regions.** When the left and right slopes are similar, the centered candidate can be selected and preserve the local gradient.
5. **Automatic protection near extrema.** Opposite slope signs make the limited slope zero, reducing the chance of a new numerical oscillation.
6. **A parameter that can be tested quantitatively.** `theta` can be swept while holding the grid, relaxation time, CFL number, final time, and initial conditions fixed.

## Validation measurements

- **Global L1 error:** average absolute difference between DBM and the exact Sod solution over the entire domain.
- **Wave-region L1 error:** the same average measured only where the rarefaction, contact, and shock are located. This prevents undisturbed constant regions from making the error look artificially small.
- **Temperature L-infinity error:** the single largest absolute temperature difference in the wave region.
- **Temperature peak excess:** `max(T_DBM) - max(T_exact)` in the wave region. A positive value means overshoot; a negative value means the DBM peak is too low.
- **Positive temperature excess:** the peak excess clipped at zero. It reports only actual overshoot.
- **Temperature peak excess percentage:** the signed peak excess divided by the exact peak and multiplied by 100.

## Important limitation

Generalized minmod improves the spatial reconstruction, but it does not by itself guarantee positivity or eliminate every source of error. Results still depend on the KT-D1V5 equilibrium, relaxation time, CFL number, grid spacing, boundary placement, SSP-RK3 time integration, and the exact region used for comparison.
