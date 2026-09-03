# Research sources

## Core model and implementation sources

1. Kataoka, T., & Tsutahara, M. (2004). *Lattice Boltzmann method for the
   compressible Euler equations*. Physical Review E, 69, 056702.
   https://doi.org/10.1103/PhysRevE.69.056702
   - Primary source for the D1Q5/D1V5 KT equilibrium form and coefficient
     construction used in the current implementation.

2. Kataoka, T., & Tsutahara, M. (2004). *Lattice Boltzmann model for the
   compressible Navier–Stokes equations with flexible specific-heat ratio*.
   Physical Review E, 69, 035701(R).
   https://doi.org/10.1103/PhysRevE.69.035701
   - Related source for flexible specific heat and internal degrees of freedom.

3. Xu et al. *Discrete Boltzmann Modeling of Compressible Flows*.
   https://arxiv.org/abs/1708.09187
   - Foundational DBM overview, hydrodynamic recovery, and TNE framing.

4. Gan, Y., Xu, A., Zhang, G., & Lai, H. (2018). *Three-dimensional discrete
   Boltzmann models for compressible flows in and out of equilibrium*.
   https://doi.org/10.1177/0954406217742181
   - Euler/Navier–Stokes-level construction and multidimensional validation.

5. Gan et al. *Discrete Boltzmann trans-scale modeling of high-speed
   compressible flows*.
   https://arxiv.org/abs/1801.04522
   - D2V26 Burnett-level model for stronger nonequilibrium.

6. Zhang et al. *Discrete Boltzmann modeling of high-speed compressible flows
   with various depths of non-equilibrium*.
   https://arxiv.org/abs/2205.13809
   - Framework for matching DBM complexity to nonequilibrium depth.

7. Lin, C., Sun, X., Su, X., Lai, H., & Fang, X. (2023). *A discrete Boltzmann
   model with symmetric velocity discretization for compressible flow*.
   https://doi.org/10.1088/1674-1056/acea6b
   - Symmetric velocity discretization and moment-matching construction.

8. Lin, Su, Fei, & Luo. *Central-moment-based discrete Boltzmann modeling of
   compressible flows*.
   https://doi.org/10.1142/S0129183126500397
   - Central-moment alternative; includes shock-tube and flow benchmarks.

9. Zhang et al. *Discrete Boltzmann method for non-equilibrium flows based on
   the Shakhov model*.
   https://doi.org/10.1016/j.cpc.2018.12.018
   - Adjustable Prandtl number and broader continuum-to-transition regimes.

## Later-stage high-speed and reactive references

10. Guo et al. *Thermodynamic nonequilibrium effects in three-dimensional
    high-speed compressible flows*.
    https://arxiv.org/abs/2502.01446

11. Lin et al. *MRT discrete Boltzmann method for compressible exothermic
    reactive flows*.
    https://doi.org/10.1016/j.compfluid.2018.02.012

12. Lin et al. *Burnett-level multi-relaxation-time central-moment discrete
    Boltzmann modeling of reactive flows*.
    https://doi.org/10.1016/j.combustflame.2025.114481

## Flux-scheme comparison references

13. *Accuracy Assessment of Upwind Algorithms for Steady-state Computations*.
    https://doi.org/10.1016/S0045-7930(98)00011-5
    - Compares several compressible-flow upwind approaches in a common
      accuracy study, including Steger–Warming, van Leer, Roe, Osher, AUSM,
      and hybrid upwind splitting.

14. *A Review on the Numerical Solution of the 1D Euler Equations*.
    https://eprints.maths.manchester.ac.uk/150/
    - Review-oriented comparison of high-resolution Euler solvers, including
      Roe, HLLE, and AUSM+ variants, with shock-tube-relevant tests.

15. Toro, E. F. *Riemann Solvers and Numerical Methods for Fluid Dynamics*.
    https://link.springer.com/book/10.1007/b79761
    - Reference for exact and approximate Riemann solvers, HLL-family methods,
      Roe, MUSCL, TVD reconstruction, and limiter concepts.

16. Kurganov, A., & Tadmor, E. (2000). *New high-resolution central schemes for
    nonlinear conservation laws and convection–diffusion equations*.
    https://doi.org/10.1006/jcph.2000.6459
    - Published source for the parameterized generalized-minmod slope with
      `theta` between the minmod and MC limits. This project adapts that
      reconstruction to KT-D1V5 populations.

17. Sweby, P. K. (1984). *High resolution schemes using flux limiters for
    hyperbolic conservation laws*.
    https://doi.org/10.1137/0721062
    - Foundational TVD flux-limiter analysis underlying the
      accuracy-versus-monotonicity tradeoff.

## Validation benchmark

18. Sod, G. A. (1978). *A survey of several finite difference methods for
    systems of nonlinear hyperbolic conservation laws*. Journal of
    Computational Physics, 27(1), 1–31.
    https://doi.org/10.1016/0021-9991(78)90023-2
    - Original shock-tube benchmark used for the exact Euler comparison.

## Source notes

- “DBM” is used inconsistently across kinetic-method communities. This list
  follows the Xu–Gan–Zhang–Lin DBM lineage plus the KT model used by the
  working reference implementation.
- Generic thermal or compressible LBM papers should only be added when their
  relationship to this finite-volume DBM program is stated.
- Equations and coefficients should be verified against the paper before a
  paraphrase is treated as implementation authority.
- The Euler flux references provide comparison concepts; they do not by
  themselves justify inserting a conservative-variable flux directly into the
  population-based DBM transport equation.
- The generalized-minmod sources establish the limiter for conservation-law
  numerics. They do not establish that this exact KT-D1V5 population
  reconstruction has previously appeared in DBM literature.
