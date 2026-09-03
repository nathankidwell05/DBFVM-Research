from dataclasses import dataclass


@dataclass(frozen=True)
class SolverConfig:
    """Numerical parameters for the finite-volume DBM solver."""

    nx: int = 400
    ny: int = 8
    lx: float = 1.0
    ly: float = 0.02
    gamma: float = 1.4
    tau: float = 1.0e-4
    cfl: float = 0.35
    t_end: float = 0.15
    limiter: str = "minmod"


@dataclass(frozen=True)
class SodConfig:
    """Left/right primitive states for Sod's shock-tube problem."""

    rho_left: float = 1.0
    u_left: float = 0.0
    v_left: float = 0.0
    p_left: float = 1.0
    rho_right: float = 0.125
    u_right: float = 0.0
    v_right: float = 0.0
    p_right: float = 0.1
    diaphragm_x: float = 0.5
