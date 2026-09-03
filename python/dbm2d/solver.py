from dataclasses import dataclass

import jax
import jax.numpy as jnp

from dbm2d.boundary import apply_transmissive_boundaries
from dbm2d.config import SolverConfig
from dbm2d.equilibrium import feq_kt_d2v9
from dbm2d.flux import dbm_advection_divergence
from dbm2d.moments import primitive_from_f
from dbm2d.velocity_models import VelocityModel


@dataclass(frozen=True)
class RunResult:
    f: jnp.ndarray
    time: float
    steps: int
    dt: float


def stable_dt(config: SolverConfig, model: VelocityModel, dx: float, dy: float) -> float:
    """Choose an explicit time step from advection and BGK-relaxation limits."""

    cx_max = float(jnp.max(jnp.abs(model.cx)))
    cy_max = float(jnp.max(jnp.abs(model.cy)))
    dt_adv = config.cfl * min(dx / max(cx_max, 1.0e-14), dy / max(cy_max, 1.0e-14))

    # Explicit source integration becomes fragile when dt >> tau.
    dt_relax = 0.25 * config.tau
    return min(dt_adv, dt_relax)


def rhs(
    f: jnp.ndarray,
    config: SolverConfig,
    model: VelocityModel,
    dx: float,
    dy: float,
) -> jnp.ndarray:
    rho, ux, uy, temperature, _ = primitive_from_f(f, model, config.gamma)
    feq = feq_kt_d2v9(rho, ux, uy, temperature, config.gamma, model)
    div = dbm_advection_divergence(f, model.cx, model.cy, dx, dy, config.limiter)
    return -div + (feq - f) / config.tau


def step_forward(
    f: jnp.ndarray,
    config: SolverConfig,
    model: VelocityModel,
    dx: float,
    dy: float,
    dt: float,
) -> jnp.ndarray:
    """SSP-RK2 update for the discrete Boltzmann equation."""

    f = apply_transmissive_boundaries(f)
    k1 = rhs(f, config, model, dx, dy)
    f1 = apply_transmissive_boundaries(f + dt * k1)
    k2 = rhs(f1, config, model, dx, dy)
    return apply_transmissive_boundaries(0.5 * f + 0.5 * (f1 + dt * k2))


def run(
    f0: jnp.ndarray,
    config: SolverConfig,
    model: VelocityModel,
    dx: float,
    dy: float,
    jit: bool = True,
) -> RunResult:
    dt0 = stable_dt(config, model, dx, dy)
    steps = int(jnp.ceil(config.t_end / dt0))
    dt = config.t_end / steps

    def body(_, f):
        return step_forward(f, config, model, dx, dy, dt)

    if jit:
        body_jit = jax.jit(body)
        f_final = jax.lax.fori_loop(0, steps, body_jit, f0)
    else:
        f_final = f0
        for n in range(steps):
            f_final = body(n, f_final)

    return RunResult(f=f_final, time=config.t_end, steps=steps, dt=float(dt))
