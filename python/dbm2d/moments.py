import jax.numpy as jnp

from dbm2d.velocity_models import VelocityModel


def primitive_from_f(
    f: jnp.ndarray,
    model: VelocityModel,
    gamma: float,
    rho_floor: float = 1.0e-12,
    temperature_floor: float = 1.0e-12,
) -> tuple[jnp.ndarray, jnp.ndarray, jnp.ndarray, jnp.ndarray, jnp.ndarray]:
    """Recover rho, ux, uy, temperature, and pressure from distributions."""

    b = 2.0 / (gamma - 1.0)
    rho = jnp.maximum(jnp.sum(f, axis=-1), rho_floor)
    mx = jnp.sum(f * model.cx, axis=-1)
    my = jnp.sum(f * model.cy, axis=-1)
    ux = mx / rho
    uy = my / rho

    c2_eta2 = model.cx**2 + model.cy**2 + model.eta2
    total_energy_moment = jnp.sum(f * c2_eta2, axis=-1)
    temperature = (total_energy_moment / rho - (ux**2 + uy**2)) / b
    temperature = jnp.maximum(temperature, temperature_floor)
    pressure = rho * temperature
    return rho, ux, uy, temperature, pressure


def nonequilibrium_stress_heat_flux(
    f: jnp.ndarray,
    feq: jnp.ndarray,
    model: VelocityModel,
    ux: jnp.ndarray,
    uy: jnp.ndarray,
) -> dict[str, jnp.ndarray]:
    """Compute simple central nonequilibrium stress and heat-flux moments.

    These are starter diagnostics. They are useful for tracking where the
    distribution departs from local equilibrium, but the final interpretation
    depends on the chosen DBM velocity model and recovered hydrodynamic order.
    """

    df = f - feq
    vx = model.cx - ux[..., None]
    vy = model.cy - uy[..., None]
    v2 = vx**2 + vy**2 + model.eta2

    return {
        "pi_xx": jnp.sum(df * vx * vx, axis=-1),
        "pi_xy": jnp.sum(df * vx * vy, axis=-1),
        "pi_yy": jnp.sum(df * vy * vy, axis=-1),
        "q_x": 0.5 * jnp.sum(df * v2 * vx, axis=-1),
        "q_y": 0.5 * jnp.sum(df * v2 * vy, axis=-1),
    }
