import jax.numpy as jnp

from dbm2d.config import SodConfig, SolverConfig
from dbm2d.equilibrium import feq_kt_d2v9
from dbm2d.velocity_models import VelocityModel


def make_grid(config: SolverConfig) -> tuple[jnp.ndarray, jnp.ndarray, float, float]:
    dx = config.lx / config.nx
    dy = config.ly / config.ny
    x = (jnp.arange(config.nx) + 0.5) * dx
    y = (jnp.arange(config.ny) + 0.5) * dy
    return x, y, dx, dy


def sod_initial_state(
    solver_config: SolverConfig,
    sod_config: SodConfig,
    model: VelocityModel,
) -> tuple[jnp.ndarray, dict[str, jnp.ndarray]]:
    x, y, _, _ = make_grid(solver_config)
    xx = x[:, None] * jnp.ones((1, solver_config.ny))
    _ = y

    left = xx < sod_config.diaphragm_x
    rho = jnp.where(left, sod_config.rho_left, sod_config.rho_right)
    ux = jnp.where(left, sod_config.u_left, sod_config.u_right)
    uy = jnp.where(left, sod_config.v_left, sod_config.v_right)
    pressure = jnp.where(left, sod_config.p_left, sod_config.p_right)
    temperature = pressure / rho

    f0 = feq_kt_d2v9(rho, ux, uy, temperature, solver_config.gamma, model)
    primitive = {
        "rho": rho,
        "ux": ux,
        "uy": uy,
        "temperature": temperature,
        "pressure": pressure,
    }
    return f0, primitive
