import jax.numpy as jnp

from dbm2d.velocity_models import VelocityModel


def feq_kt_d2v9(
    rho: jnp.ndarray,
    ux: jnp.ndarray,
    uy: jnp.ndarray,
    temperature: jnp.ndarray,
    gamma: float,
    model: VelocityModel,
    v1: float = 1.0,
    v2: float = 3.0,
    eta0: float = 2.0,
) -> jnp.ndarray:
    """Compressible D2V9 equilibrium used by the starter prototype.

    The expression follows the compact A/B/D polynomial equilibrium used in
    many KT-style compressible lattice/discrete Boltzmann shock-tube examples.
    The output has shape `(nx, ny, 9)`.
    """

    b = 2.0 / (gamma - 1.0)
    u2 = ux**2 + uy**2

    a_rest = ((b - 2.0) / eta0**2) * temperature
    a_card = (
        (((b - 2.0) * v2**2 / eta0**2 + 2.0) * temperature)
        + (v2**2 / v1**2) * u2
        - v2**2
    ) / (4.0 * (v1**2 - v2**2))
    a_diag = (
        (((b - 2.0) * v1**2 / eta0**2 + 2.0) * temperature)
        + (v1**2 / v2**2) * u2
        - v1**2
    ) / (4.0 * (v2**2 - v1**2))

    b_card = (-v2**2 + (b + 2.0) * temperature + u2) / (
        2.0 * v1**2 * (v1**2 - v2**2)
    )
    b_diag = (-v1**2 + (b + 2.0) * temperature + u2) / (
        2.0 * v2**2 * (v2**2 - v1**2)
    )

    d_card = 1.0 / (2.0 * v1**4)
    d_diag = 1.0 / (2.0 * v2**4)

    s = ux[..., None] * model.cx + uy[..., None] * model.cy

    a = jnp.stack(
        [a_rest, a_card, a_card, a_card, a_card, a_diag, a_diag, a_diag, a_diag],
        axis=-1,
    )
    bcoef = jnp.stack(
        [
            jnp.zeros_like(rho),
            b_card,
            b_card,
            b_card,
            b_card,
            b_diag,
            b_diag,
            b_diag,
            b_diag,
        ],
        axis=-1,
    )
    dcoef = jnp.array([0.0, d_card, d_card, d_card, d_card, d_diag, d_diag, d_diag, d_diag])
    return rho[..., None] * (a + bcoef * s + dcoef * s**2)
