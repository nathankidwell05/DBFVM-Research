import jax.numpy as jnp


def apply_transmissive_boundaries(f: jnp.ndarray) -> jnp.ndarray:
    """Zero-gradient boundaries on all sides.

    This is adequate for first Sod-shock validation when waves do not reach the
    boundaries by `t_end`. Wedge flow will need supersonic inflow/outflow and
    wall treatment.
    """

    f = f.at[0, :, :].set(f[1, :, :])
    f = f.at[-1, :, :].set(f[-2, :, :])
    f = f.at[:, 0, :].set(f[:, 1, :])
    f = f.at[:, -1, :].set(f[:, -2, :])
    return f
