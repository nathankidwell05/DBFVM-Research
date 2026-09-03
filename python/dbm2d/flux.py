import jax.numpy as jnp


def minmod(a: jnp.ndarray, b: jnp.ndarray) -> jnp.ndarray:
    same_sign = jnp.sign(a) == jnp.sign(b)
    return jnp.where(same_sign, jnp.sign(a) * jnp.minimum(jnp.abs(a), jnp.abs(b)), 0.0)


def limited_slope(u: jnp.ndarray, axis: int, limiter: str = "minmod") -> jnp.ndarray:
    """Return a slope-limited reconstruction slope along `axis`."""

    u_left = jnp.roll(u, 1, axis=axis)
    u_right = jnp.roll(u, -1, axis=axis)
    dm = u - u_left
    dp = u_right - u

    if limiter.lower() == "mc":
        return minmod(minmod(2.0 * dm, 2.0 * dp), 0.5 * (dm + dp))
    if limiter.lower() == "minmod":
        return minmod(dm, dp)
    raise ValueError(f"Unknown limiter: {limiter}")


def upwind_flux_divergence(
    f: jnp.ndarray,
    velocity: jnp.ndarray | float,
    dx: float,
    axis: int,
    limiter: str = "minmod",
) -> jnp.ndarray:
    """Finite-volume divergence of `velocity * f` along one physical axis.

    Periodic roll operations are used internally, then boundary rows/columns are
    overwritten by the boundary-condition routine before each step. This first
    prototype uses transmissive behavior at domain edges.
    """

    slope = limited_slope(f, axis=axis, limiter=limiter)

    left_state_at_face = f + 0.5 * slope
    right_state_at_face = jnp.roll(f - 0.5 * slope, -1, axis=axis)
    flux_right = jnp.where(
        velocity >= 0.0,
        velocity * left_state_at_face,
        velocity * right_state_at_face,
    )

    flux_left = jnp.roll(flux_right, 1, axis=axis)
    return (flux_right - flux_left) / dx


def dbm_advection_divergence(
    f: jnp.ndarray,
    cx: jnp.ndarray,
    cy: jnp.ndarray,
    dx: float,
    dy: float,
    limiter: str = "minmod",
) -> jnp.ndarray:
    """Compute sum_i c_ix d_x f_i + c_iy d_y f_i for all populations."""

    pieces = []
    for q in range(int(cx.shape[0])):
        fq = f[..., q]
        div_x = upwind_flux_divergence(fq, cx[q], dx, axis=0, limiter=limiter)
        div_y = upwind_flux_divergence(fq, cy[q], dy, axis=1, limiter=limiter)
        pieces.append(div_x + div_y)
    return jnp.stack(pieces, axis=-1)
