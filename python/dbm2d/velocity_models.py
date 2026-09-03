from dataclasses import dataclass

import jax.numpy as jnp


@dataclass(frozen=True)
class VelocityModel:
    """Discrete molecular velocity model.

    `cx` and `cy` are the translational discrete velocities.
    `eta2` stores extra internal-energy degrees associated with each population.
    """

    name: str
    cx: jnp.ndarray
    cy: jnp.ndarray
    eta2: jnp.ndarray

    @property
    def q(self) -> int:
        return int(self.cx.shape[0])

    @property
    def speed_max(self) -> float:
        return float(jnp.max(jnp.sqrt(self.cx**2 + self.cy**2)))


def kt_d2v9(v1: float = 1.0, v2: float = 3.0, eta0: float = 2.0) -> VelocityModel:
    """Kataoka-Tsutahara-style D2V9 model used as the first prototype.

    The diagonal speed magnitude is `v2`; each diagonal component is
    `v2 / sqrt(2)`. A nonzero internal-energy parameter is assigned only to the
    rest particle, following the common compact thermal D2V9 construction.

    This model is a useful starter for Sod-shock validation. It should not be
    treated as the final model for wedge-flow nonequilibrium studies.
    """

    inv_sqrt2 = 1.0 / jnp.sqrt(2.0)
    cx = jnp.array(
        [
            0.0,
            v1,
            -v1,
            0.0,
            0.0,
            v2 * inv_sqrt2,
            -v2 * inv_sqrt2,
            v2 * inv_sqrt2,
            -v2 * inv_sqrt2,
        ]
    )
    cy = jnp.array(
        [
            0.0,
            0.0,
            0.0,
            v1,
            -v1,
            v2 * inv_sqrt2,
            v2 * inv_sqrt2,
            -v2 * inv_sqrt2,
            -v2 * inv_sqrt2,
        ]
    )
    eta2 = jnp.array([eta0**2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    return VelocityModel(name="KT-D2V9", cx=cx, cy=cy, eta2=eta2)
