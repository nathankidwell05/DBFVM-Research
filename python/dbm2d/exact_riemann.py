from __future__ import annotations

import math

import numpy as np


def _pressure_function(p: float, rho_k: float, p_k: float, a_k: float, gamma: float) -> tuple[float, float]:
    if p > p_k:
        a = 2.0 / ((gamma + 1.0) * rho_k)
        b = (gamma - 1.0) / (gamma + 1.0) * p_k
        root = math.sqrt(a / (p + b))
        f = (p - p_k) * root
        df = root * (1.0 - 0.5 * (p - p_k) / (p + b))
        return f, df

    exponent = (gamma - 1.0) / (2.0 * gamma)
    f = (2.0 * a_k / (gamma - 1.0)) * ((p / p_k) ** exponent - 1.0)
    df = (1.0 / (rho_k * a_k)) * (p / p_k) ** (-(gamma + 1.0) / (2.0 * gamma))
    return f, df


def sod_exact(
    x: np.ndarray,
    t: float,
    x0: float,
    left: tuple[float, float, float],
    right: tuple[float, float, float],
    gamma: float = 1.4,
) -> dict[str, np.ndarray]:
    """Exact Euler solution for Sod's shock tube.

    `left` and `right` are `(rho, u, p)`.
    """

    x = np.asarray(x)
    rho_l, u_l, p_l = left
    rho_r, u_r, p_r = right

    if t <= 0.0:
        mask = x < x0
        rho = np.where(mask, rho_l, rho_r)
        u = np.where(mask, u_l, u_r)
        p = np.where(mask, p_l, p_r)
        return {"rho": rho, "u": u, "p": p, "temperature": p / rho}

    a_l = math.sqrt(gamma * p_l / rho_l)
    a_r = math.sqrt(gamma * p_r / rho_r)

    p_star = max(1.0e-12, 0.5 * (p_l + p_r) - 0.125 * (u_r - u_l) * (rho_l + rho_r) * (a_l + a_r))
    for _ in range(80):
        f_l, df_l = _pressure_function(p_star, rho_l, p_l, a_l, gamma)
        f_r, df_r = _pressure_function(p_star, rho_r, p_r, a_r, gamma)
        update = -(f_l + f_r + u_r - u_l) / (df_l + df_r)
        p_star = max(1.0e-12, p_star + update)
        if abs(update) / (p_star + 1.0e-12) < 1.0e-12:
            break

    f_l, _ = _pressure_function(p_star, rho_l, p_l, a_l, gamma)
    f_r, _ = _pressure_function(p_star, rho_r, p_r, a_r, gamma)
    u_star = 0.5 * (u_l + u_r + f_r - f_l)

    rho = np.empty_like(x, dtype=float)
    u = np.empty_like(x, dtype=float)
    p = np.empty_like(x, dtype=float)
    s = (x - x0) / t

    left_of_contact = s <= u_star

    # Left wave.
    if p_star > p_l:
        alpha = (gamma - 1.0) / (gamma + 1.0)
        rho_star_l = rho_l * ((p_star / p_l + alpha) / (alpha * p_star / p_l + 1.0))
        shock_speed_l = u_l - a_l * math.sqrt((gamma + 1.0) / (2.0 * gamma) * p_star / p_l + (gamma - 1.0) / (2.0 * gamma))
        region = left_of_contact & (s >= shock_speed_l)
        original = left_of_contact & ~region
        rho[region], u[region], p[region] = rho_star_l, u_star, p_star
        rho[original], u[original], p[original] = rho_l, u_l, p_l
    else:
        rho_star_l = rho_l * (p_star / p_l) ** (1.0 / gamma)
        a_star_l = a_l * (p_star / p_l) ** ((gamma - 1.0) / (2.0 * gamma))
        head = u_l - a_l
        tail = u_star - a_star_l
        original = left_of_contact & (s < head)
        star = left_of_contact & (s > tail)
        fan = left_of_contact & ~(original | star)
        rho[original], u[original], p[original] = rho_l, u_l, p_l
        rho[star], u[star], p[star] = rho_star_l, u_star, p_star
        a_fan = (gamma - 1.0) / (gamma + 1.0) * (u_l + 2.0 * a_l / (gamma - 1.0) - s[fan])
        u_fan = s[fan] + a_fan
        rho[fan] = rho_l * (a_fan / a_l) ** (2.0 / (gamma - 1.0))
        u[fan] = u_fan
        p[fan] = p_l * (a_fan / a_l) ** (2.0 * gamma / (gamma - 1.0))

    # Right wave.
    right_of_contact = ~left_of_contact
    if p_star > p_r:
        alpha = (gamma - 1.0) / (gamma + 1.0)
        rho_star_r = rho_r * ((p_star / p_r + alpha) / (alpha * p_star / p_r + 1.0))
        shock_speed_r = u_r + a_r * math.sqrt((gamma + 1.0) / (2.0 * gamma) * p_star / p_r + (gamma - 1.0) / (2.0 * gamma))
        region = right_of_contact & (s <= shock_speed_r)
        original = right_of_contact & ~region
        rho[region], u[region], p[region] = rho_star_r, u_star, p_star
        rho[original], u[original], p[original] = rho_r, u_r, p_r
    else:
        rho_star_r = rho_r * (p_star / p_r) ** (1.0 / gamma)
        a_star_r = a_r * (p_star / p_r) ** ((gamma - 1.0) / (2.0 * gamma))
        head = u_r + a_r
        tail = u_star + a_star_r
        original = right_of_contact & (s > head)
        star = right_of_contact & (s < tail)
        fan = right_of_contact & ~(original | star)
        rho[original], u[original], p[original] = rho_r, u_r, p_r
        rho[star], u[star], p[star] = rho_star_r, u_star, p_star
        a_fan = (gamma - 1.0) / (gamma + 1.0) * (s[fan] - u_r + 2.0 * a_r / (gamma - 1.0))
        u_fan = s[fan] - a_fan
        rho[fan] = rho_r * (a_fan / a_r) ** (2.0 / (gamma - 1.0))
        u[fan] = u_fan
        p[fan] = p_r * (a_fan / a_r) ** (2.0 * gamma / (gamma - 1.0))

    return {"rho": rho, "u": u, "p": p, "temperature": p / rho}
