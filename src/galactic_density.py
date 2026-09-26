"""Galactic line-of-sight density model for COSMOS-Web UCD counts.

This module intentionally contains no luminosity-function numbers yet.
Those will be added only after source values/units are verified.
"""

from __future__ import annotations
import numpy as np
from scipy.integrate import quad

DEG2_TO_SR = (np.pi / 180.0) ** 2


def galactocentric_cylindrical(d_pc, l_rad, b_rad, r_sun_pc=8300.0, z_sun_pc=27.0):
    x = r_sun_pc - d_pc * np.cos(b_rad) * np.cos(l_rad)
    y = -d_pc * np.cos(b_rad) * np.sin(l_rad)
    z = z_sun_pc + d_pc * np.sin(b_rad)
    return np.hypot(x, y), z


def disk_density(d_pc, l_rad, b_rad, rho_sun, L_pc, H_pc,
                 r_sun_pc=8300.0, z_sun_pc=27.0):
    R, Z = galactocentric_cylindrical(d_pc, l_rad, b_rad, r_sun_pc, z_sun_pc)
    return rho_sun * np.exp(-(R-r_sun_pc)/L_pc) * np.exp(-abs(Z-z_sun_pc)/H_pc)


def halo_density(d_pc, l_rad, b_rad, rho_sun,
                 r_sun_pc=8300.0, z_sun_pc=27.0, q=0.64, n=2.77):
    R, Z = galactocentric_cylindrical(d_pc, l_rad, b_rad, r_sun_pc, z_sun_pc)
    ell = np.sqrt(R**2 + (Z/q)**2)
    return rho_sun * (r_sun_pc/ell)**n


def expected_count_component(dmin_pc, dmax_pc, area_deg2, density_fn):
    """Integrate Omega * rho(d) * d^2 dd for one Galactic component."""
    omega = area_deg2 * DEG2_TO_SR
    val, err = quad(lambda d: density_fn(d) * d*d, dmin_pc, dmax_pc)
    return omega * val, omega * err
