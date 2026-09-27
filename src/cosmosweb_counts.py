"""COSMOS-Web brown-dwarf count model.

Implements line-of-sight integration by Teff bin. The default local densities
are the observational 150-K bins from Kirkpatrick et al. (2024, ApJS 271, 55,
Table 17), while synthetic JWST reach is supplied separately from FC-ENZO
Elf Owl photometry.

Population density, survey depth, and selection function remain separate.
Kirkpatrick local densities are treated as the TOTAL local population by
default, then partitioned into thin disk, thick disk, and halo components.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import quad
from astropy.coordinates import SkyCoord
import astropy.units as u

TEFF_DENSITY_2024 = {
    (300, 450): (3.09e-3, None, "lower_limit"),
    (450, 600): (4.20e-3, 0.67e-3, "measurement"),
    (600, 750): (3.10e-3, 0.37e-3, "measurement"),
    (750, 900): (1.99e-3, 0.31e-3, "measurement"),
    (900, 1050): (1.58e-3, 0.27e-3, "measurement"),
    (1050, 1200): (1.38e-3, 0.26e-3, "measurement"),
    (1200, 1350): (2.12e-3, 0.30e-3, "measurement"),
    (1350, 1500): (1.04e-3, 0.23e-3, "measurement"),
    (1500, 1650): (0.78e-3, 0.19e-3, "measurement"),
    (1650, 1800): (0.72e-3, 0.19e-3, "measurement"),
    (1800, 1950): (0.48e-3, 0.16e-3, "measurement"),
    (1950, 2100): (0.81e-3, 0.19e-3, "measurement"),
    (2100, 2250): (0.36e-3, None, "lower_limit"),
}


def galactic_lb(ra_deg=150.12, dec_deg=2.20):
    c = SkyCoord(ra_deg * u.deg, dec_deg * u.deg, frame="icrs").galactic
    return c.l.radian, c.b.radian


def disk_rho(d, l, b, rho0, L, H, R0=8300.0, Z0=27.0):
    X = R0 - d * np.cos(b) * np.cos(l)
    Y = -d * np.cos(b) * np.sin(l)
    Z = Z0 + d * np.sin(b)
    R = np.hypot(X, Y)
    return rho0 * np.exp(-(R - R0) / L) * np.exp(-np.abs(Z - Z0) / H)


def halo_rho(d, l, b, rho0, R0=8300.0, Z0=27.0, q=0.64, n=2.77):
    X = R0 - d * np.cos(b) * np.cos(l)
    Y = -d * np.cos(b) * np.sin(l)
    Z = Z0 + d * np.sin(b)
    R = np.hypot(X, Y)
    ell = np.sqrt(R**2 + (Z / q) ** 2)
    ell0 = np.sqrt(R0**2 + (Z0 / q) ** 2)
    return rho0 * (ell0 / ell) ** n


def local_component_densities(
    rho_local_total,
    thick_fraction=0.02,
    halo_fraction=0.0025,
):
    """Partition a measured total local density without double counting.

    Fractions are thick/thin and halo/thin ratios. Returned components sum
    exactly to rho_local_total.
    """
    thin = rho_local_total / (1.0 + thick_fraction + halo_fraction)
    return thin, thin * thick_fraction, thin * halo_fraction


def count_to_distance(
    dmax_pc,
    area_arcmin2,
    rho0,
    H_pc=187.0,
    thick_fraction=0.02,
    halo_fraction=0.0025,
    rho0_is_total=True,
):
    """Expected count to dmax for a COSMOS-like sightline.

    By default rho0 is the measured TOTAL local density. Set rho0_is_total=False
    only for a legacy sensitivity test in which rho0 explicitly anchors the
    thin disk and thick/halo populations are added on top.
    """
    l, b = galactic_lb()
    omega = (area_arcmin2 / 3600.0) * (np.pi / 180.0) ** 2

    if rho0_is_total:
        rho_thin, rho_thick, rho_halo = local_component_densities(
            rho0, thick_fraction, halo_fraction
        )
    else:
        rho_thin = rho0
        rho_thick = rho0 * thick_fraction
        rho_halo = rho0 * halo_fraction

    def integrand(d):
        return (
            disk_rho(d, l, b, rho_thin, 2600.0, H_pc)
            + disk_rho(d, l, b, rho_thick, 3600.0, 900.0)
            + halo_rho(d, l, b, rho_halo)
        ) * d * d

    return omega * quad(integrand, 0.0, dmax_pc, epsabs=1e-9)[0]


def distance_at_mag(d25_pc, mag):
    return d25_pc * 10 ** (0.2 * (mag - 25.0))
