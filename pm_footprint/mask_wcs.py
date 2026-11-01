"""FITS PC/TAN WCS helpers for COSMOS-Web JWST star masks.

FITS convention: CD[i,j] = CDELT[i] * PC[i,j], i.e. row scaling.
Scaling columns instead changes the off-diagonal terms if axis signs differ.
Input coordinates are zero-indexed pixel centers, as in NumPy.
"""
import numpy as np

RAD = np.pi / 180.0


def cd_from_header(header):
    """Return the 2x2 degree-per-pixel matrix from explicit PC/CDELT."""
    pc = np.array([[float(header[f"PC{i}_{j}"]) for j in (1, 2)]
                   for i in (1, 2)], dtype=float)
    scales = np.array([float(header[f"CDELT{i}"]) for i in (1, 2)])
    cd = scales[:, None] * pc
    if not np.isfinite(cd).all() or abs(np.linalg.det(cd)) < 1e-25:
        raise ValueError("Invalid or singular FITS PC/CDELT matrix")
    return cd


def pixel_to_sky(x, y, header):
    """Inverse gnomonic TAN projection, returning RA and Dec in radians."""
    x, y = np.broadcast_arrays(np.asarray(x, dtype=float),
                               np.asarray(y, dtype=float))
    cd = cd_from_header(header)
    dx = x + 1 - float(header["CRPIX1"])
    dy = y + 1 - float(header["CRPIX2"])
    xi = (cd[0, 0] * dx + cd[0, 1] * dy) * RAD
    eta = (cd[1, 0] * dx + cd[1, 1] * dy) * RAD
    a = float(header["CRVAL1"]) * RAD
    d = float(header["CRVAL2"]) * RAD
    sd, cd0 = np.sin(d), np.cos(d)
    ra = a + np.arctan2(xi, cd0 - eta * sd)
    dec = np.arctan2(sd + eta * cd0, np.hypot(cd0 - eta * sd, xi))
    return ra, dec


def pixel_to_lambert(x, y, header, center_ra_deg=150.12,
                     center_dec_deg=2.2):
    """Lambert azimuthal equal-area sky coordinates in degrees."""
    ra, dec = pixel_to_sky(x, y, header)
    a0, d0 = center_ra_deg * RAD, center_dec_deg * RAD
    delta = (ra - a0 + np.pi) % (2 * np.pi) - np.pi
    sd0, cd0 = np.sin(d0), np.cos(d0)
    k = np.sqrt(2 / (1 + sd0 * np.sin(dec)
                     + cd0 * np.cos(dec) * np.cos(delta)))
    return (k * np.cos(dec) * np.sin(delta) / RAD,
            k * (cd0 * np.sin(dec)
                 - sd0 * np.cos(dec) * np.cos(delta)) / RAD)


def tan_pixel_area_deg2(x, y, header):
    """Local TAN solid-angle Jacobian in square degrees per native pixel."""
    cd = cd_from_header(header)
    dx = np.asarray(x, dtype=float) + 1 - float(header["CRPIX1"])
    dy = np.asarray(y, dtype=float) + 1 - float(header["CRPIX2"])
    xi = (cd[0, 0] * dx + cd[0, 1] * dy) * RAD
    eta = (cd[1, 0] * dx + cd[1, 1] * dy) * RAD
    return abs(np.linalg.det(cd)) / (1 + xi * xi + eta * eta) ** 1.5


def legacy_column_scaling_error_arcsec(header, x, y):
    """Exact tangent-plane shift caused by the former column-scaled PC code.

    This is a diagnostic of coordinate error, NOT a sky-area correction.
    """
    dx = np.asarray(x, dtype=float) + 1 - float(header["CRPIX1"])
    dy = np.asarray(y, dtype=float) + 1 - float(header["CRPIX2"])
    s1, s2 = float(header["CDELT1"]), float(header["CDELT2"])
    e1 = (s2 - s1) * float(header["PC1_2"]) * dy
    e2 = (s1 - s2) * float(header["PC2_1"]) * dx
    return np.hypot(e1, e2) * 3600.0
