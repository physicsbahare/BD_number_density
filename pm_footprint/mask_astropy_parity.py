"""Independent FITS WCS parity audit for the twenty COSMOS2025 starmasks.

The production rasterizer uses a lightweight analytic inverse TAN transform.
Before accepting any new union-area product, compare it with Astropy WCS
using the REAL per-tile FITS header, including rotated PC/CDELT axes.
No remote data are fetched; no FITS image array is opened by this helper.
"""
from __future__ import annotations

import numpy as np
from .mask_wcs import cd_from_header, pixel_to_sky, tan_pixel_area_deg2

ARCSEC_PER_RAD = 180.0 * 3600.0 / np.pi


def angular_distance_arcsec(ra1, dec1, ra2, dec2):
    """Numerically stable great-circle distance, all arguments in radians."""
    dra = (np.asarray(ra1) - np.asarray(ra2) + np.pi) % (2*np.pi) - np.pi
    ddec = np.asarray(dec1) - np.asarray(dec2)
    a = np.sin(ddec/2)**2 + np.cos(dec1)*np.cos(dec2)*np.sin(dra/2)**2
    return 2*np.arcsin(np.sqrt(np.clip(a, 0, 1))) * ARCSEC_PER_RAD


def validate_fits_header_wcs(header, nx, ny, tolerance_arcsec=0.001):
    """Verify analytical TAN WCS against independent Astropy at 9 locations.

    A numerical failure stops publication of the union mask; this audit
    does not by itself validate mask-code semantics or exact sky union.
    """
    from astropy.io import fits
    from astropy.wcs import WCS

    if int(nx) < 2 or int(ny) < 2:
        raise ValueError("Invalid native mask dimensions")
    if str(header.get("CTYPE1", "")).strip() != "RA---TAN":
        raise ValueError("Analytic raster only supports undistorted RA---TAN")
    if str(header.get("CTYPE2", "")).strip() != "DEC--TAN":
        raise ValueError("Analytic raster only supports undistorted DEC--TAN")
    if any(k.startswith(("A_", "B_", "AP_", "BP_")) for k in header):
        raise ValueError("SIP distortion keys require a distortion-aware raster")

    keys = ("NAXIS", "NAXIS1", "NAXIS2", "WCSAXES",
            "CTYPE1", "CTYPE2", "CUNIT1", "CUNIT2",
            "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2",
            "CDELT1", "CDELT2", "PC1_1", "PC1_2", "PC2_1", "PC2_2",
            "RADESYS", "EQUINOX")
    ah = fits.Header()
    for key in keys:
        if key in header:
            ah[key] = header[key]
    w = WCS(ah, relax=False)
    if not w.has_celestial or w.pixel_n_dim != 2 or w.world_n_dim != 2:
        raise ValueError("Expected 2D celestial mask WCS")

    xx, yy = np.meshgrid(np.array([0, (nx - 1)/2, nx - 1], float),
                         np.array([0, (ny - 1)/2, ny - 1], float))
    xx, yy = xx.ravel(), yy.ravel()
    ara, adec = pixel_to_sky(xx, yy, header)
    ref_ra, ref_dec = w.all_pix2world(xx, yy, 0)
    sep = angular_distance_arcsec(
        ara, adec, np.deg2rad(ref_ra), np.deg2rad(ref_dec)
    )
    measured = float(np.max(sep))
    if not np.isfinite(measured) or measured > tolerance_arcsec:
        raise ValueError(
            f"Astropy vs analytic TAN WCS disagreement {measured:.6g} arcsec "
            f"(tolerance {tolerance_arcsec})"
        )
    cd = cd_from_header(header)
    cd_error = float(np.max(np.abs(w.pixel_scale_matrix - cd)))
    if cd_error > 1e-12:
        raise ValueError(f"Astropy/CD matrix disagreement: {cd_error}")
    center_area = float(tan_pixel_area_deg2((nx-1)/2, (ny-1)/2, header))
    return {
        "reference": "astropy.wcs.WCS.all_pix2world, origin=0",
        "positions_checked": 9,
        "max_angular_sep_arcsec": measured,
        "tolerance_arcsec": float(tolerance_arcsec),
        "max_cd_matrix_disagreement_deg_per_pix": cd_error,
        "center_solid_angle_deg2_per_native_pixel": center_area,
    }
