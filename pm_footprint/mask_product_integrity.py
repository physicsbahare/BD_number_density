"""Fail-closed provenance gate for scientific JWST star-mask intersections.

The 2026-10-31 raster and its weighted points were made before correction
of FITS PC/CDELT axis scaling. They must NOT enter WHT/epoch estimates.
This module requires a fresh 20-tile mask product and matching sampled points.
It does not claim that 0.3/0.6 arcsec raster areas are exact sky unions.
"""
from __future__ import annotations

import hashlib
import math
from pathlib import Path

OFFICIAL_ARCHIVE_SHA256 = "a06107b1c654387bf3ab5284792201cb299b7942f52dc95e41089b0e2967cd52"
WCS_CONVENTION = "FITS CD[i,j]=CDELT[i]*PC[i,j]; earlier mask union/MC files require WCS validation before scientific use"


def _finite_positive(x, name):
    try:
        value = float(x)
    except (TypeError, ValueError) as e:
        raise ValueError(f"Missing/invalid {name}") from e
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"Nonpositive/nonfinite {name}")
    return value


def validate_mask_product(summary, points_path=None, n_points=None):
    """Validate provenance and basic area bookkeeping; raise on legacy data.

    This is an *integrity gate*, not validation of mask semantics or true
    WHT/epoch intersection. Sample points must match the fresh WCS product.
    """
    if summary.get("wcs_convention") != WCS_CONVENTION:
        raise ValueError("Stale/unverified JWST mask WCS; regenerate from the original 20 FITS tiles with corrected row-scaled PC/CDELT")
    if summary.get("source_archive_sha256") != OFFICIAL_ARCHIVE_SHA256:
        raise ValueError("JWST mask archive SHA256 is missing or differs from frozen official input")
    tiles = summary.get("tiles")
    if not isinstance(tiles, list) or len(tiles) != 20:
        raise ValueError("Expected 20 audited JWST mask tiles")
    nominal_sum = _finite_positive(summary.get("sum_tile_native_masked_area_deg2_NOT_UNION"), "native tile-sum area")
    tile_sum = 0.0
    for t in tiles:
        tile_sum += _finite_positive(t.get("native_area_deg2"), "tile native area")
        _finite_positive(t.get("wcs_center_area_ratio_to_nominal_30mas"), "TAN native-pixel area ratio")
        if len(t.get("wcs_pc_off_diagonal", [])) != 2:
            raise ValueError("Missing measured FITS PC off-diagonal terms")
        for key in ("grid_0.3_area_ratio_to_native", "grid_0.6_area_ratio_to_native"):
            r = _finite_positive(t.get(key), key)
            if not 0.95 <= r <= 1.01:
                raise ValueError(f"Unexpected raster area conservation: {key}={r}")
    if abs(tile_sum - nominal_sum) > max(1e-11, 1e-7 * nominal_sum):
        raise ValueError("Native masked tile sum does not conserve per-tile area")
    grids = summary.get("grid_results", {})
    for g in ("0.3", "0.6"):
        v = grids.get(g, {})
        lower = _finite_positive(v.get("union_lower_bound_deg2"), f"{g} union lower")
        upper = _finite_positive(v.get("union_upper_bound_deg2"), f"{g} union upper")
        if lower > upper + 1e-12 or upper > nominal_sum * 1.01:
            raise ValueError(f"Invalid {g} raster union area ordering/normalization")
    if points_path is not None:
        sample = summary.get("weighted_mc_sample", {})
        if sample.get("source_grid_arcsec") != 0.6:
            raise ValueError("Weighted sample must originate from 0.6 arcsec mask grid")
        if not isinstance(sample.get("n"), int) or sample["n"] < 100:
            raise ValueError("Invalid weighted mask sample count")
        if n_points is not None and sample["n"] != n_points:
            raise ValueError("Weighted mask sample row count mismatch")
        expected = sample.get("sha256")
        if not isinstance(expected, str) or len(expected) != 64:
            raise ValueError("Weighted mask sample missing SHA256: possibly pre-WCS-fix points")
        actual = hashlib.sha256(Path(points_path).read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError("Weighted mask sample SHA256 mismatch; cannot pair with this WCS mask union")
    return True
