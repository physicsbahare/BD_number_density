"""Intersect exact JWST multi-epoch polygons with COSMOS-Web DR1 F444W WHT.

Streams the 20 official *60mas* WHT tiles sequentially. Produces a
coarse, equal-area, fractionally occupied grid (not a claimed 60mas mask).
The three WHT thresholds are sensitivity experiments, not a validated
F444W 26.2-mag selection completeness map.

Requires: astropy, numpy, pandas, shapely>=2, scipy, curl, disk >=3 GB.
No local workstation data and no duplicate MAST query.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import tempfile
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from shapely import contains_xy, prepare

from .exact_geometry import compute, lambert, DEG2_PER_SR

BASE = "https://cosmos2025.iap.fr/data/nircam/extensions"
TILES = [f"{row}{j}" for row in "AB" for j in range(1, 11)]
THRESHOLDS = (0.0, 0.1, 0.5)
GRID_HALF_DEG = 0.62


def url_for(tile):
    if tile not in TILES:
        raise ValueError(f"Unknown COSMOS-Web tile {tile}")
    return (f"{BASE}/mosaic_nircam_f444w_COSMOS-Web_"
            f"60mas_{tile}_v1.0_wht.fits.gz")


def download(url, path):
    subprocess.run([
        "curl", "--fail", "--location", "--retry", "4",
        "--retry-all-errors", "--connect-timeout", "40",
        "--max-time", "2400", "--silent", "--show-error",
        "--output", str(path), url
    ], check=True)
    if path.stat().st_size < 1024 * 1024:
        raise RuntimeError(f"Suspiciously small FITS tile: {url}")


def select_celestial_wcs(hdul, image_hdu, tile):
    """Find WCS in the image, primary, or a combined header.

    Official extracted FITS extensions need not copy the SCI WCS into
    their own image header. Fail closed with header diagnostics if neither
    header provides a usable celestial transform.
    """
    headers = [("image", image_hdu.header), ("primary", hdul[0].header)]
    merged = hdul[0].header.copy()
    merged.extend(image_hdu.header, update=True, strip=False)
    headers.append(("primary_plus_image", merged))
    for label, header in headers:
        try:
            candidate = WCS(header, relax=True)
            if candidate.has_celestial and candidate.pixel_n_dim >= 2:
                return candidate.celestial, label
        except Exception:
            continue
    diagnostic = []
    for i, hdu in enumerate(hdul):
        header = hdu.header
        diagnostic.append({
            "hdu": i, "name": hdu.name,
            "naxis": header.get("NAXIS"),
            "naxis1": header.get("NAXIS1"),
            "naxis2": header.get("NAXIS2"),
            "wcs_keys": {k: str(header[k]) for k in
                         ("WCSAXES", "CTYPE1", "CTYPE2", "CRVAL1", "CRVAL2",
                          "CRPIX1", "CRPIX2", "CD1_1", "CDELT1", "PC1_1")
                         if k in header},
        })
    raise ValueError(
        f"No celestial WCS in extracted WHT {tile}; "
        f"HDU header diagnostics: {json.dumps(diagnostic)}. "
        "Do not infer coordinates from pixel indices."
    )


def wht_tile(tile, grid, step, lo, block, work):
    gz = work / f"{tile}.fits.gz"
    fits_path = work / f"{tile}.fits"
    url = url_for(tile)
    try:
        download(url, gz)
        digest = hashlib.sha256(gz.read_bytes()).hexdigest()
        with gzip.open(gz, "rb") as fin, open(fits_path, "wb") as fout:
            shutil.copyfileobj(fin, fout, length=16 * 1024 * 1024)
        gz.unlink()
        with fits.open(fits_path, memmap=True) as hdul:
            hdu = next((h for h in hdul if h.data is not None
                        and getattr(h.data, "ndim", None) == 2), None)
            if hdu is None:
                raise ValueError(f"No 2D image in WHT tile {tile}")
            data = hdu.data
            wcs, wcs_source = select_celestial_wcs(hdul, hdu, tile)
            ny, nx = data.shape
            by, bx = ny // block, nx // block
            sample = np.asarray(data[::48, ::48], dtype=float)
            positive = sample[np.isfinite(sample) & (sample > 0)]
            if not len(positive):
                raise ValueError(f"No positive WHT pixels in {tile}")
            median = float(np.median(positive))
            n = grid.shape[1]
            lost = 0
            used = 0
            for by0 in range(0, by, 64):
                by1 = min(by, by0 + 64)
                raw = np.asarray(data[by0 * block:by1 * block, :bx * block])
                cube = raw.reshape(by1 - by0, block, bx, block)
                finite = np.isfinite(cube)
                block_fractions = [
                    np.mean(finite & (cube > (median * t)), axis=(1, 3))
                    for t in THRESHOLDS
                ]
                ypix = (np.arange(by0, by1) * block + (block - 1) / 2)
                xpix = (np.arange(bx) * block + (block - 1) / 2)
                xx, yy = np.meshgrid(xpix, ypix)
                ra, dec = wcs.all_pix2world(xx, yy, 0)
                xx_sky, yy_sky = lambert(ra, dec)
                ix = np.floor((xx_sky - lo) / step).astype(np.int64)
                iy = np.floor((yy_sky - lo) / step).astype(np.int64)
                in_grid = ((ix >= 0) & (ix < n)
                           & (iy >= 0) & (iy < n))
                if np.any((block_fractions[0] > 0) & ~in_grid):
                    lost += int(np.count_nonzero((block_fractions[0] > 0) & ~in_grid))
                good = in_grid & (block_fractions[0] > 0)
                if not np.any(good):
                    continue
                flat_idx = (iy[good] * n + ix[good]).ravel()
                used += int(len(flat_idx))
                for k, fractions in enumerate(block_fractions):
                    occupancy = np.rint(100 * fractions[good]).astype(np.uint8)
                    np.maximum.at(grid[k].ravel(), flat_idx, occupancy)
            if lost:
                raise RuntimeError(f"{tile}: {lost} valid WHT coarse cells outside grid")
            return {
                "tile": tile, "url": url, "gz_sha256": digest,
                "wcs_source": wcs_source,
                "median_positive_wht_sample": median,
                "fits_shape": [ny, nx], "coarse_cells_mapped": used,
                "trimmed_native_edge_pixels": [ny % block, nx % block],
            }
    finally:
        for f in (gz, fits_path):
            f.unlink(missing_ok=True)


def exact_area_limit(survey):
    # Rasterization may differ from exact polygon area by a small edge term.
    return survey.area * DEG2_PER_SR + 0.003


def evaluate(grid, lo, step, summary, roster, survey, masks):
    # A cell contains a 0..100 occupancy fraction from the WHT 60mas tile.
    # Fractional pixel coverage reduces boundary bias vs any-positive binning.
    cell_area = step * step * DEG2_PER_SR
    n = grid.shape[1]
    geoms = {"archive_f444w": survey}
    for gap, d in masks.items():
        for key, shape in d.items():
            geoms[f"{key}_{gap}d"] = shape
    for shape in geoms.values():
        prepare(shape)
    totals = {name: np.zeros(len(THRESHOLDS), dtype=float) for name in geoms}
    # Only process coarse grid cells with at least some WHT data.
    for y0 in range(0, n, 96):
        sl = grid[:, y0:min(y0 + 96, n), :]
        yy, xx = np.nonzero(sl[0] > 0)
        if len(xx) == 0:
            continue
        xxsky = lo + (xx.astype(float) + 0.5) * step
        yysky = lo + (yy.astype(float) + y0 + 0.5) * step
        weights = sl[:, yy, xx].astype(float) / 100.0
        for name, shape in geoms.items():
            inside = contains_xy(shape, xxsky, yysky)
            totals[name] += weights[:, inside].sum(axis=1) * cell_area

    px, py = lambert(roster.ra_deg.to_numpy(), roster.dec_deg.to_numpy())
    ix = np.floor((px - lo) / step).astype(int)
    iy = np.floor((py - lo) / step).astype(int)
    if np.any((ix < 0) | (iy < 0) | (ix >= n) | (iy >= n)):
        raise RuntimeError("Dave125 coordinates outside WHT sampling grid")
    membership = roster[["candidate_id", "ra_deg", "dec_deg"]].copy()
    for k, t in enumerate(THRESHOLDS):
        membership[f"wht_{t:g}_fraction"] = grid[k, iy, ix].astype(float) / 100
    output = {"method": "fractional_60mas_WHT_binned_on_equal_area_grid",
              "wht_thresholds_relative_to_each_tile_positive_median": list(THRESHOLDS),
              "grid_cell_arcsec": math.sqrt(cell_area) * 3600,
              "area_deg2": {}, "fractions": {}, "candidate_counts": {}}
    for k, t in enumerate(THRESHOLDS):
        label = f"{t:g}"
        base = totals["archive_f444w"][k]
        minimum = 0.3 if t == 0 else 0.01
        if not (minimum < base <= exact_area_limit(survey)):
            raise RuntimeError(
                f"Implausible WHT/archival intersection area {base:.6f} deg2 "
                f"for threshold {t}; inspect WCS, tile coverage and mask")
        output["area_deg2"][label] = {name: float(val[k]) for name, val in totals.items()}
        output["fractions"][label] = {name: float(val[k] / base) for name, val in totals.items()}
        wht_good = grid[k, iy, ix] > 0
        output["candidate_counts"][label] = {}
        for name, geom in geoms.items():
            inside = contains_xy(geom, px, py)
            membership[f"in_{name}"] = inside
            output["candidate_counts"][label][name] = int(np.sum(wht_good & inside))
    return output, membership


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--metadata", default="pm_footprint/results/mast_observations.csv")
    p.add_argument("--candidates", default="pm_footprint/dave125_positions.csv")
    p.add_argument("--out", default="pm_footprint/results/wht60")
    p.add_argument("--block", type=int, default=10,
                   help="Native 60mas pixels per sampled block; default 10=0.6 arcsec")
    p.add_argument("--tiles", nargs="*", default=TILES)
    a = p.parse_args()
    if not 2 <= a.block <= 40:
        p.error("--block must be 2..40")
    if sorted(a.tiles) != sorted(TILES):
        p.error("All 20 DR1 tiles required for a complete-area claim")
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    meta_path, roster_path = Path(a.metadata), Path(a.candidates)
    meta = pd.read_csv(meta_path)
    roster = pd.read_csv(roster_path)
    # This reuses the saved MAST snapshot; no expensive archive query.
    exact, _, survey, masks = compute(meta, roster, return_masks=True)
    step = math.radians(a.block * 0.06 / 3600)
    lo = math.radians(-GRID_HALF_DEG)
    n = int(math.ceil(2 * GRID_HALF_DEG * 3600 / (a.block * 0.06)))
    grid = np.zeros((len(THRESHOLDS), n, n), dtype=np.uint8)
    tiles = []
    try:
        with tempfile.TemporaryDirectory(prefix="cosmosweb_wht_") as temp:
            work = Path(temp)
            for tile in TILES:
                record = wht_tile(tile, grid, step, lo, a.block, work)
                tiles.append(record)
                print(f"WHT {tile}: {record['fits_shape']} median={record['median_positive_wht_sample']:.4g}", flush=True)
                (out / "tile_progress.json").write_text(json.dumps(tiles, indent=2) + "\n")
        result, membership = evaluate(grid, lo, step, exact, roster, survey, masks)
        result.update({
            "status": "WHT_GEOMETRIC_INTERSECTION_PROVISIONAL_REQUIRES_DEPTH_AND_RESOLUTION_QC",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "official_release": "https://cosmos2025.iap.fr/nircam.html",
            "exact_archive_area_deg2": exact["archive_f444w_area_deg2"],
            "wht_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "exact_source_sha256": hashlib.sha256(Path(__import__(
                "pm_footprint.exact_geometry", fromlist=["__file__"]).__file__).read_bytes()).hexdigest(),
            "mast_metadata_sha256": hashlib.sha256(meta_path.read_bytes()).hexdigest(),
            "candidate_catalog_sha256": hashlib.sha256(roster_path.read_bytes()).hexdigest(),
            "tiles": tiles,
            "caveats": [
                "F444W DR1 60mas WHT>0 is an image-coverage mask, not a selection-matched F444W 26.2 AB depth mask.",
                "Relative per-tile WHT median thresholds 0, 0.1, 0.5 are sensitivity tests, not photometric completeness limits.",
                "WHT blocks are binned to an equal-area grid, using maximum fractional occupancy in overlapping tiles.",
                "Spatial-grid resolution and offsets require convergence tests before publication.",
                "F444W+red cross-filter coverage is not proof of red-source detection or astrometric precision.",
                "Close candidate pairs are not deblended; counts are catalogue entries, not guaranteed independent objects.",
                "Only publicly released JWST metadata used; no exclusive-access data.",
            ],
        })
        (out / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
        membership.to_csv(out / "dave125_wht_membership.csv", index=False)
        print(json.dumps({"status": result["status"], "area_deg2": result["area_deg2"],
                          "fractions": result["fractions"], "counts": result["candidate_counts"]}, indent=2))
    except Exception as exc:
        (out / "failure.json").write_text(json.dumps({
            "status": "FAILED_NO_USABLE_AREA_CLAIM",
            "error": repr(exc),
            "tiles_completed": len(tiles),
            "utc": datetime.now(timezone.utc).isoformat(),
        }, indent=2) + "\n")
        raise


if __name__ == "__main__":
    main()
