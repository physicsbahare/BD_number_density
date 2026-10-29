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


def official_corner_wcs(tile, nx, ny):
    """Reconstruct a TAN transform from the official DR1 tile-corner table.

    The extracted standalone WHT FITS has no WCS cards (confirmed for A1).
    This transform is *derived*, not the original SCI header: require
    independent full-i2d SCI WCS comparison before a final paper claim.
    """
    corners_file = Path(__file__).with_name("dr1_tile_corners.csv")
    table = pd.read_csv(corners_file).set_index("tile")
    if tile not in table.index:
        raise ValueError(f"No official DR1 tile corners for {tile}")
    row = table.loc[tile]
    ra = np.array([row[f"ra{i}"] for i in range(1, 5)], dtype=float)
    dec = np.array([row[f"dec{i}"] for i in range(1, 5)], dtype=float)
    ra0, dec0 = np.deg2rad([150.1163213, 2.200973097])
    lon, lat = np.deg2rad(ra), np.deg2rad(dec)
    delta = lon - ra0
    denominator = (np.sin(dec0) * np.sin(lat)
                   + np.cos(dec0) * np.cos(lat) * np.cos(delta))
    if np.any(denominator <= 0):
        raise ValueError(f"{tile}: tile corners outside TAN hemisphere")
    xi = np.cos(lat) * np.sin(delta) / denominator
    eta = (np.cos(dec0) * np.sin(lat)
           - np.sin(dec0) * np.cos(lat) * np.cos(delta)) / denominator
    plane = np.rad2deg(np.column_stack((xi, eta)))
    # Published corner order is bottom-left, bottom-right, top-right, top-left
    # in FITS array pixel (x,y) coordinates; infer both axes from actual shape.
    pixels = np.array([[0., 0.], [nx - 1., 0.],
                       [nx - 1., ny - 1.], [0., ny - 1.]])
    design = np.column_stack((np.ones(4), pixels))
    coeff, _, _, _ = np.linalg.lstsq(design, plane, rcond=None)
    cd = coeff[1:].T
    scale = np.linalg.norm(cd, axis=0) * 3600
    if not np.all((scale > 0.055) & (scale < 0.065)):
        raise ValueError(f"{tile}: inferred WCS pixel scales {scale} not ~60 mas")
    if np.linalg.det(cd) >= 0:
        raise ValueError(f"{tile}: unexpected sky orientation / WCS parity")
    fit_resid = np.linalg.norm(design @ coeff - plane, axis=1) * 3600
    if float(fit_resid.max()) > 1.0:
        raise ValueError(f"{tile}: official-corner TAN fit residual {fit_resid} arcsec")
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.crval = [150.1163213, 2.200973097]
    wcs.wcs.cd = cd
    wcs.wcs.crpix = np.ones(2) - np.linalg.solve(cd, coeff[0])
    ra_fit, dec_fit = wcs.all_pix2world(pixels[:, 0], pixels[:, 1], 0)
    residual = np.hypot((ra_fit - ra) * np.cos(np.deg2rad(dec)),
                        dec_fit - dec) * 3600
    if float(residual.max()) > 1.0:
        raise ValueError(f"{tile}: WCS does not reproduce official corners: {residual}")
    return wcs, {
        "wcs_source": "reconstructed_from_official_DR1_tile_corners_NOT_original_header",
        "wcs_corner_residual_arcsec_max": float(residual.max()),
        "wcs_pixel_scales_arcsec": [float(v) for v in scale],
        "wcs_reference": "https://cosmos2025.iap.fr/nircam.html",
        "wcs_independent_i2d_validation": "PENDING",
    }


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
            ny, nx = data.shape
            try:
                wcs, wcs_source = select_celestial_wcs(hdul, hdu, tile)
                wcs_meta = {"wcs_source": wcs_source}
            except ValueError as err:
                # Standalone DR1 WHT has no celestial WCS cards. The
                # authoritative original i2d SCI header *does* have WCS.
                # Read only 4 MiB of its gzip prefix, never the 2-GB image.
                from .probe_i2d_wcs import find_i2d_url, i2d_sci_header
                print(f"{tile}: {err}; reading original i2d SCI WCS prefix",
                      flush=True)
                i2d_url = find_i2d_url(tile)
                sci, wcs, prefix_meta = i2d_sci_header(i2d_url)
                if (int(sci.get("NAXIS1", -1)) != nx
                        or int(sci.get("NAXIS2", -1)) != ny):
                    raise RuntimeError(
                        f"{tile}: original SCI and extracted WHT dimensions differ")
                corner_wcs, corner_meta = official_corner_wcs(tile, nx, ny)
                pts = np.array([[0., 0.], [nx - 1., 0.],
                                [nx - 1., ny - 1.], [0., ny - 1.],
                                [(nx - 1.) / 2, (ny - 1.) / 2]])
                ra0, dec0 = wcs.all_pix2world(pts[:, 0], pts[:, 1], 0)
                ra1, dec1 = corner_wcs.all_pix2world(pts[:, 0], pts[:, 1], 0)
                diff_arcsec = np.hypot(
                    (ra0 - ra1) * np.cos(np.deg2rad(dec0)),
                    dec0 - dec1) * 3600
                if not np.all(np.isfinite(diff_arcsec)) or diff_arcsec.max() > 2.0:
                    raise RuntimeError(
                        f"{tile}: original SCI WCS disagrees with published "
                        f"tile corners by {diff_arcsec} arcsec")
                wcs_meta = {
                    "wcs_source": "ORIGINAL_DR1_i2d_SCI_header_prefix",
                    "wcs_i2d_url": i2d_url,
                    "wcs_vs_published_corners_arcsec_max":
                        float(diff_arcsec.max()),
                    "wcs_corner_fit": corner_meta,
                    "wcs_prefix_http_status": prefix_meta["http_status"],
                    "wcs_prefix_bytes": prefix_meta["compressed_prefix_bytes"],
                }
            by, bx = ny // block, nx // block
            sample = np.asarray(data[::48, ::48], dtype=float)
            positive = sample[np.isfinite(sample) & (sample > 0)]
            if not len(positive):
                raise ValueError(f"No positive WHT pixels in {tile}")
            median = float(np.median(positive))
            n = grid.shape[1]
            # A rotated native 0.6-arcsec block cannot be assigned to a
            # single 0.6-arcsec sky cell by nearest-neighbour sampling:
            # that creates artificial holes and underestimates area.
            # Accumulate each block's fractional area bilinearly into its
            # four surrounding equal-area grid centres, *within this tile*.
            # The final cross-tile maximum then approximates the union.
            tile_grid = np.zeros(grid.shape, dtype=np.float32)
            native_fraction_sums = np.zeros(len(THRESHOLDS), dtype=float)
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
                xfloat = (xx_sky - lo) / step - 0.5
                yfloat = (yy_sky - lo) / step - 0.5
                ix0 = np.floor(xfloat).astype(np.int64)
                iy0 = np.floor(yfloat).astype(np.int64)
                fx, fy = xfloat - ix0, yfloat - iy0
                active = block_fractions[0] > 0
                used += int(np.count_nonzero(active))
                for k, frac in enumerate(block_fractions):
                    native_fraction_sums[k] += float(frac.sum())
                # Detect positive-WHT blocks outside the *entire* grid.
                if np.any(active & ((ix0 + 1 < 0) | (ix0 >= n)
                                    | (iy0 + 1 < 0) | (iy0 >= n))):
                    lost += int(np.count_nonzero(
                        active & ((ix0 + 1 < 0) | (ix0 >= n)
                                  | (iy0 + 1 < 0) | (iy0 >= n))))
                for dx in (0, 1):
                    for dy in (0, 1):
                        ix, iy = ix0 + dx, iy0 + dy
                        good = (active & (ix >= 0) & (ix < n)
                                & (iy >= 0) & (iy < n))
                        if not np.any(good):
                            continue
                        wx = fx if dx else 1 - fx
                        wy = fy if dy else 1 - fy
                        weights = (wx * wy)[good]
                        flat_idx = (iy[good] * n + ix[good]).ravel()
                        for k, fractions in enumerate(block_fractions):
                            np.add.at(tile_grid[k].ravel(), flat_idx,
                                      (fractions[good] * weights).astype(np.float32))
            if lost:
                raise RuntimeError(f"{tile}: {lost} valid WHT coarse blocks outside grid")
            cell_area = step * step * DEG2_PER_SR
            native_area = native_fraction_sums * cell_area
            # Clipping at 1 represents fractional occupancy of a sky cell;
            # area should remain close to the native positive-pixel area.
            np.minimum(tile_grid, 1.0, out=tile_grid)
            gridded_area = tile_grid.sum(axis=(1, 2), dtype=np.float64) * cell_area
            area_error = np.divide(gridded_area - native_area, native_area,
                                   out=np.zeros_like(native_area),
                                   where=native_area > 0)
            if np.any(np.abs(area_error) > 0.03):
                raise RuntimeError(
                    f"{tile}: fractional-grid area conservation failed "
                    f"(native={native_area}, gridded={gridded_area}, "
                    f"fractional_error={area_error}); do not quote area")
            np.maximum(grid, np.rint(tile_grid * 100).astype(np.uint8), out=grid)
            return {
                "tile": tile, "url": url, "gz_sha256": digest,
                **wcs_meta,
                "median_positive_wht_sample": median,
                "fits_shape": [ny, nx], "coarse_cells_mapped": used,
                "native_tile_area_deg2": [float(v) for v in native_area],
                "gridded_tile_area_deg2": [float(v) for v in gridded_area],
                "grid_area_conservation_error_fraction": [float(v) for v in area_error],
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


def grid_resolution_qc(grid, lo, step, reference, roster, survey, masks):
    """Compare WHT-polygon integration under coarsening and grid translation.

    This tests *integration* sensitivity without downloading the 7GB FITS
    again; it is not an independent 60mas resampling of the native images.
    """
    n = grid.shape[1]
    if n % 2:
        raise ValueError("Resolution-QC grid requires even dimensions")
    coarse = np.rint(grid.reshape(len(THRESHOLDS), n // 2, 2,
                                  n // 2, 2).mean(axis=(2, 4))).astype(np.uint8)
    coarse_result, _ = evaluate(coarse, lo, step * 2, {}, roster, survey, masks)
    del coarse
    # Half-cell shift: interpolate occupancy rather than shifting the sky.
    shifted = np.rint(
        (grid[:, :-1, :-1].astype(np.float32)
         + grid[:, 1:, :-1].astype(np.float32)
         + grid[:, :-1, 1:].astype(np.float32)
         + grid[:, 1:, 1:].astype(np.float32)) / 4
    ).astype(np.uint8)
    shifted_result, _ = evaluate(shifted, lo + step / 2, step,
                                 {}, roster, survey, masks)
    del shifted
    output = {"method": "post_rasterization_occupancy_grid_sensitivity",
              "limitations": [
                  "1.2arcsec grid is a 2x2 mean of the 0.6arcsec grid, not an independent 60mas FITS resample.",
                  "Half-cell translation uses bilinear occupancy interpolation.",
                  "Grid convergence does not establish photometric depth or astrometric completeness."
              ],
              "comparisons": {}}
    max_delta = 0.0
    max_count_delta = 0
    for label in reference["area_deg2"]:
        output["comparisons"][label] = {}
        for name, base in reference["area_deg2"][label].items():
            row = {}
            for key, alt in (("coarsened_1p2arcsec", coarse_result),
                             ("shifted_0p3arcsec", shifted_result)):
                new = alt["area_deg2"][label][name]
                delta = float((new - base) / base) if base > 0 else 0.0
                count_delta = (alt["candidate_counts"][label][name]
                               - reference["candidate_counts"][label][name])
                row[key] = {
                    "area_deg2": float(new),
                    "relative_area_delta": delta,
                    "candidate_count_delta": int(count_delta),
                }
                max_delta = max(max_delta, abs(delta))
                max_count_delta = max(max_count_delta, abs(count_delta))
            output["comparisons"][label][name] = row
    output["max_abs_relative_area_delta"] = max_delta
    output["max_abs_candidate_count_delta"] = max_count_delta
    output["status"] = ("GRID_INTEGRATION_STABLE_WITHIN_2PCT_AND_2_CANDIDATES"
                        if max_delta <= 0.02 and max_count_delta <= 2
                        else "GRID_INTEGRATION_REVIEW_REQUIRED")
    return output


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
        grid_qc = grid_resolution_qc(grid, lo, step, result, roster, survey, masks)
        (out / "grid_resolution_qc.json").write_text(
            json.dumps(grid_qc, indent=2) + "\n")
        result.update({
            "grid_resolution_qc_status": grid_qc["status"],
            "grid_resolution_qc_max_relative_area_delta":
                grid_qc["max_abs_relative_area_delta"],
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
