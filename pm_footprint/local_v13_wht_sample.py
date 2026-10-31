"""Bounded-memory local WHT>0 spatial audit for Zhen-Kai v13.

Uncompressed 2-D float32 FITS with plain TAN celestial WCS only. Uses
read-only os.pread row spans, not Astropy image loading or memmap.
Results are *sampled geometric diagnostics*, NOT calibrated depth,
dropout-selection effective area, astrometric sensitivity or completeness.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

PHASES = (0.25, 0.75)


def read_header_and_offset(path):
    with open(path, "rb") as stream:
        header = fits.Header.fromfile(stream, endcard=True, padding=True)
        offset = stream.tell()
    if int(header.get("NAXIS", 0)) != 2 or int(header.get("BITPIX", 0)) != -32:
        raise ValueError("Expected 2D primary uncompressed float32 FITS")
    if float(header.get("BSCALE", 1)) != 1 or float(header.get("BZERO", 0)) != 0:
        raise ValueError("Scaled FITS data require a different reader")
    if not (str(header.get("CTYPE1", "")).endswith("TAN")
            and str(header.get("CTYPE2", "")).endswith("TAN")):
        raise ValueError("Plain TAN WCS required; check FITS header first")
    if any(k.startswith(("A_", "B_", "AP_", "BP_")) for k in header.keys()):
        raise ValueError("Distortion keywords present; no plain TAN area approximation")
    nx, ny = int(header["NAXIS1"]), int(header["NAXIS2"])
    if offset + 4 * nx * ny > path.stat().st_size:
        raise ValueError("FITS data are truncated")
    wcs = WCS(header).celestial
    if not wcs.has_celestial:
        raise ValueError("No celestial WCS")
    return header, offset, wcs


def sample_raw(fd, offset, nx, row, columns):
    """FITS primary image is row-major, big-endian float32."""
    columns = np.asarray(columns, dtype=np.int64)
    if not columns.size:
        return np.empty(0, dtype=np.float32)
    start, end = int(columns.min()), int(columns.max())
    raw = os.pread(fd, (end - start + 1) * 4,
                   int(offset + (int(row) * nx + start) * 4))
    if len(raw) != (end - start + 1) * 4:
        raise IOError("Incomplete FITS row read")
    return np.frombuffer(raw, dtype=">f4")[columns - start]


def tan_cell_areas(x, y, widths, height, matrix, crpix):
    """Solid-angle corrected pixel areas in deg^2 for undistorted TAN WCS."""
    dx, dy = np.asarray(x, dtype=float) + 1 - crpix[0], float(y) + 1 - crpix[1]
    xi = np.deg2rad(matrix[0, 0] * dx + matrix[0, 1] * dy)
    eta = np.deg2rad(matrix[1, 0] * dx + matrix[1, 1] * dy)
    return (abs(np.linalg.det(matrix)) * widths * height
            / np.power(1 + xi * xi + eta * eta, 1.5))


def sample_phase(fd, offset, wcs, nx, ny, stride, phase):
    left = np.arange(0, nx, stride, dtype=np.int64)
    top = np.arange(0, ny, stride, dtype=np.int64)
    widths = np.minimum(stride, nx - left)
    heights = np.minimum(stride, ny - top)
    xs = left + np.floor(phase * widths).astype(np.int64)
    ys = top + np.floor(phase * heights).astype(np.int64)
    occupancy = np.empty((len(ys), len(xs)), dtype=np.bool_)
    cd = np.asarray(wcs.pixel_scale_matrix, dtype=float)
    crpix = np.asarray(wcs.wcs.crpix, dtype=float)
    total = positive = finite = 0.0
    for j, (row, height) in enumerate(zip(ys, heights)):
        vals = sample_raw(fd, offset, nx, row, xs)
        isfinite = np.isfinite(vals)
        ispositive = isfinite & (vals > 0)
        occupancy[j] = ispositive
        area = tan_cell_areas(xs, row, widths, height, cd, crpix)
        total += float(area.sum())
        finite += float(area[isfinite].sum())
        positive += float(area[ispositive].sum())
        if (j + 1) % 200 == 0:
            print(f"phase {phase}: {j+1}/{len(ys)} rows", flush=True)
    return {
        "phase": phase, "grid_shape": list(occupancy.shape),
        "sampled_cells": int(occupancy.size),
        "positive_cells": int(occupancy.sum()),
        "rectangle_area_deg2": total,
        "finite_sampled_area_deg2": finite,
        "positive_wht_sampled_area_deg2": positive
    }, occupancy


def check_candidates(fd, offset, wcs, nx, ny, candidates):
    data = []
    for c in candidates:
        x, y = wcs.world_to_pixel_values(float(c["ra_deg"]), float(c["dec_deg"]))
        ix, iy = int(np.rint(x)), int(np.rint(y))
        inside = 0 <= ix < nx and 0 <= iy < ny
        value = float(sample_raw(fd, offset, nx, iy, [ix])[0]) if inside else None
        if value is not None and not np.isfinite(value):
            value = None
        data.append({
            "candidate_id": c["candidate_id"], "ra_deg": c["ra_deg"],
            "dec_deg": c["dec_deg"], "pixel_x": float(x), "pixel_y": float(y),
            "inside_fits_rectangle": bool(inside), "wht": value,
            "wht_positive": bool(value is not None and value > 0)
        })
    return data


def run(v13, candidates_csv, stride, prefix):
    if stride < 2:
        raise ValueError("stride must be >=2")
    header, offset, wcs = read_header_and_offset(v13)
    nx, ny = int(header["NAXIS1"]), int(header["NAXIS2"])
    with open(candidates_csv, newline="") as f:
        candidates = list(csv.DictReader(f))
    if not candidates or len(set(c["candidate_id"] for c in candidates)) != len(candidates):
        raise ValueError("Candidate IDs must be nonempty and unique")
    roster_sha = hashlib.sha256(candidates_csv.read_bytes()).hexdigest()
    fd = os.open(v13, os.O_RDONLY)
    try:
        stats, grids = [], {}
        for phase in PHASES:
            s, grid = sample_phase(fd, offset, wcs, nx, ny, stride, phase)
            stats.append(s)
            grids[f"wht_positive_{int(phase*100)}"] = grid
        points = check_candidates(fd, offset, wcs, nx, ny, candidates)
    finally:
        os.close(fd)
    areas = [s["positive_wht_sampled_area_deg2"] for s in stats]
    result = {
        "status": "PROVISIONAL_SAMPLED_WHT_GEOMETRY_NOT_SELECTION_AREA",
        "source_basename": v13.name, "source_file_bytes": v13.stat().st_size,
        "source_header_sha256": hashlib.sha256(header.tostring(sep="\n").encode()).hexdigest(),
        "roster_sha256": roster_sha, "stride_pixels": stride,
        "native_shape_yx": [ny, nx], "phase_stats": stats,
        "sampled_mean_area_deg2": float(np.mean(areas)),
        "phase_difference_area_deg2": abs(areas[0] - areas[1]),
        "candidate_count": len(points),
        "candidate_centers_inside_image": sum(p["inside_fits_rectangle"] for p in points),
        "candidate_centers_wht_positive": sum(p["wht_positive"] for p in points),
        "caveats": [
            "Only a sampled WHT>0 map, not an exact area or formal error interval.",
            "Subgrid star masks, edges and holes may be missed.",
            "WHT positive is not calibrated F444W=26.2 selection completeness.",
            "Star-mask and independent-epoch intersections are not done here.",
            "No cross-catalog WHT amplitude comparison without calibration."
        ]
    }
    prefix.parent.mkdir(parents=True, exist_ok=True)
    prefix.with_suffix(".json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    with open(prefix.parent / (prefix.name + "_candidates.csv"), "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(points[0]))
        writer.writeheader()
        writer.writerows(points)
    np.savez_compressed(prefix.parent / (prefix.name + "_grid.npz"),
                        stride=np.array(stride), ncols=np.array(nx), nrows=np.array(ny),
                        wcs_cd=wcs.pixel_scale_matrix, wcs_crpix=wcs.wcs.crpix,
                        wcs_crval=wcs.wcs.crval, **grids)
    print("Summary:", json.dumps({
        "sampled_mean_area_deg2": result["sampled_mean_area_deg2"],
        "phase_difference_area_deg2": result["phase_difference_area_deg2"],
        "candidate_centers_wht_positive": result["candidate_centers_wht_positive"]
    }))
    print("Saved JSON, candidate CSV and compressed grid under:", prefix)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--v13", type=Path, required=True)
    p.add_argument("--candidates", type=Path,
                   default=Path("pm_footprint/dave125_positions.csv"))
    p.add_argument("--stride", type=int, default=120)
    p.add_argument("--out-prefix", type=Path,
                   default=Path("pm_footprint/results/local_v13_sample"))
    a = p.parse_args()
    run(a.v13, a.candidates, a.stride, a.out_prefix)


if __name__ == "__main__":
    main()
