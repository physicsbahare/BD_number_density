"""Overlay local sampled v13 WHT>0 cells with frozen exact >=365-day JWST epochs.

No raw 51GB WHT reads, no MAST queries, no 7GB redownload. Requires the
already-committed small MAST snapshot and uploaded compact v13 sample.
Result is GEOMETRIC ONLY, not star-mask or dropout completeness.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.wcs import WCS
from shapely import contains_xy
from .exact_geometry import compute, lambert
from .local_v13_wht_sample import tan_cell_areas


def cell_centers(stride, nx, ny, phase):
    xx, yy = np.arange(0, nx, stride), np.arange(0, ny, stride)
    wx, wy = np.minimum(stride, nx-xx), np.minimum(stride, ny-yy)
    return (xx+np.floor(phase*wx).astype(int),
            yy+np.floor(phase*wy).astype(int), wx, wy)


def reconstruct_wcs(grid):
    w = WCS(naxis=2)
    w.wcs.cd = grid["wcs_cd"]
    w.wcs.crpix = grid["wcs_crpix"]
    w.wcs.crval = grid["wcs_crval"]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    return w


def evaluate_phase(grid, phase, polygons, chunk_size=25000):
    stride = int(grid["stride"])
    nx, ny = int(grid["ncols"]), int(grid["nrows"])
    wht = grid["wht_positive_"+str(int(100*phase))]
    xs, ys, widths, heights = cell_centers(stride, nx, ny, phase)
    if wht.shape != (len(ys), len(xs)):
        raise ValueError("Sample grid shape differs from header")
    cells = np.flatnonzero(wht.ravel())
    wcs = reconstruct_wcs(grid)
    sums = {name: 0.0 for name in polygons}
    total, outside_count = 0.0, 0
    for start in range(0, len(cells), chunk_size):
        flat = cells[start:start+chunk_size]
        iy, ix = np.divmod(flat, len(xs))
        x, y = xs[ix], ys[iy]
        ra, dec = wcs.all_pix2world(x, y, 0)
        px, py = lambert(ra, dec)
        area = tan_cell_areas(x, y, widths[ix], heights[iy],
                              grid["wcs_cd"], grid["wcs_crpix"])
        total += float(area.sum())
        in_archive = contains_xy(polygons["archive"], px, py)
        outside_count += int((~in_archive).sum())
        for key, geom in polygons.items():
            selected = in_archive if key == "archive" else contains_xy(geom, px, py)
            sums[key] += float(area[selected].sum())
        if (start//chunk_size) % 8 == 0:
            print("Phase",phase,"processed",min(start+chunk_size,len(cells)),
                  "/",len(cells),"positive cells",flush=True)
    return {"phase":phase,"wht_positive_cells":int(len(cells)),
            "wht_positive_area_deg2":total,
            "wht_positive_centers_outside_archive":outside_count,
            "area_intersections_deg2":sums,
            "red_epoch_fraction_of_v13_wht":sums["red_365d"]/total}


def run(grid_path, sample_path, metadata_path, roster_path, saved_exact_path):
    def sha(path):
        return hashlib.sha256(Path(path).read_bytes()).hexdigest()
    saved = json.loads(Path(saved_exact_path).read_text())
    sample = json.loads(Path(sample_path).read_text())
    if sha(metadata_path) != saved["mast_metadata_sha256"]:
        raise ValueError("Frozen MAST snapshot SHA changed")
    if sha(roster_path) != saved["candidate_catalog_sha256"] or sha(roster_path) != sample["roster_sha256"]:
        raise ValueError("Frozen candidate roster SHA changed")
    print("Building exact JWST polygons from committed CSV (no network calls)",flush=True)
    computed, membership, survey, masks = compute(
        pd.read_csv(metadata_path), pd.read_csv(roster_path),
        baselines=(365,), return_masks=True)
    if abs(computed["archive_f444w_area_deg2"]-saved["archive_f444w_area_deg2"]) > 1e-8:
        raise ValueError("Archive geometry area mismatch")
    red_area = computed["baseline_results"]["365"]["area_deg2"]["f444w_plus_red"]
    ref_red = saved["baseline_results"]["365"]["area_deg2"]["f444w_plus_red"]
    if abs(red_area-ref_red) > 1e-8:
        raise ValueError("Red-epoch geometry area mismatch")
    if int(membership["in_f444w_plus_red_365d"].sum()) != 62:
        raise ValueError("Frozen red-epoch candidate membership mismatch")
    polygons = {"archive":survey,"red_365d":masks["365"]["f444w_plus_red"],
                "same_f444w_365d":masks["365"]["same_f444w"],
                "any_nircam_365d":masks["365"]["any_nircam"]}
    with np.load(grid_path, allow_pickle=False) as g:
        results = [evaluate_phase(g, phase, polygons) for phase in (.25,.75)]
    for actual, original in zip(results, sample["phase_stats"]):
        if actual["phase"] != original["phase"]:
            raise ValueError("Phase mismatch")
        if abs(actual["wht_positive_area_deg2"]-original["positive_wht_sampled_area_deg2"]) > 1e-7:
            raise ValueError("Recomputed TAN pixel area differs from original upload")
    overlaps = [r["area_intersections_deg2"]["red_365d"] for r in results]
    fractions = [r["red_epoch_fraction_of_v13_wht"] for r in results]
    return {"status":"PROVISIONAL_SAMPLED_V13_WHT_EXACT_ARCHIVE_EPOCH_OVERLAY",
            "input_sha256": {"grid":sha(grid_path),"sample":sha(sample_path),
                             "metadata":sha(metadata_path),"roster":sha(roster_path)},
            "exact_archive_deg2":computed["archive_f444w_area_deg2"],
            "exact_red_epoch_deg2":red_area,
            "phase_results":results,
            "red_wht_intersection_phase_range_deg2":[min(overlaps),max(overlaps)],
            "red_epoch_fraction_phase_range":[min(fractions),max(fractions)],
            "limitations":[
                "3.6-arcsec sampling with two offsets: phase differences are NOT confidence intervals.",
                "JWST star-mask intersection remains outstanding.",
                "Cross-filter epoch geometry does NOT establish F356W source detection or PM precision.",
                "WHT>0 is not F444W=26.2 AB dropout selection completeness.",
                "TAN WCS reconstructed from saved grid metadata; no full FITS is loaded."]}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--grid", default="pm_footprint/results/local_v13_sample_grid.npz")
    p.add_argument("--sample", default="pm_footprint/results/local_v13_sample.json")
    p.add_argument("--metadata", default="pm_footprint/results/mast_observations.csv")
    p.add_argument("--roster", default="pm_footprint/dave125_positions.csv")
    p.add_argument("--exact", default="pm_footprint/results/exact/summary.json")
    p.add_argument("--out", default="pm_footprint/results/local_v13_epoch_overlay.json")
    a = p.parse_args()
    result = run(a.grid,a.sample,a.metadata,a.roster,a.exact)
    path = Path(a.out); path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"status":result["status"],
                      "red_wht_intersection_deg2":result["red_wht_intersection_phase_range_deg2"],
                      "red_epoch_fraction":result["red_epoch_fraction_phase_range"]},indent=2))


if __name__ == "__main__":
    main()
