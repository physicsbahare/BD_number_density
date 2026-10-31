"""Worst-case WHT geometric epoch-fraction bounds without WHT-mask pixels.

For WHT-positive area W, red-epoch subset WR, and a spatial exclusion
mask with unique sky area <= M, distribute at most M between red/nonred
regions. No Monte Carlo or mask raster is needed for these bounds.
"""
import argparse
import json
from pathlib import Path

def bounds(W, WR, M):
    W, WR, M = map(float, (W, WR, M))
    if not (0 < WR < W and 0 <= M < min(WR, W-WR)):
        raise ValueError("Require 0<=M<min(WR,W-WR) and 0<WR<W")
    return {
        "lower": (WR-M)/(W-M),
        "upper": WR/(W-M),
        "uncorrected": WR/W,
        "max_mask_deg2": M,
        "meaning": "Conditional on WHT geometric areas; not selection or PM completeness.",
    }

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--wht", default="pm_footprint/results/wht60/summary.json")
    p.add_argument("--mask", default="pm_footprint/results/jwst_starmask/summary.json")
    p.add_argument("--output", required=True)
    a = p.parse_args()
    w = json.loads(Path(a.wht).read_text())
    mask = json.loads(Path(a.mask).read_text())
    W = w["area_deg2"]["0"]["archive_f444w"]
    WR = w["area_deg2"]["0"]["f444w_plus_red_365d"]
    M = mask["sum_per_tile_mask_code1_area_deg2_NOT_UNION"]
    out = {
        "status": "DETERMINISTIC_WORST_CASE_STAR_MASK_GEOMETRIC_BOUNDS",
        "source": "Committed WHT-positive areas and all 20 user-supplied JWST mask code1 counts; no WHT or MAST rerun",
        "wht_area_deg2": W,
        "wht_red365_area_deg2": WR,
        "native_mask_pixel_sum_deg2_upper_bound_unique_sky": M,
        "why_upper_bound": "Union area <= sum of per-tile flagged pixel areas; TAN solid-angle Jacobian <= nominal tangent-plane pixel area for unit-absolute-determinant PC and 0.03 arcsec CDELT.",
        "bounds": bounds(W, WR, M),
        "caveats": [
            "Requires code1=excluded; confirm official binary-mask semantics.",
            "Mathematical bounds conditional on previously gridded WHT area estimates with ~0.96% deposition loss.",
            "No WHT-star-mask spatial intersection, local blue/F444W depth, source quality or F356W SNR measured.",
            "No F444W=26.2 selection or proper-motion completeness.",
        ],
    }
    Path(a.output).write_text(json.dumps(out, indent=2) + "\n")

if __name__ == "__main__":
    main()
