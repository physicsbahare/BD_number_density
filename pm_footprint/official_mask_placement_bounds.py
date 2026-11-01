"""Worst-case mask-placement bounds for official COSMOS2025 WHT geometry.

Reads only frozen small JSON summaries; never downloads or rereads FITS.
These are nominal set-theoretic bounds, NOT a spatial intersection.
"""
import argparse
import json
from pathlib import Path


def fraction_bounds(survey_area, red_area, mask_cap):
    """Exact extremal fraction given ONLY an upper cap on excluded area.

    Assumes red subset of survey, mask may occur anywhere, and the cap
    is a true upper bound in the same units.
    """
    s, r, m = map(float, (survey_area, red_area, mask_cap))
    if not (s > 0 and 0 <= r <= s and 0 <= m < s):
        raise ValueError("Expected 0 <= red <= survey and 0 <= mask_cap < survey")
    red_removed = min(m, r)
    nonred_removed = min(m, s-r)
    lower = (r-red_removed)/(s-red_removed)
    upper = r/(s-nonred_removed)
    return lower, upper


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--wht-summary",
                   default="pm_footprint/results/wht60/summary.json")
    p.add_argument("--mask-summary",
                   default="pm_footprint/results/jwst_starmask/summary.json")
    p.add_argument("--output",
                   default="pm_footprint/results/jwst_starmask/official_wht_mask_placement_bounds_20261101.json")
    a = p.parse_args()
    wht = json.loads(Path(a.wht_summary).read_text())
    mask = json.loads(Path(a.mask_summary).read_text())
    m = mask["sum_per_tile_mask_code1_area_deg2_NOT_UNION"]
    s = wht["area_deg2"]["0"]["archive_f444w"]
    r = wht["area_deg2"]["0"]["f444w_plus_red_365d"]
    low, high = fraction_bounds(s, r, m)
    result = {
        "status": "NOMINAL_PIXEL_AREA_SET_BOUNDS_NOT_DIRECT_SPATIAL_INTERSECTION",
        "wht_survey_deg2": s,
        "wht_red365_deg2": r,
        "mask_per_tile_sum_nominal_deg2_upper_cap_if_code1_excluded": m,
        "unmasked_fraction": r/s,
        "fraction_lower": low,
        "fraction_upper": high,
        "assumptions": [
            "Mask code 1 denotes exclusion.",
            "Per-tile nominal 30mas pixel-area sum upper-bounds unique masked area; verify FITS PC determinant and TAN solid-angle correction.",
            "Red epoch region is a subset of official WHT>0 footprint."
        ],
        "interpretation": "Mask placement alone cannot move nominal geometric red fraction outside this interval, conditional on area cap; no actual mask-WHT intersection or dropout/PM completeness measured.",
        "independent_of_mask_sky_WCS": True,
        "repro": "python -m pm_footprint.official_mask_placement_bounds"
    }
    Path(a.output).write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
