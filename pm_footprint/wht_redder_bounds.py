"""Bound WHT-filter overlaps using only *committed* small COSMOS-Web products.

No FITS downloads, MAST queries, or assumptions about missing WHT pixels.
The existing WHT summary measures the all-red union, not separate filters.
Set theory bounds each redder-filter intersection with that gridded WHT map.
Candidate counts are directly joined from the two committed 125-row CSVs.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

FILTERS = ("F277W", "F356W", "F410M", "F444W", "red_without_f277w")
BASELINES = (30, 180, 365, 730)


def intersection_bounds(union_area, subset_area, union_inside_area):
    """For S subset B, bound |S intersect W| from |B|, |S| and |B intersect W|."""
    eps = 1e-9
    if (min(union_area, subset_area, union_inside_area) < -eps
            or subset_area > union_area + eps
            or union_inside_area > union_area + eps):
        raise ValueError("Invalid subset/union/WHT area relationship")
    return (max(0.0, subset_area - (union_area - union_inside_area)),
            min(subset_area, union_inside_area))


def load_rows(path):
    with Path(path).open(newline="") as f:
        rows = list(csv.DictReader(f))
    ids = [int(r["candidate_id"]) for r in rows]
    if len(ids) != 125 or len(set(ids)) != 125:
        raise ValueError(f"Expected exactly 125 unique candidate IDs in {path}")
    return {int(r["candidate_id"]): r for r in rows}


def flag(x):
    if x not in ("True", "False"):
        raise ValueError(f"Invalid boolean value: {x!r}")
    return x == "True"


def analyze(wht, filt, wht_rows, filter_rows, external_usable_area=0.49672):
    if set(wht_rows) != set(filter_rows):
        raise ValueError("WHT and filter candidate rosters differ")
    result = {}
    for threshold in ("0", "0.1", "0.5"):
        A = float(wht["area_deg2"][threshold]["archive_f444w"])
        if A <= 0:
            raise ValueError("WHT area must be positive")
        threshold_result = {}
        for gap in BASELINES:
            exact = filt["baselines"][str(gap)]
            B = float(exact["area_deg2"]["red_union"])
            W_B = float(wht["area_deg2"][threshold][f"f444w_plus_red_{gap}d"])
            expected_wht_count = int(
                wht["candidate_counts"][threshold][f"f444w_plus_red_{gap}d"])
            observed_red = 0
            subset_count = {f: 0 for f in FILTERS}
            exact_red_ids = set(exact["candidate_ids"]["red_union"])
            csv_red_ids = set()
            for id_, row in wht_rows.items():
                fr = filter_rows[id_]
                in_wht = float(row[f"wht_{threshold}_fraction"]) > 0
                in_red = flag(row[f"in_f444w_plus_red_{gap}d"])
                if flag(fr[f"in_red_union_{gap}d"]):
                    csv_red_ids.add(id_)
                if in_wht and in_red:
                    observed_red += 1
                for f in FILTERS:
                    if in_wht and in_red and flag(fr[f"in_{f.lower()}_{gap}d"]):
                        subset_count[f] += 1
            if csv_red_ids != exact_red_ids:
                raise ValueError(f"{gap}d filter CSV red IDs mismatch exact summary")
            if observed_red != expected_wht_count:
                raise ValueError(f"{threshold} {gap}d WHT red count mismatch")
            if len(exact_red_ids) < observed_red:
                raise ValueError("WHT red count cannot exceed exact red count")
            masks = {}
            for f in FILTERS:
                area = float(exact["area_deg2"][f])
                lo, hi = intersection_bounds(B, area, W_B)
                masks[f] = {
                    "wht_intersection_area_bounds_deg2": [lo, hi],
                    "fraction_of_wht_area_bounds": [lo / A, hi / A],
                    "candidate_count_from_existing_memberships": subset_count[f],
                }
            threshold_result[str(gap)] = {
                "all_red_wht_intersection_deg2": W_B,
                "all_red_wht_candidate_count": observed_red,
                "filter_intersection_bounds": masks,
            }
            if threshold == "0" and gap == 365:
                lo, hi = masks["red_without_f277w"]["wht_intersection_area_bounds_deg2"]
                U = float(external_usable_area)
                if not 0 < U <= A:
                    raise ValueError("External usable area must be subset of WHT area")
                threshold_result[str(gap)]["conditional_external_usable_area"] = {
                    "assumed_area_deg2": U,
                    "assumption": "External usable area is an unknown spatial subset of gridded WHT>0",
                    "redder_fraction_bounds": [
                        max(0, lo - (A - U)) / U,
                        min(hi, U) / U,
                    ],
                }
        result[threshold] = threshold_result
    return result


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--wht-summary", default="pm_footprint/results/wht60/summary.json")
    p.add_argument("--filter-summary", default="pm_footprint/results/filter_epoch_qc/summary.json")
    p.add_argument("--wht-candidates", default="pm_footprint/results/wht60/dave125_wht_membership.csv")
    p.add_argument("--filter-candidates", default="pm_footprint/results/filter_epoch_qc/candidate_filter_membership.csv")
    p.add_argument("--output", default="pm_footprint/results/wht60/redder_filter_set_bounds.json")
    args = p.parse_args()
    paths = {
        "wht_summary": Path(args.wht_summary),
        "filter_summary": Path(args.filter_summary),
        "wht_candidates": Path(args.wht_candidates),
        "filter_candidates": Path(args.filter_candidates),
    }
    wht = json.loads(paths["wht_summary"].read_text())
    filt = json.loads(paths["filter_summary"].read_text())
    results = analyze(wht, filt, load_rows(paths["wht_candidates"]),
                      load_rows(paths["filter_candidates"]))
    report = {
        "status": "SET_THEORY_BOUNDS_CONDITIONAL_ON_EXISTING_RASTER_WHT_NOT_EXACT_INTERSECTION",
        "method": "Set-inclusion area bounds; exact 125-source CSV membership join",
        "inputs_sha256": {
            k: hashlib.sha256(v.read_bytes()).hexdigest() for k, v in paths.items()
        },
        "results": results,
        "limitations": [
            "No filter-specific WHT spatial intersection was computed; report intervals, not point estimates.",
            "The 0.6arcsec WHT gridding loses about 0.96% native tile area; bounds are conditional on stored rasterized areas.",
            "External 0.49672deg2 usable area is not a validated spatial mask; its bounds require a subset assumption.",
            "WHT>0 and cross-filter imaging are not F444W=26.2 dropout completeness or PM measurement completeness.",
            "Candidate counts refer to catalog entries, not deblended independent objects.",
        ],
    }
    dest = Path(args.output)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["results"]["0"]["365"], indent=2))


if __name__ == "__main__":
    main()
