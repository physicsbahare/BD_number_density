"""Filter-resolved *geometrical* JWST epoch opportunities for Dave125.

Uses ONLY the frozen PUBLIC-only MAST metadata and Dave125 positions.
F444W+F277W, F444W+F356W, F444W+F410M, F444W+F444W
are computed separately: a very red F444W-only dropout need not be
detected in F277W. None of these masks measures PM completeness.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from shapely import STRtree
from shapely.geometry import Point, Polygon
from shapely.ops import unary_union

from .exact_geometry import DEG2_PER_SR, BASELINES, bands, lambert, polygon_region

FILTERS = ("F277W", "F356W", "F410M", "F444W")
NONBLUE = ("F356W", "F410M", "F444W")


def filter_epoch_masks(raw, roster, baselines=BASELINES):
    """Return filter-resolved exact masks and a positional membership table."""
    required = {"obs_collection", "instrument_name", "dataproduct_type",
                "intentType", "dataRights", "proposal_id", "filters",
                "t_min", "t_max", "s_region"}
    if required - set(raw.columns):
        raise ValueError("Missing PUBLIC MAST columns: " + repr(sorted(required - set(raw.columns))))
    if len(roster) != 125 or roster.candidate_id.nunique() != 125:
        raise ValueError("Expected frozen 125 unique Dave candidate IDs")
    if not {387132, 1651050}.issubset(set(roster.candidate_id)):
        raise ValueError("Missing confirmed movers")
    frame = raw.loc[
        raw.obs_collection.astype(str).str.upper().eq("JWST")
        & raw.instrument_name.astype(str).str.upper().eq("NIRCAM/IMAGE")
        & raw.dataproduct_type.astype(str).str.lower().eq("image")
        & raw.intentType.astype(str).str.lower().eq("science")
        & raw.dataRights.astype(str).str.upper().eq("PUBLIC")
    ].copy()
    frame["t_min"] = pd.to_numeric(frame.t_min, errors="coerce")
    frame["t_max"] = pd.to_numeric(frame.t_max, errors="coerce")
    frame = frame.loc[np.isfinite(frame.t_min) & np.isfinite(frame.t_max)
                      & (frame.t_max >= frame.t_min)].reset_index(drop=True)
    if frame.empty:
        raise ValueError("No valid PUBLIC NIRCam imaging metadata")
    polys = [polygon_region(x) for x in frame.s_region]
    fsets = [bands(v) for v in frame.filters]
    f444 = np.array(["F444W" in b for b in fsets])
    base = (frame.proposal_id.astype(str).str.replace(r"\.0$", "", regex=True)
            .eq("1727").to_numpy() & f444)
    if not base.any():
        raise ValueError("Missing COSMOS-Web program 1727 F444W baseline")
    survey = unary_union([polys[k] for k in np.where(base)[0]])
    tree = STRtree(polys)
    ii, jj = tree.query(polys, predicate="intersects")
    unique = ii < jj
    ii, jj = ii[unique], jj[unique]
    t0, t1 = frame.t_min.to_numpy(), frame.t_max.to_numpy()
    gap_days = np.maximum(t0[ii] - t1[jj], t0[jj] - t1[ii])
    xx, yy = lambert(roster.ra_deg.to_numpy(), roster.dec_deg.to_numpy())
    pts = [Point(float(x), float(y)) for x, y in zip(xx, yy)]
    if not all(survey.covers(p) for p in pts):
        raise ValueError("Dave roster contains positions outside archive F444W")
    members = roster[["candidate_id", "ra_deg", "dec_deg"]].copy()
    result = {}
    for gap in baselines:
        eligible = gap_days >= gap
        masks = {}
        for f in FILTERS:
            has = np.array([f in b for b in fsets])
            sel = eligible & ((f444[ii] & has[jj]) | (f444[jj] & has[ii]))
            pieces = [polys[int(i)].intersection(polys[int(j)])
                      for i, j in zip(ii[sel], jj[sel])]
            masks[f] = unary_union(pieces).intersection(survey) if pieces else Polygon()
        red_union = unary_union([masks[f] for f in FILTERS])
        nonblue = unary_union([masks[f] for f in NONBLUE])
        masks["red_union"] = red_union
        masks["red_without_f277w"] = nonblue
        counts = {}
        area = {}
        ids = {}
        for f, mask in masks.items():
            col = f"in_{f.lower()}_{gap}d"
            inside = np.array([mask.covers(p) for p in pts], dtype=bool)
            members[col] = inside
            counts[f] = int(inside.sum())
            ids[f] = [int(x) for x in roster.candidate_id[inside]]
            area[f] = float(mask.area * DEG2_PER_SR)
        if area["red_without_f277w"] > area["red_union"] + 1e-10:
            raise AssertionError("Red nonblue union exceeds red union")
        if counts["red_without_f277w"] > counts["red_union"]:
            raise AssertionError("Nonblue candidate count exceeds red union")
        # The all-red set includes F277W, a *non-detection* band in Dave125.
        # Count catalogue entries that have NO redder second-epoch option.
        only_f277w = sorted(set(ids["red_union"]) - set(ids["red_without_f277w"]))
        f277w_only_area = max(0.0, area["red_union"] - area["red_without_f277w"])
        result[str(gap)] = {
            "only_f277w_area_deg2": f277w_only_area,
            "only_f277w_candidate_count": len(only_f277w),
            "only_f277w_candidate_ids": only_f277w,
            "confirmed_movers_only_f277w": [
                x for x in (387132, 1651050) if x in only_f277w
            ],
            "area_deg2": area,
            "fraction_of_archive_f444w": {f: a / (survey.area * DEG2_PER_SR)
                                           for f, a in area.items()},
            "candidate_counts": counts,
            "candidate_ids": ids,
            "confirmed_movers_by_filter": {
                f: [x for x in (387132, 1651050) if x in ids[f]] for f in masks
            },
        }
    return result, members, float(survey.area * DEG2_PER_SR)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", default="pm_footprint/results/mast_observations.csv")
    parser.add_argument("--candidates", default="pm_footprint/dave125_positions.csv")
    parser.add_argument("--exact-summary", default="pm_footprint/results/exact/summary.json")
    parser.add_argument("--out", default="pm_footprint/results/filter_epoch_qc")
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    meta_path, cand_path = Path(args.metadata), Path(args.candidates)
    exact_path = Path(args.exact_summary)
    try:
        meta = pd.read_csv(meta_path)
        roster = pd.read_csv(cand_path)
        result, members, base_area = filter_epoch_masks(meta, roster)
        previous = json.loads(exact_path.read_text())
        # A matching area alone is insufficient if the public CAOM or Dave
        # roster snapshot changed since the reference exact-geometry run.
        if hashlib.sha256(meta_path.read_bytes()).hexdigest() != previous["mast_metadata_sha256"]:
            raise RuntimeError("PUBLIC MAST snapshot differs from frozen exact run")
        if hashlib.sha256(cand_path.read_bytes()).hexdigest() != previous["candidate_catalog_sha256"]:
            raise RuntimeError("Dave125 roster differs from frozen exact run")
        if abs(base_area - previous["archive_f444w_area_deg2"]) > 1e-8:
            raise RuntimeError("Archival footprint changed: source/metadata mismatch")
        for gap, entry in result.items():
            earlier = previous["baseline_results"][gap]
            area = entry["area_deg2"]["red_union"]
            expected = earlier["area_deg2"]["f444w_plus_red"]
            n = entry["candidate_counts"]["red_union"]
            n_expected = earlier["candidate_counts"]["f444w_plus_red"]
            expected_ids = set(earlier["candidate_ids"]["f444w_plus_red"])
            actual_ids = set(entry["candidate_ids"]["red_union"])
            if abs(area - expected) > 1e-8 or n != n_expected or actual_ids != expected_ids:
                raise RuntimeError(f"{gap}d filter union mismatch with frozen exact mask: "
                                   f"area={area}, expected={expected}; counts={n},{n_expected}; "
                                   f"IDs only new={sorted(actual_ids - expected_ids)}, "
                                   f"only old={sorted(expected_ids - actual_ids)}")
        report = {
            "status": "PUBLIC_FILTER_RESOLVED_GEOMETRIC_OPPORTUNITY_ONLY",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "archive_area_deg2": base_area,
            "baselines": result,
            "metadata_sha256": hashlib.sha256(meta_path.read_bytes()).hexdigest(),
            "candidates_sha256": hashlib.sha256(cand_path.read_bytes()).hexdigest(),
            "exact_summary_sha256": hashlib.sha256(exact_path.read_bytes()).hexdigest(),
            "warnings": [
                "F277W is a blue non-detection band in Dave's F444W-only selection: its geometric overlap may be astrophysically unusable for PM.",
                "F356W/F410M/F444W geometric overlap is a more relevant but still permissive red-source imaging opportunity, not an astrometric completeness estimate.",
                "These are PUBLIC-only exact archival polygons, not a 26.2 AB selection-weighted WHT area.",
                "Some close candidate pairs may be blends; do not interpret catalog entries as independent brown dwarfs.",
                "Cross-filter astrometry additionally requires adequate SNR and distortion registration."
            ],
        }
        members.to_csv(out / "candidate_filter_membership.csv", index=False)
        (out / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps({
            "status": report["status"],
            "365d_area_deg2": result["365"]["area_deg2"],
            "365d_candidate_counts": result["365"]["candidate_counts"]
        }, indent=2))
    except Exception as exc:
        (out / "failure.json").write_text(json.dumps({
            "status": "FAILED_NO_FILTER_COVERAGE_CLAIM", "error": repr(exc),
            "utc": datetime.now(timezone.utc).isoformat(),
        }, indent=2) + "\n")
        raise


if __name__ == "__main__":
    main()
