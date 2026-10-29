"""Independent exact polygon-intersection coverage audit (no MOC rasterization).

Uses a unit-sphere Lambert azimuthal equal-area projection centered on COSMOS.
This yields steradian areas from Shapely polygon unions and is independent of
the HEALPix/MOC implementation. CAOM image footprints are geometric proxies,
not proof of detectable or precise proper motion.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from shapely import STRtree, make_valid
from shapely.geometry import Point, Polygon
from shapely.ops import unary_union

CENTER_RA, CENTER_DEC = 150.12, 2.20
COS_DEC0 = math.cos(math.radians(CENTER_DEC))
SIN_DEC0 = math.sin(math.radians(CENTER_DEC))
RA0 = math.radians(CENTER_RA)
DEG2_PER_SR = (180.0 / math.pi) ** 2
BASELINES = (30, 180, 365, 730)
RED = {"F277W", "F356W", "F410M", "F444W"}


def lambert(lon_deg, lat_deg):
    """Lambert azimuthal equal-area coordinates on the unit celestial sphere."""
    phi = np.deg2rad(np.asarray(lat_deg, dtype=float))
    dlon = np.deg2rad(np.asarray(lon_deg, dtype=float)) - RA0
    sinphi, cosphi = np.sin(phi), np.cos(phi)
    c = SIN_DEC0 * sinphi + COS_DEC0 * cosphi * np.cos(dlon)
    if np.any(c <= -0.99):
        raise ValueError("Projection far outside COSMOS field")
    k = np.sqrt(2.0 / (1.0 + c))
    return (k * cosphi * np.sin(dlon),
            k * (COS_DEC0 * sinphi - SIN_DEC0 * cosphi * np.cos(dlon)))


def polygon_region(s):
    """Parse all disjoint POLYGON groups in a MAST CAOM s_region."""
    if not isinstance(s, str) or not s.strip():
        raise ValueError("Missing CAOM footprint")
    chunks = re.split(r"\bPOLYGON\b", s, flags=re.I)
    if chunks[0].strip() or len(chunks) < 2:
        raise ValueError("Unsupported CAOM region: " + s[:90])
    polys = []
    for chunk in chunks[1:]:
        chunk = re.sub(r"^\s*ICRS\s+", "", chunk, flags=re.I).strip()
        if not re.fullmatch(r"[+\-\d.eE\s]+", chunk):
            raise ValueError("Non-numeric CAOM polygon: " + chunk[:90])
        values = np.fromstring(chunk, sep=" ")
        if len(values) < 6 or len(values) % 2:
            raise ValueError("Invalid polygon vertex count")
        x, y = lambert(values[::2], values[1::2])
        poly = Polygon(np.column_stack([x, y]))
        if not poly.is_valid:
            poly = make_valid(poly)
        if poly.is_empty or poly.area <= 0:
            raise ValueError("Empty CAOM polygon")
        polys.append(poly)
    return unary_union(polys)


def bands(value):
    return set(re.findall(r"F\d{3,4}[A-Z]", str(value).upper()))


def compute(raw, roster, baselines=BASELINES):
    if len(roster) != 125 or roster.candidate_id.nunique() != 125:
        raise ValueError("Expected 125 distinct Dave candidate IDs")
    if not {387132, 1651050}.issubset(set(roster.candidate_id)):
        raise ValueError("Both known moving-source IDs must be present")
    required = {"obs_collection", "instrument_name", "dataproduct_type",
                "intentType", "proposal_id", "filters", "t_min", "t_max", "s_region"}
    if required - set(raw.columns):
        raise ValueError("Missing MAST columns: " + str(required - set(raw.columns)))
    frame = raw.loc[
        raw.obs_collection.astype(str).str.upper().eq("JWST")
        & raw.instrument_name.astype(str).str.upper().eq("NIRCAM/IMAGE")
        & raw.dataproduct_type.astype(str).str.lower().eq("image")
        & raw.intentType.astype(str).str.lower().eq("science")
    ].copy()
    frame["t_min"] = pd.to_numeric(frame.t_min, errors="coerce")
    frame["t_max"] = pd.to_numeric(frame.t_max, errors="coerce")
    frame = frame[np.isfinite(frame.t_min) & np.isfinite(frame.t_max)
                  & (frame.t_max >= frame.t_min)].reset_index(drop=True)
    if len(frame) == 0:
        raise ValueError("No valid imaging records")
    polygons = [polygon_region(s) for s in frame.s_region]
    f444 = np.array(["F444W" in bands(v) for v in frame.filters], dtype=bool)
    red = np.array([bool(bands(v) & RED) for v in frame.filters], dtype=bool)
    primary = (frame.proposal_id.astype(str).str.replace(r"\.0$", "", regex=True)
               .eq("1727").to_numpy() & f444)
    if primary.sum() < 10:
        raise ValueError("COSMOS-Web F444W program-1727 footprint not found")
    survey = unary_union([polygons[i] for i in np.where(primary)[0]])
    survey_area = float(survey.area * DEG2_PER_SR)
    if not 0.1 < survey_area < 1.5:
        raise ValueError("Implausible F444W survey area")

    tree = STRtree(polygons)
    i, j = tree.query(polygons, predicate="intersects")
    distinct = i < j
    i, j = i[distinct], j[distinct]
    starts = frame.t_min.to_numpy()
    ends = frame.t_max.to_numpy()
    # Positive only when observation time intervals are disjoint.
    separation = np.maximum(starts[i] - ends[j], starts[j] - ends[i])
    xs, ys = lambert(roster.ra_deg.to_numpy(), roster.dec_deg.to_numpy())
    pts = [Point(float(x), float(y)) for x, y in zip(xs, ys)]
    survey_members = np.array([survey.covers(p) for p in pts], dtype=bool)
    membership = roster[["candidate_id", "ra_deg", "dec_deg"]].copy()
    membership["in_archive_f444w"] = survey_members
    if not survey_members.all():
        raise ValueError("Dave candidate coordinates outside archival survey mask")

    summary = {
        "method": "unit_sphere_lambert_equal_area_shapely_polygon_intersections",
        "status": "EXACT_GEOMETRY_PRELIMINARY_REQUIRES_USABLE_MASK_AND_ASTROMETRY_QC",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "mast_rows": int(len(raw)),
        "imaging_rows_with_valid_time": int(len(frame)),
        "cweb_f444w_rows": int(primary.sum()),
        "archive_f444w_area_deg2": survey_area,
        "dave_candidate_count": int(len(roster)),
        "baseline_results": {},
        "warnings": [
            "MAST observation polygons approximate geometric coverage, not calibrated detection completeness.",
            "NIRCAM/GRISM spectra excluded from imaging coverage.",
            "Conservative temporal gap uses nonoverlapping [t_min,t_max] intervals.",
            "F444W+red and any-NIRCam do not guarantee detection of F444W-only dropout sources.",
            "A real COSMOS-Web F444W usable-depth/masking footprint is needed for selection-matched fractions.",
            "Candidate deblending and source-level astrometric S/N require independent checks.",
        ],
    }
    for gap in baselines:
        enough = separation >= gap
        selections = {
            "same_f444w": enough & f444[i] & f444[j],
            "f444w_plus_red": enough & ((f444[i] & red[j]) | (f444[j] & red[i])),
            "any_nircam": enough,
        }
        area = {}
        members = {}
        ids = {}
        for key, sel in selections.items():
            parts = [polygons[int(ii)].intersection(polygons[int(jj)])
                     for ii, jj in zip(i[sel], j[sel])]
            mask = unary_union(parts).intersection(survey) if parts else Polygon()
            area[key] = float(mask.area * DEG2_PER_SR)
            contained = np.array([mask.covers(p) for p in pts], dtype=bool)
            members[key] = int(contained.sum())
            ids[key] = [int(v) for v in roster.candidate_id[contained]]
            membership[f"in_{key}_{gap}d"] = contained
        if not area["same_f444w"] <= area["f444w_plus_red"] + 1e-10 <= area["any_nircam"] + 1e-10:
            raise AssertionError("Masks not nested")
        summary["baseline_results"][str(gap)] = {
            "area_deg2": area,
            "fraction_of_archive_f444w": {key: v / survey_area for key, v in area.items()},
            "candidate_counts": members,
            "candidate_ids": ids,
            "confirmed_pm_ids_in_red_mask": [
                x for x in (387132, 1651050) if x in ids["f444w_plus_red"]
            ],
        }
    return summary, membership


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--metadata-csv", default="pm_footprint/results/mast_observations.csv")
    p.add_argument("--candidates", default="pm_footprint/dave125_positions.csv")
    p.add_argument("--out", default="pm_footprint/results/exact")
    a = p.parse_args()
    meta, cand = Path(a.metadata_csv), Path(a.candidates)
    raw, roster = pd.read_csv(meta), pd.read_csv(cand)
    summary, membership = compute(raw, roster)
    summary["exact_source_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    summary["mast_metadata_sha256"] = hashlib.sha256(meta.read_bytes()).hexdigest()
    summary["candidate_catalog_sha256"] = hashlib.sha256(cand.read_bytes()).hexdigest()
    dest = Path(a.out)
    dest.mkdir(parents=True, exist_ok=True)
    (dest / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    membership.to_csv(dest / "dave125_exact_membership.csv", index=False)
    print(json.dumps({
        "status": summary["status"],
        "archive_f444w_area_deg2": summary["archive_f444w_area_deg2"],
        "baselines": {
            k: {"area_deg2": v["area_deg2"], "candidate_counts": v["candidate_counts"]}
            for k, v in summary["baseline_results"].items()
        }
    }, indent=2))


if __name__ == "__main__":
    main()
