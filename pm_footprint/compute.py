"""Estimate multi-epoch JWST coverage in the COSMOS-Web F444W footprint.

This is geometric coverage, not PM detection or astrometric completeness.
Only public MAST metadata are queried; no FITS pixels are downloaded.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import re

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
from mocpy import MOC

PROGRAM = "1727"
DEPTH = 14
GAP_DAYS = 30.0
CENTER = (150.12, 2.20)
RED = {"F277W", "F356W", "F410M", "F444W"}


def empty(depth):
    return MOC.new_empty(max_depth=depth)


def deg2(moc):
    return float(moc.sky_fraction * (4 * math.pi) * (180 / math.pi) ** 2)


def filters(value):
    if value is None or pd.isna(value):
        return set()
    return set(re.findall(r"F\d{3,4}[A-Z]", str(value).upper()))


def polygon(region, depth):
    """Parse MAST CAOM s_region, including adjacent disjoint POLYGON groups.

    JWST CAOM observation footprints frequently encode multiple detector
    polygons as 'POLYGON <coords> POLYGON <coords>' without STC-S UNION.
    Treating that as a single polygon or discarding it biases sky coverage.
    """
    if region is None or pd.isna(region) or not str(region).strip():
        raise ValueError("missing CAOM s_region")
    source = str(region).strip()
    if not re.match(r"^POLYGON\\b", source, flags=re.I):
        return MOC.from_stcs(source, max_depth=depth)

    pieces = re.split(r"\\bPOLYGON\\b", source, flags=re.I)
    if pieces[0].strip() or len(pieces) < 2:
        raise ValueError("malformed multi-polygon CAOM region")
    combined = empty(depth)
    for body in pieces[1:]:
        body = re.sub(r"^\\s*ICRS\\s+", "", body, flags=re.I).strip()
        if not re.fullmatch(r"[+\\-\\d.eE\\s]+", body):
            raise ValueError("non-numeric polygon vertices: " + body[:90])
        v = [float(z) for z in body.split()]
        if len(v) < 6 or len(v) % 2:
            raise ValueError("bad CAOM polygon vertices")
        poly = MOC.from_polygon(
            lon=np.array(v[0::2]) * u.deg,
            lat=np.array(v[1::2]) * u.deg,
            max_depth=depth,
        )
        combined = combined.union(poly)
    return combined


def overlap(rows, depth, gap):
    """Union of regions imaged at least twice, separated by >= gap days."""
    rows = sorted(rows, key=lambda row: row[0])
    old = empty(depth)
    found = empty(depth)
    j = 0
    for t, region in rows:
        while j < len(rows) and rows[j][0] <= t - gap:
            old = old.union(rows[j][1])
            j += 1
        if not old.empty():
            found = found.union(region.intersection(old))
    return found


def f444_red_overlap(reds, f444s, depth, gap):
    """F444W at one epoch and a red NIRCam filter at another."""
    result = empty(depth)
    for direction in (1, -1):
        a = sorted(f444s, key=lambda r: direction * r["t"])
        b = sorted(reds, key=lambda r: direction * r["t"])
        older = empty(depth)
        j = 0
        for item in a:
            while j < len(b) and direction * (item["t"] - b[j]["t"]) >= gap:
                older = older.union(b[j]["moc"])
                j += 1
            if not older.empty():
                result = result.union(item["moc"].intersection(older))
    return result


def retrieve():
    from astroquery.mast import Observations
    obs = Observations.query_criteria(
        coordinates=f"{CENTER[0]} {CENTER[1]}",
        radius="0.9 deg",
        obs_collection="JWST",
        instrument_name="NIRCAM*",
        dataproduct_type="image",
        intentType="science",
    )
    if len(obs) == 0:
        raise RuntimeError("MAST returned no JWST/NIRCam observations")
    return obs.to_pandas()


def calculate(raw, candidate_df, depth=DEPTH, gap=GAP_DAYS):
    expected = {"obs_collection", "instrument_name", "filters", "proposal_id",
                "s_region", "t_min"}
    missing = expected - set(raw.columns)
    if missing:
        raise ValueError("Missing MAST CAOM columns: " + repr(sorted(missing)))
    if set(("candidate_id", "ra_deg", "dec_deg")) - set(candidate_df.columns):
        raise ValueError("Bad Dave candidate table")
    if len(candidate_df) != 125 or candidate_df.candidate_id.nunique() != 125:
        raise ValueError("Dave roster must have 125 distinct IDs")
    if not {387132, 1651050}.issubset(set(candidate_df.candidate_id)):
        raise ValueError("Missing the two established movers")

    frame = raw.copy()
    frame = frame[
        frame.obs_collection.astype(str).str.upper().eq("JWST")
        & frame.instrument_name.astype(str).str.upper().eq("NIRCAM/IMAGE")
    ].copy()
    if "intentType" in frame:
        frame = frame[frame.intentType.astype(str).str.lower().eq("science")]
    if "dataproduct_type" in frame:
        frame = frame[frame.dataproduct_type.astype(str).str.lower().eq("image")]
    frame["t_min"] = pd.to_numeric(frame["t_min"], errors="coerce")
    if "t_max" in frame:
        frame["t_max"] = pd.to_numeric(frame["t_max"], errors="coerce")
    else:
        frame["t_max"] = frame.t_min
    frame = frame[np.isfinite(frame.t_min)].copy()
    # A multi-date stacked mosaic is not one independent imaging epoch.
    multi_date = frame.t_max - frame.t_min > 2.0
    skipped_long = int(multi_date.sum())
    frame = frame[~multi_date]
    if frame.empty:
        raise RuntimeError("No usable timestamped single-visit image metadata")
    records = []
    errors = []
    for i, r in frame.iterrows():
        bands = filters(r["filters"])
        if not bands:
            continue
        try:
            footprint = polygon(r["s_region"], depth)
        except Exception as ex:
            errors.append({"index": str(i), "obs_id": str(r.get("obs_id", "")),
                           "reason": str(ex)})
            continue
        records.append({
            "t": float(r["t_min"]), "bands": bands,
            "pid": str(r["proposal_id"]).removesuffix(".0"),
            "moc": footprint,
        })
    # No silent geometrical undercount from dropped, unparsable footprints.
    if errors:
        raise RuntimeError("Unparsed footprint(s): " + json.dumps(errors[:8])
                           + f" [total={len(errors)}]")
    primary = [r for r in records if r["pid"] == PROGRAM and "F444W" in r["bands"]]
    if not primary:
        raise RuntimeError("No COSMOS-Web program 1727 F444W archive footprint")
    survey = empty(depth)
    for r in primary:
        survey = survey.union(r["moc"])
    if survey.empty() or not 0.1 < deg2(survey) < 1.5:
        raise RuntimeError(f"Implausible survey area {deg2(survey)} deg2")

    f444 = [r for r in records if "F444W" in r["bands"]]
    reds = [r for r in records if r["bands"] & RED]
    masks = {
        "cweb_f444w": survey,
        "same_f444w": overlap([(r["t"], r["moc"]) for r in f444], depth, gap).intersection(survey),
        "f444w_plus_red": f444_red_overlap(reds, f444, depth, gap).intersection(survey),
        "any_nircam": overlap([(r["t"], r["moc"]) for r in records], depth, gap).intersection(survey),
    }
    coords = SkyCoord(candidate_df.ra_deg.to_numpy() * u.deg,
                      candidate_df.dec_deg.to_numpy() * u.deg)
    members = candidate_df[["candidate_id", "ra_deg", "dec_deg"]].copy()
    for key, footprint in masks.items():
        members["in_" + key] = footprint.contains_skycoords(coords)
    area = {key: deg2(val) for key, val in masks.items()}
    base = area["cweb_f444w"]
    for key in masks:
        if area[key] > base * (1 + 1e-6):
            raise AssertionError("Mask exceeds survey base: " + key)
    if not (area["same_f444w"] <= area["f444w_plus_red"] + 1e-8
            <= area["any_nircam"] + 1e-8):
        raise AssertionError("Unexpected nested mask ordering")
    result = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "status": "GEOMETRIC_PRELIMINARY_REQUIRES_QC",
        "baseline_days": gap,
        "moc_max_depth": depth,
        "cosmosweb_program": PROGRAM,
        "mast_rows": len(raw),
        "single_visit_rows": len(frame),
        "excluded_multi_date_rows": skipped_long,
        "converted_footprints": len(records),
        "cosmosweb_f444w_rows": len(primary),
        "dave_candidates": len(candidate_df),
        "area_deg2": area,
        "fraction_of_archive_f444w_area": {key: area[key] / base for key in masks},
        "candidates_with_geometric_coverage": {
            key: int(members["in_" + key].sum()) for key in masks
        },
        "warnings": [
            "CAOM observation polygons approximate geometric coverage only.",
            "NIRCAM/GRISM spectral products are excluded from imaging opportunity.",
            "Actual F444W usable depth/masking denominator may differ from archive area.",
            "Other-filter overlap does NOT imply a red F444W-only candidate is detectable.",
            "No proper-motion detection/completeness efficiency is modeled.",
            "Compare HEALPix depth 14 to depth 15 and inspect visit-level metadata before quoting."
        ],
    }
    return result, members, masks


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--candidates", default="pm_footprint/dave125_positions.csv")
    p.add_argument("--out", default="pm_footprint/results")
    p.add_argument("--metadata-csv", default=None)
    p.add_argument("--depth", type=int, default=DEPTH)
    p.add_argument("--baseline-days", type=float, default=GAP_DAYS)
    opt = p.parse_args()
    if not 8 <= opt.depth <= 16 or opt.baseline_days < 1:
        p.error("Require 8 <= HEALPix depth <= 16 and baseline >= 1 day")
    out = Path(opt.out)
    out.mkdir(parents=True, exist_ok=True)
    try:
        table = pd.read_csv(opt.metadata_csv) if opt.metadata_csv else retrieve()
        cols = [c for c in ("obs_id", "obsid", "obs_collection", "instrument_name",
                            "dataproduct_type", "intentType", "filters", "proposal_id",
                            "t_min", "t_max", "s_region") if c in table]
        table[cols].to_csv(out / "mast_observations.csv", index=False)
        data = pd.read_csv(opt.candidates)
        result, members, masks = calculate(table, data, opt.depth, opt.baseline_days)
        members.to_csv(out / "dave125_coverage.csv", index=False)
        for key, mask in masks.items():
            mask.save(out / (key + ".fits"), format="fits", overwrite=True)
        (out / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
        print(json.dumps(result, indent=2))
    except Exception as e:
        failure = {"status": "FAILED_NO_AREA_CLAIM",
                   "utc": datetime.now(timezone.utc).isoformat(),
                   "error": repr(e)}
        (out / "failure.json").write_text(json.dumps(failure, indent=2) + "\n")
        raise


if __name__ == "__main__":
    main()
