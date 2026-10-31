"""Audit conservative F444W/F356W imaging baselines for frozen Dave125 candidates.

Run from repository root:
 python -m pm_footprint.pm_f356w_baseline_audit
Uses frozen PUBLIC CAOM metadata, not new MAST/WHT downloads. These are
geometric opportunities, NOT F356W detections or PM completeness.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd

FILTERS = ("F277W", "F356W", "F410M", "F444W")
BASELINES = (30, 180, 365, 730, 1000)


def clean_observations(raw):
    f = raw.loc[
        raw.obs_collection.astype(str).str.upper().eq("JWST")
        & raw.instrument_name.astype(str).str.upper().eq("NIRCAM/IMAGE")
        & raw.dataproduct_type.astype(str).str.lower().eq("image")
        & raw.intentType.astype(str).str.lower().eq("science")
        & raw.dataRights.astype(str).str.upper().eq("PUBLIC")
    ].copy()
    f["t_min"] = pd.to_numeric(f.t_min, errors="coerce")
    f["t_max"] = pd.to_numeric(f.t_max, errors="coerce")
    return f[np.isfinite(f.t_min) & np.isfinite(f.t_max) & (f.t_max >= f.t_min)].reset_index(drop=True)


def best_interval_pair(records, other_filter):
    """Maximize disjoint-interval gap, not midpoint separation (MJD days)."""
    from .exact_geometry import bands
    f444 = [(i, r) for i, r in enumerate(records) if "F444W" in bands(r["filters"])]
    other = [(j, r) for j, r in enumerate(records) if other_filter in bands(r["filters"])]
    best = None
    for i, a in f444:
        for j, b in other:
            if i == j or str(a["obsid"]) == str(b["obsid"]):
                continue
            gap = max(float(a["t_min"]) - float(b["t_max"]),
                      float(b["t_min"]) - float(a["t_max"]))
            if gap < 0:
                continue
            if best is None or gap > best["gap_days"]:
                best = {
                    "gap_days": float(gap),
                    "midpoint_days": abs((float(a["t_min"]) + float(a["t_max"])) / 2
                                         - (float(b["t_min"]) + float(b["t_max"])) / 2),
                    "f444_obsid": str(a["obsid"]),
                    "other_obsid": str(b["obsid"]),
                    "f444_proposal": str(a["proposal_id"]),
                    "other_proposal": str(b["proposal_id"]),
                }
    return best


def candidate_report(candidate_id, records):
    from .exact_geometry import bands
    primary = [r for r in records if "F444W" in bands(r["filters"])
               and str(r["proposal_id"]).removesuffix(".0") == "1727"]
    best = {f: best_interval_pair(records, f) for f in FILTERS}
    return {"candidate_id": int(candidate_id), "n_primary_cosmosweb_f444w": len(primary),
            "best": best}


def collect_records(raw, roster):
    from shapely import STRtree, points
    from .exact_geometry import polygon_region, lambert
    f = clean_observations(raw)
    polygons = [polygon_region(s) for s in f.s_region]
    x, y = lambert(roster.ra_deg.to_numpy(), roster.dec_deg.to_numpy())
    i, j = STRtree(polygons).query(points(x, y), predicate="intersects")
    matches = [[] for _ in range(len(roster))]
    for ii, jj in zip(i, j):
        matches[int(ii)].append(int(jj))
    return f, [[f.iloc[j].to_dict() for j in js] for js in matches]


def analyze(raw, roster, wht_membership=None):
    if len(roster) != 125 or roster.candidate_id.nunique() != 125:
        raise ValueError("Expected exactly 125 distinct Dave IDs")
    f, matches = collect_records(raw, roster)
    reports = [candidate_report(cid, rows) for cid, rows in zip(roster.candidate_id, matches)]
    if any(r["n_primary_cosmosweb_f444w"] == 0 for r in reports):
        raise ValueError("Candidate without primary COSMOS-Web F444W coverage")
    counts = {
        str(g): {
            **{band: sum(bool(r["best"][band] and r["best"][band]["gap_days"] >= g)
                         for r in reports) for band in FILTERS},
            "red_union": sum(any(r["best"][band] and r["best"][band]["gap_days"] >= g
                                 for band in FILTERS) for r in reports),
        } for g in BASELINES
    }
    selected = [r for r in reports if r["best"]["F356W"]
                and r["best"]["F356W"]["gap_days"] >= 365]
    gaps = np.array([r["best"]["F356W"]["gap_days"] for r in selected], dtype=float)
    if wht_membership is not None:
        for gap in (365, 730):
            col = f"in_f444w_plus_red_{gap}d"
            expected = set(wht_membership.loc[wht_membership[col].astype(str).eq("True"),
                                               "candidate_id"].astype(int))
            observed = {r["candidate_id"] for r in reports if any(
                r["best"][band] and r["best"][band]["gap_days"] >= gap for band in FILTERS)}
            if observed != expected:
                raise AssertionError(f"{gap}d WHT membership mismatch: {sorted(observed ^ expected)}")
    movers = {str(r["candidate_id"]): r["best"]["F356W"] for r in reports
              if r["candidate_id"] in (387132, 1651050)}
    summary = {
        "status": "F356W_EPOCH_BASELINES_VERIFIED_GEOMETRIC_ONLY",
        "n_frozen_candidates": len(reports),
        "n_valid_public_nircam_imaging_rows": len(f),
        "counts_by_conservative_gap_days": counts,
        "f356w_365d_n_candidates": len(selected),
        "f356w_365d_gap_days_quantiles": {
            "min": float(gaps.min()), "q25": float(np.quantile(gaps, .25)),
            "median": float(np.median(gaps)), "q75": float(np.quantile(gaps, .75)),
            "max": float(gaps.max()),
        },
        "f356w_365d_other_proposal_counts": {
            p: sum(r["best"]["F356W"]["other_proposal"] == p for r in selected)
            for p in sorted({r["best"]["F356W"]["other_proposal"] for r in selected})
        },
        "confirmed_movers_f356w": movers,
        "wht_candidate_membership_crosschecked": wht_membership is not None,
        "illustrative_3sigma_pm_threshold_mas_per_year_for_20mas_independent_per_epoch_error": {
            "median_baseline": float(3*np.sqrt(2)*20/(np.median(gaps)/365.25)),
            "mover_387132": float(3*np.sqrt(2)*20/(movers["387132"]["gap_days"]/365.25)),
            "mover_1651050": float(3*np.sqrt(2)*20/(movers["1651050"]["gap_days"]/365.25)),
        },
        "caveats": [
            "PUBLIC MAST CAOM polygons/times are geometric opportunity, not F356W flux/SNR, centroid accuracy, or PM completeness.",
            "Longest available temporal gap need not be highest-quality astrometric epoch.",
            "F444W/F356W are different filters; chromatic centroid biases require calibration.",
            "3sigma threshold is illustrative for independent 20mas errors per epoch, no systematics; not measured completeness.",
            "No WHT/starmask spatial rerun, no new archive download; user starmask binary not in GitHub.",
        ],
    }
    rows = []
    for r in reports:
        b = r["best"]["F356W"]
        rows.append({
            "candidate_id": r["candidate_id"], "n_primary_f444w": r["n_primary_cosmosweb_f444w"],
            "best_f356w_gap_days": b["gap_days"] if b else None,
            "best_f356w_midpoint_days": b["midpoint_days"] if b else None,
            "f444_obsid": b["f444_obsid"] if b else None,
            "f356_obsid": b["other_obsid"] if b else None,
            "f444_proposal": b["f444_proposal"] if b else None,
            "f356_proposal": b["other_proposal"] if b else None,
            "f356w_ge365": bool(b and b["gap_days"] >= 365),
            "f356w_ge730": bool(b and b["gap_days"] >= 730),
        })
    return summary, pd.DataFrame(rows)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--metadata", default="pm_footprint/results/mast_observations.csv")
    p.add_argument("--candidates", default="pm_footprint/dave125_positions.csv")
    p.add_argument("--wht-membership", default="pm_footprint/results/wht60/dave125_wht_membership.csv")
    p.add_argument("--output-dir", default="pm_footprint/results/f356w_baselines")
    a = p.parse_args()
    wht = pd.read_csv(a.wht_membership) if Path(a.wht_membership).exists() else None
    summary, table = analyze(pd.read_csv(a.metadata), pd.read_csv(a.candidates), wht)
    out = Path(a.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    table.to_csv(out / "candidate_f356w_baselines.csv", index=False)
    print(json.dumps({"n": summary["f356w_365d_n_candidates"],
                      "median_gap_days": summary["f356w_365d_gap_days_quantiles"]["median"]}, indent=2))


if __name__ == "__main__":
    main()
