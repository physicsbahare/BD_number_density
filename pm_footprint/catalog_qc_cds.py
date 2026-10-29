"""Small public CDS/VizieR COSMOS2025 catalog audit of Dave's 125 positions.

This queries ONLY candidate rows, not the 784,016-row full catalog and not
the multi-GB WHT mosaics. It audits catalog flags, not spatial completeness.
CDS catalog J/A+A/704/A339/phot is a separate public mirror of COSMOS2025;
check the catalog release/version and coordinates before accepting ID joins.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from io import StringIO
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd
import requests

TABLE = '"J/A+A/704/A339/phot"'
TAP_URLS = [
    "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync",
    "https://tapvizier.u-strasbg.fr/TAPVizieR/tap/sync",
]
REQUIRED = ["id", "raj2000", "dej2000"]
INTERESTING = [
    "snrf444w", "whtf444w", "snrf115w", "snrf150w", "snrf277w",
    "whtf115w", "whtf150w", "whtf277w", "flagstar", "flagblend",
    "warnflag", "magautof444w", "fluxautof444w", "efluxautof444w",
    "tile", "segmentid",
]


def normal(name):
    return re.sub(r"[^a-z0-9]", "", str(name).lower())


def choose_columns(columns):
    lookup = {}
    for column in columns:
        key = normal(column)
        if key in lookup:
            raise ValueError(f"Ambiguous normalized CDS column {key}")
        lookup[key] = column
    missing = set(REQUIRED) - lookup.keys()
    if missing:
        raise ValueError(f"CDS catalog lacks columns: {sorted(missing)}")
    return lookup


def tap_query(adql, session=requests, timeout=120):
    errors = []
    for endpoint in TAP_URLS:
        try:
            r = session.get(
                endpoint,
                params={"REQUEST": "doQuery", "LANG": "ADQL",
                        "FORMAT": "csv", "QUERY": adql, "MAXREC": 1000},
                timeout=timeout,
            )
            r.raise_for_status()
            head = r.text[:300].lower()
            if "<votable" in head or "<html" in head or "query_status" in head:
                raise RuntimeError(f"Unexpected TAP response: {r.text[:400]}")
            table = pd.read_csv(StringIO(r.text), comment="#")
            if not len(table.columns):
                raise RuntimeError("TAP returned no CSV columns")
            return table, endpoint
        except Exception as exc:
            errors.append(f"{endpoint}: {type(exc).__name__}: {str(exc)[:300]}")
    raise RuntimeError("All CDS TAP endpoints failed: " + "; ".join(errors))


def validate_matches(roster, rows, max_sep_arcsec=0.5):
    cols = choose_columns(rows.columns)
    selected = rows.copy()
    selected["candidate_id"] = pd.to_numeric(selected[cols["id"]], errors="raise").astype("int64")
    if selected.candidate_id.duplicated().any():
        raise ValueError("CDS returned duplicate catalog IDs")
    merged = roster.merge(selected, on="candidate_id", how="left", indicator=True, validate="one_to_one")
    ra = pd.to_numeric(merged[cols["raj2000"]], errors="coerce")
    dec = pd.to_numeric(merged[cols["dej2000"]], errors="coerce")
    sep = np.hypot((ra - merged.ra_deg) * np.cos(np.deg2rad(merged.dec_deg)),
                   dec - merged.dec_deg) * 3600
    merged["cds_position_offset_arcsec"] = sep
    merged["cds_id_found"] = merged["_merge"].eq("both")
    merged["cds_coordinate_consistent"] = merged.cds_id_found & (sep <= max_sep_arcsec)
    for key in INTERESTING:
        if key in cols:
            merged["cds_" + key] = merged[cols[key]]
    output_cols = ["candidate_id", "ra_deg", "dec_deg",
                   "cds_id_found", "cds_position_offset_arcsec",
                   "cds_coordinate_consistent"] + [
                       "cds_" + key for key in INTERESTING if key in cols]
    return merged[output_cols].copy(), cols


def summarize(membership, source, columns):
    good = membership.cds_coordinate_consistent
    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "status": "CDS_CANDIDATE_FLAG_AUDIT_NOT_AREA_COMPLETENESS",
        "cds_table": TABLE,
        "tap_endpoint": source,
        "roster_entries": len(membership),
        "matched_by_id": int(membership.cds_id_found.sum()),
        "matched_by_id_and_position_within_0p5arcsec": int(good.sum()),
        "missing_or_mismatched_ids": membership.loc[~good, "candidate_id"].astype(int).tolist(),
        "available_cds_columns": [key for key in INTERESTING if key in columns],
        "not_available_cds_columns": [key for key in INTERESTING if key not in columns],
        "warning": [
            "CDS copy may differ from the frozen Dave selection catalog; compare source IDs and positions.",
            "flag_star is a graded overlap flag, not automatically a binary area mask.",
            "flag_blend is source-level blending, not a sky-area exclusion.",
            "SNR of already-detected catalog objects cannot define sky completeness at m_F444W=26.2.",
            "Blue-band non-detection selection is spectral-shape and noise dependent.",
            "Geometric multi-epoch coverage is not astrometric PM detectability.",
        ],
    }
    for col in ["cds_flagstar", "cds_flagblend", "cds_warnflag"]:
        if col in membership:
            report[col + "_distribution_valid_position_matches"] = {
                str(k): int(v) for k, v in
                membership.loc[good, col].astype("string").fillna("MISSING").value_counts().items()
            }
    for col in ["cds_snrf444w", "cds_snrf115w", "cds_snrf150w", "cds_snrf277w"]:
        if col in membership:
            val = pd.to_numeric(membership.loc[good, col], errors="coerce")
            report[col + "_n_finite"] = int(np.isfinite(val).sum())
            report[col + "_n_snr_gt_5"] = int((val > 5).sum())
            report[col + "_n_snr_lt_2"] = int((val < 2).sum())
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", default="pm_footprint/dave125_positions.csv")
    parser.add_argument("--membership", default="pm_footprint/results/wht60/dave125_wht_membership.csv")
    parser.add_argument("--out", default="pm_footprint/results/catalog_qc")
    args = parser.parse_args()
    output = Path(args.out)
    output.mkdir(parents=True, exist_ok=True)
    try:
        roster = pd.read_csv(args.candidates)
        if len(roster) != 125 or roster.candidate_id.nunique() != 125:
            raise ValueError("Frozen Dave roster is not 125 unique IDs")
        # A single small indexed ADQL query: only 125 known source IDs.
        ids = ",".join(map(str, sorted(roster.candidate_id.astype(int))))
        adql = f'SELECT * FROM {TABLE} WHERE "Id" IN ({ids})'
        rows, endpoint = tap_query(adql)
        membership, columns = validate_matches(roster, rows)
        wht = pd.read_csv(args.membership)
        if len(wht) != 125 or wht.candidate_id.nunique() != 125:
            raise ValueError("Existing WHT candidate membership roster is invalid")
        membership = membership.merge(
            wht[["candidate_id", "in_f444w_plus_red_365d", "in_any_nircam_365d"]],
            on="candidate_id", how="left", validate="one_to_one")
        report = summarize(membership, endpoint, columns)
        # Preserve provenance of the returned source; do not overwrite full MAST data.
        membership.to_csv(output / "dave125_cds_flags.csv", index=False)
        (output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report, indent=2))
    except Exception as exc:
        (output / "failure.json").write_text(json.dumps({
            "status": "NO_CATALOG_FLAG_CLAIM",
            "error": repr(exc),
            "utc": datetime.now(timezone.utc).isoformat(),
        }, indent=2) + "\n")
        raise


if __name__ == "__main__":
    main()
