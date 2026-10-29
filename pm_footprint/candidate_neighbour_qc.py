"""Audit close candidate positions without assuming they are duplicate sources.

Connected-component grouping is only a sensitivity diagnostic for number counts;
no source is automatically removed from the published 125-candidate roster.
"""
from __future__ import annotations
import json
import math
from pathlib import Path

import pandas as pd


def cluster(df, threshold_arcsec):
    rows = df.reset_index(drop=True)
    n = len(rows)
    parent = list(range(n))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def merge(i, j):
        a, b = root(i), root(j)
        if a != b:
            parent[b] = a

    pairs = []
    for i in range(n):
        a = rows.iloc[i]
        for j in range(i + 1, n):
            b = rows.iloc[j]
            mean_dec = math.radians((float(a.dec_deg) + float(b.dec_deg)) / 2)
            dx = (float(a.ra_deg) - float(b.ra_deg)) * 3600 * math.cos(mean_dec)
            dy = (float(a.dec_deg) - float(b.dec_deg)) * 3600
            distance = math.hypot(dx, dy)
            if distance < threshold_arcsec:
                merge(i, j)
                pairs.append({
                    "id1": int(a.candidate_id),
                    "id2": int(b.candidate_id),
                    "separation_arcsec": round(distance, 4),
                })
    groups = {}
    for i in range(n):
        groups.setdefault(root(i), []).append(int(rows.iloc[i].candidate_id))
    return {
        "raw_catalog_entries": n,
        "positional_clusters": len(groups),
        "close_pairs": pairs,
        "multi_entry_clusters": [x for x in groups.values() if len(x) > 1],
        "warning": "Angular proximity is not evidence that sources are duplicates; inspect deblending and imaging.",
    }


def main():
    root = Path("pm_footprint/results/exact")
    df = pd.read_csv(root / "dave125_exact_membership.csv")
    result = {"status": "POSITIONAL_QC_ONLY_NOT_DEDUPLICATED",
              "thresholds_arcsec": [0.5, 1.0], "samples": {}}
    for key in ["all_125"] + [c for c in df.columns if c.startswith("in_") and c != "in_archive_f444w"]:
        subset = df if key == "all_125" else df[df[key].astype(bool)]
        result["samples"][key] = {
            str(thresh): cluster(subset, thresh) for thresh in [0.5, 1.0]
        }
    (root / "neighbour_qc.json").write_text(json.dumps(result, indent=2) + "\n")
    print("Near-neighbour QC:", {
        key: (v["0.5"]["raw_catalog_entries"], v["0.5"]["positional_clusters"])
        for key, v in result["samples"].items()
    })


if __name__ == "__main__":
    main()
