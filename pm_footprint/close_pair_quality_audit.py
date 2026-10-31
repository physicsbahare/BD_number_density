"""Source-level 0.2--1.0 arcsec proximity and Kron-flag audit of Dave's 125.

Proximity components are not proven duplicates or contaminants. This uses
only the frozen, source-verified compact v13 fixture; no WHT/FITS downloads.
"""
import argparse
import csv
import json
import math
from pathlib import Path

from .brightness_stratified_density import fisher_two_sided

RED_AREA = 0.31462445722188803
TOTAL_AREA = 0.528823684998808
RADII = (0.2, 0.3, 0.4, 0.5, 0.6, 1.0)


def separation_arcsec(a, b):
    """Small-angle sky separation (sufficient at <=1 arcsec)."""
    dx = (a["ra"]-b["ra"])*3600*math.cos(math.radians((a["dec"]+b["dec"])/2))
    dy = (a["dec"]-b["dec"])*3600
    return math.hypot(dx, dy)


def parse_rows(rows):
    data = []
    for row in rows:
        obj = {"id": int(row["candidate_id"]),
               "ra": float(row["ra_deg"]), "dec": float(row["dec_deg"]),
               "red": int(row["red365"]), "mag": float(row["mag_F444W"]),
               "bad_kron": int(row["FLAG_BAD_KRON_LW"])}
        if obj["red"] not in (0, 1) or obj["bad_kron"] not in (0, 1):
            raise ValueError("Invalid binary coverage/Kron flag")
        if not all(math.isfinite(obj[k]) for k in ("ra", "dec", "mag")):
            raise ValueError("Invalid sky coordinate/magnitude")
        data.append(obj)
    if len(data) != 125 or len({x["id"] for x in data}) != 125:
        raise ValueError("Expected frozen 125 unique source IDs")
    return data


def connected_components(data, radius):
    if not 0 < radius <= 2:
        raise ValueError("Radius must be in (0,2] arcsec")
    n = len(data)
    parent = list(range(n))
    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    edges = []
    for i in range(n):
        for j in range(i+1, n):
            d = separation_arcsec(data[i], data[j])
            if d <= radius:
                parent[root(i)] = root(j)
                edges.append((i, j, d))
    groups = {}
    for i in range(n):
        groups.setdefault(root(i), []).append(i)
    return list(groups.values()), edges


def density_ratio(n_red, n_nonred):
    return ((n_red/RED_AREA)/(n_nonred/(TOTAL_AREA-RED_AREA))
            if n_red and n_nonred else None)


def audit(rows):
    data = parse_rows(rows)
    sweeps = []
    main_groups = None
    main_edges = None
    for radius in RADII:
        groups, edges = connected_components(data, radius)
        members = {i for g in groups if len(g) > 1 for i in g}
        subset = [data[i] for i in sorted(members)]
        sweeps.append({"radius_arcsec": radius,
                       "n_components": len(groups),
                       "n_multi_member_components": sum(len(g)>1 for g in groups),
                       "n_close_members": len(members),
                       "red_close_members": sum(x["red"] for x in subset),
                       "nonred_close_members": sum(1-x["red"] for x in subset),
                       "bad_kron_close_members": sum(x["bad_kron"] for x in subset)})
        if radius == 0.5:
            main_groups, main_edges = groups, edges
    members = {i for g in main_groups if len(g)>1 for i in g}
    close = [x for i,x in enumerate(data) if i in members]
    isolated = [x for i,x in enumerate(data) if i not in members]
    bad_close = sum(x["bad_kron"] for x in close)
    bad_isolated = sum(x["bad_kron"] for x in isolated)
    red_close = sum(x["red"] for x in close)
    red_isolated = sum(x["red"] for x in isolated)
    details = []
    for group in main_groups:
        if len(group) <= 1:
            continue
        ids = sorted(data[i]["id"] for i in group)
        if len({data[i]["red"] for i in group}) != 1:
            raise ValueError("Close component has mixed red epoch geometry")
        ds = [d for i,j,d in main_edges if i in group and j in group]
        details.append({"candidate_ids": ids,
                        "n": len(ids),
                        "red365": data[group[0]]["red"],
                        "bad_kron_count": sum(data[i]["bad_kron"] for i in group),
                        "mag_min": min(data[i]["mag"] for i in group),
                        "mag_max": max(data[i]["mag"] for i in group),
                        "link_separations_arcsec": [round(d, 6) for d in sorted(ds)]})
    details.sort(key=lambda x:x["candidate_ids"][0])
    bright_all = [x for x in data if x["mag"]<25]
    bright_isolated = [x for x in isolated if x["mag"]<25]
    def counts(sub):
        nr = sum(x["red"] for x in sub)
        return {"red":nr,"nonred":len(sub)-nr,
                "red_nonred_area_normalized_density_ratio":density_ratio(nr,len(sub)-nr)}
    return {
        "status":"EXPLORATORY_SOURCE_PROXIMITY_KRON_ASSOCIATION_NOT_DUPLICATE_CLASSIFICATION",
        "source":"frozen 125-row independently ID-verified v13 photometry fixture",
        "source_original_125_rows_sha256":"93400baee3daccfc5020b01f93d05e9b65ef238afa62b5b691f4aa90728e2a40",
        "radius_sweep":sweeps,
        "close_0p5arcsec":{
            "n_multi_components":len(details),
            "n_close_members":len(close),
            "n_isolated":len(isolated),
            "bad_kron_close":bad_close,
            "bad_kron_isolated":bad_isolated,
            "bad_kron_odds_ratio_close_vs_isolated":
                (bad_close*(len(isolated)-bad_isolated)) /
                ((len(close)-bad_close)*bad_isolated),
            "bad_kron_fisher_two_sided_exploratory_p":
                fisher_two_sided([[bad_close,len(close)-bad_close],
                                  [bad_isolated,len(isolated)-bad_isolated]]),
            "red_close":red_close,"nonred_close":len(close)-red_close,
            "red_isolated":red_isolated,
            "nonred_isolated":len(isolated)-red_isolated,
            "red_nonred_close_membership_fisher_exploratory_p":
                fisher_two_sided([[red_close,red_isolated],
                                  [len(close)-red_close,len(isolated)-red_isolated]]),
            "groups":details},
        "density_sensitivity":{
            "all_sources":counts(data),
            "isolated_only":counts(isolated),
            "bright_m444_lt25_all":counts(bright_all),
            "bright_m444_lt25_isolated_only":counts(bright_isolated)},
        "limitations":[
            "<=0.5 arcsec is an exploratory proximity threshold, not a duplicate or blend criterion.",
            "Kron flags may reflect photometry failures, not necessarily spurious objects; no automatic exclusion justified.",
            "Fisher tests ignore spatial dependence, threshold selection and multiple comparisons; p-values are descriptive.",
            "Removing all close-group members is an illustrative sensitivity extreme, not a recommended catalog cut.",
            "No source segmentation, F444W PSF fitting, F356W SNR or measured PM completeness.",
            "WHT areas are public 60mas gridded WHT>0 geometric areas, not dropout selection-effective area."
        ]}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--fixture",default="pm_footprint/results/v13_epoch_candidate_photometry_fixture_20261031.csv")
    p.add_argument("--out",default="pm_footprint/results/close_pair_quality_audit_20261101.json")
    a=p.parse_args()
    with open(a.fixture,newline="") as f:
        result=audit(csv.DictReader(f))
    path=Path(a.out)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"close":result["close_0p5arcsec"]["n_close_members"],
                      "groups":result["close_0p5arcsec"]["n_multi_components"],
                      "bad_kron_fisher":result["close_0p5arcsec"]["bad_kron_fisher_two_sided_exploratory_p"]}))


if __name__=="__main__":
    main()
