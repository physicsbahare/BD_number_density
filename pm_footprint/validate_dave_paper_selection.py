"""Reproduce the Dave et al. 2026-10-29 draft Section 2.2 cuts on frozen v13 candidates.

Paper recipe: signed SNR_F444W>5; signed SNR_F115W,F150W,F277W<2;
f444/f277_eff>30 with f277_eff=f277 if f277>0 else 2*e277;
m_F444W<26.2. Inputs are matched Zhen-Kai v13 catalog rows (nJy).
Does NOT independently discover the 125 sources or measure completeness.
"""
import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
import numpy as np

BANDS = ("F115W", "F150W", "F277W", "F444W")


def evaluate(rows, prefix="new_e_"):
    if not rows:
        raise ValueError("Empty catalog rows")
    f = {b: np.array([float(r["f_"+b]) for r in rows]) for b in BANDS}
    e = {b: np.array([float(r[prefix+b]) for r in rows]) for b in BANDS}
    m444 = np.array([float(r["mag_F444W"]) for r in rows])
    if (not all(np.isfinite(x).all() for x in [*f.values(),*e.values(),m444])
            or any(np.any(x<=0) for x in e.values())):
        raise ValueError("Nonfinite flux or nonpositive ERR")
    sn = {b: f[b]/e[b] for b in BANDS}
    feff = np.where(f["F277W"]>0, f["F277W"], 2*e["F277W"])
    ratio = f["F444W"]/feff
    checks = {
        "F444W_SNR_gt5": sn["F444W"]>5,
        "F115W_SNR_lt2": sn["F115W"]<2,
        "F150W_SNR_lt2": sn["F150W"]<2,
        "F277W_SNR_lt2": sn["F277W"]<2,
        "effective_ratio_gt30": ratio>30,
        "F444W_AB_mag_lt26p2": m444<26.2,
    }
    passed = np.logical_and.reduce(list(checks.values()))
    return {
        "counts_by_cut": {k:int(v.sum()) for k,v in checks.items()},
        "all_cuts_passed":int(passed.sum()),
        "fails_candidate_ids":[int(r["candidate_id"]) for r,ok in zip(rows,passed) if not ok],
        "negative_or_zero_F277W":int(np.sum(f["F277W"]<=0)),
        "min_F444W_SNR":float(np.min(sn["F444W"])),
        "max_F277W_SNR":float(np.max(sn["F277W"])),
        "min_effective_ratio":float(np.min(ratio)),
        "max_F444W_AB_mag":float(np.max(m444)),
        "negative_blue_SNR_at_most_minus2":{
            b:[int(r["candidate_id"]) for r,s in zip(rows,sn[b]) if s<=-2]
            for b in ("F115W","F150W","F277W")
        },
    }


def analyze(rows):
    ids=[int(r["candidate_id"]) for r in rows]
    if len(ids)!=125 or len(set(ids))!=125:
        raise ValueError("Must use the frozen 125 unique candidate roster")
    if max(float(r["crossmatch_separation_arcsec"]) for r in rows)>0.5:
        raise ValueError("Candidate ID matches have unacceptably large sky offsets")
    corrected=evaluate(rows,"new_e_")
    raw=evaluate(rows,"e_")
    return {
        "status":"PAPER_SELECTION_REPRODUCED_ON_125_PRESELECTED_SOURCES_ONLY",
        "paper_section":"Kumar et al. draft 2026-10-29, Section 2.2",
        "candidate_count":len(rows),
        "corrected_new_e":corrected,
        "uncorrected_e":raw,
        "corrected_matches_all_frozen_125":corrected["all_cuts_passed"]==125,
        "positive_F277W_flux":125-corrected["negative_or_zero_F277W"],
        "use_phot_counts":dict(Counter(str(r["use_phot"]) for r in rows)),
        "use_aper_counts":dict(Counter(str(r["use_aper"]) for r in rows)),
        "FLAG_AUTO_LW_counts":dict(Counter(str(r["FLAG_AUTO_LW"]) for r in rows)),
        "FLAG_BAD_KRON_LW_counts":dict(Counter(str(r["FLAG_BAD_KRON_LW"]) for r in rows)),
        "caveats":[
            "This 125-row validation does not recreate the paper's 2997->335->125 FULL catalog selection.",
            "Dave's exact discovery script is not yet available; agreement establishes a strong candidate recipe consistency check.",
            "S/N<2 in the paper is signed (not absolute), so large negative flux can pass; such cases merit artifact assessment.",
            "The 2sigma F277W denominator applies ONLY when F277W flux<=0, NOT to positive F277W values.",
            "Corrected errors are consistent with the manuscript and notebook, but photometric selection area/completeness and F356W PM sensitivity are still unmeasured.",
            "WCS-header mismatch in blue images is distinct from astrophysical nondetection; test shared image-array registration before using local depth maps."
        ]
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--rows",type=Path,required=True)
    p.add_argument("--out",type=Path,default=Path("pm_footprint/results/dave125_paper_selection_local.json"))
    a=p.parse_args()
    with a.rows.open(newline="") as f: rows=list(csv.DictReader(f))
    report=analyze(rows)
    report["input_rows_sha256"]=hashlib.sha256(a.rows.read_bytes()).hexdigest()
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(report,indent=2)+"\n")
    print("Corrected:",report["corrected_new_e"]["all_cuts_passed"],"/125; raw:",
          report["uncorrected_e"]["all_cuts_passed"],"/125")


if __name__=="__main__":
    main()
