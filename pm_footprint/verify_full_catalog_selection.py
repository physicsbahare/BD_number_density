"""Offline scientific QA of user's full v13 catalog selection-count scan.

Input is the compact JSON from v13_full_catalog_paper_stages.py; no images,
catalog download, FITS memory map, MAST calls, or extra scanning required.
Verifies the manuscript's 2997→335→125 and frozen-ID identity, compares
the essential use_phot gate and empirically corrected/new_e error convention,
and emits a concise independently inspectable evidence JSON.
"""
import argparse
import hashlib
import json
from pathlib import Path

EXPECTED=[2997,335,125]
SCENARIOS=("new_e_all","new_e_use_phot","e_all","e_use_phot")


def summarize(data):
    if data.get("status")!="FULL_LOCAL_V13_CATALOG_PAPER_SELECTION_STAGE_AUDIT":
        raise ValueError("Unexpected input provenance/status")
    if int(data["catalog_rows"])!=3280023 or int(data["FITS_row_bytes"])!=2003:
        raise ValueError("Not the verified Zhen-Kai v13 input catalog")
    if data["paper_reported_stage_counts"]!=EXPECTED:
        raise ValueError("Paper's 3 stage counts changed")
    sc=data["scenarios"]
    if set(sc)!=set(SCENARIOS):
        raise ValueError("Missing or unexpected scenarios")
    for key in SCENARIOS:
        r=sc[key]
        stages=r["stage_counts"]
        if len(stages)!=3 or not all(isinstance(x,int) for x in stages):
            raise ValueError(f"Malformed stage counts: {key}")
        if not (0<=stages[2]<=stages[1]<=stages[0]<=r["complete_numeric_rows"]<=data["catalog_rows"]):
            raise ValueError(f"Non-monotonic selection stages: {key}")
        if stages[2]!=r["final_ID_count"]:
            raise ValueError(f"Duplicate final IDs or inconsistent counting: {key}")
        if (r["frozen_125_IDs_recovered"]+r["not_in_frozen_125_count"]!=stages[2]
                or r["frozen_125_IDs_recovered"]+r["missing_frozen_125_count"]!=125):
            raise ValueError(f"Frozen roster vs final count inconsistent: {key}")
    for prefix in ("new_e","e"):
        for key in ("complete_numeric_rows",):
            if sc[prefix+"_use_phot"][key]>sc[prefix+"_all"][key]:
                raise ValueError("use_phot gate increased valid source rows")
        for i in range(3):
            if sc[prefix+"_use_phot"]["stage_counts"][i]>sc[prefix+"_all"]["stage_counts"][i]:
                raise ValueError("use_phot gate increased selected sources")
    valid=sc["new_e_use_phot"]
    reproduced=(valid["stage_counts"]==EXPECTED
                and valid["final_ID_count"]==125
                and valid["frozen_125_IDs_recovered"]==125
                and valid["not_in_frozen_125_count"]==0
                and valid["missing_frozen_125_count"]==0)
    excluded=[
        a-b for a,b in zip(sc["new_e_all"]["stage_counts"],valid["stage_counts"])
    ]
    return {
        "status":"FULL_3P28M_V13_SELECTION_REPRODUCED" if reproduced else "SELECTION_NOT_REPRODUCED",
        "reported_vs_reproduced_stage_counts":{
            "manuscript":EXPECTED,"new_e_with_use_phot":valid["stage_counts"]},
        "corrected_error_with_gate_matches_all_125_IDs":reproduced,
        "catalog_rows":data["catalog_rows"],
        "source_uploaded_full_scan_json_sha256":data.get("user_uploaded_original_json_sha256"),
        "complete_numeric_rows_new_e_without_gate":sc["new_e_all"]["complete_numeric_rows"],
        "complete_numeric_rows_new_e_with_gate":valid["complete_numeric_rows"],
        "numeric_rows_excluded_by_use_phot_gate":(
            sc["new_e_all"]["complete_numeric_rows"]-valid["complete_numeric_rows"]),
        "stage_counts_by_scenario":{k:sc[k]["stage_counts"] for k in SCENARIOS},
        "candidate_ID_recovery_by_scenario":{
            k:{"in_frozen":sc[k]["frozen_125_IDs_recovered"],
               "unexpected":sc[k]["not_in_frozen_125_count"],
               "missing_frozen":sc[k]["missing_frozen_125_count"]}
            for k in SCENARIOS
        },
        "additional_candidates_excluded_by_use_phot_gate_per_stage":excluded,
        "no_use_phot_gate_final_pass_count":sc["new_e_all"]["stage_counts"][-1],
        "fraction_of_ungated_corrected_final_retained_after_gate":(
            valid["stage_counts"][-1]/sc["new_e_all"]["stage_counts"][-1]),
        "conditional_step_retention":{
            "ratio_given_SNR":valid["stage_counts"][1]/valid["stage_counts"][0],
            "magnitude_given_ratio":valid["stage_counts"][2]/valid["stage_counts"][1],
            "all_cuts_given_SNR":valid["stage_counts"][2]/valid["stage_counts"][0]},
        "science_limits":[
            "Full *catalog* selection and final source IDs are reproducible, not the original image reduction.",
            "use_phot removes sources by upstream detection/quality logic; the 1039 extra unflagged candidates are NOT demonstrated brown dwarfs or necessarily bad data.",
            "No completeness/injection-recovery probability or spatial selection-effective area is measured.",
            "Do not divide or multiply predicted Galactic brown-dwarf counts by the geometric PM footprint fraction to infer the actual number of measurable proper motions.",
            "A separate F356W cross-filter astrometric sensitivity assessment and a matched four-band depth/spatial mask remain required."
        ]
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--input",type=Path,required=True)
    parser.add_argument("--out",type=Path,default=Path("pm_footprint/results/v13_full_catalog_selection_verified_20261031.json"))
    args=parser.parse_args()
    data=json.loads(args.input.read_text())
    result=summarize(data)
    result["input_json_sha256"]=hashlib.sha256(args.input.read_bytes()).hexdigest()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps({k:result[k] for k in (
        "status","reported_vs_reproduced_stage_counts",
        "stage_counts_by_scenario",
        "additional_candidates_excluded_by_use_phot_gate_per_stage"
    )},indent=2))
    if result["status"]!="FULL_3P28M_V13_SELECTION_REPRODUCED":
        raise SystemExit("Selection mismatch: preserve diagnostics for review")


if __name__=="__main__":
    main()
