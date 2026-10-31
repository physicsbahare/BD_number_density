"""Sparse direct-ID lookup into 6.6GB v13 BINTABLE (no full catalog load).

WARNING: Dave's 125 source IDs might belong to another catalog. All direct
ID-index matches MUST pass an independent RA/Dec tolerance, else stop
without treating the returned photometry as matched. A later separate
coordinate-index scan is required when identifiers differ.

Uses FITS big-endian fixed-column table offsets and read-only os.pread of
125 candidate rows. Calculates descriptive *scenario* selection flags for
the 125 matched rows only; does NOT derive population completeness or
selection effective area. Does not access 50GB images.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import re
import struct
from pathlib import Path
import numpy as np
from astropy.io import fits

WIDTH={"A":1,"L":1,"B":1,"I":2,"J":4,"K":8,"E":4,"D":8,"C":8,"M":16,"P":8,"Q":16}
SCALARS={"B":">B","I":">h","J":">i","K":">q","E":">f","D":">d"}

KEYS=["ID","RA","Dec","use_phot","use_aper","low_snr_flag",
      "bad_pixel_lw_flag","artifacts_near_badpixels_flag",
      "combined_artifact_flag","FLAG_AUTO_LW","FLAG_BAD_KRON_LW"]
for band in ("F115W","F150W","F277W","F444W"):
    KEYS.extend([f"f_{band}",f"e_{band}",f"new_e_{band}",
                 f"mag_{band}",f"FLAG_APER_032_{band}"])
KEYS+=["FLAG_APER_032_LW","FLUX_APER_032_LW","FLUXERR_APER_032_LW"]


def table_specs(path):
    with open(path,"rb") as f:
        primary=fits.Header.fromfile(f,endcard=True,padding=True)
        if int(primary.get("NAXIS",-1))!=0:
            raise ValueError("Expected empty primary HDU; no large image skip")
        h=fits.Header.fromfile(f,endcard=True,padding=True)
        offset=f.tell()
    if str(h.get("XTENSION","")).strip().upper()!="BINTABLE":
        raise ValueError("First extension must be BINTABLE")
    fields={}; pos=0
    for i in range(1,int(h["TFIELDS"])+1):
        name=str(h[f"TTYPE{i}"]).strip()
        form=str(h[f"TFORM{i}"]).strip().upper()
        m=re.fullmatch(r"(\d*)([A-Z])(?:\(\d+\))?",form)
        if not m or m.group(2) not in WIDTH:
            raise ValueError(f"Unsupported FITS type: {name}: {form}")
        n=int(m.group(1) or 1);letter=m.group(2)
        width=(n+7)//8 if letter=="X" else n*WIDTH[letter]
        fields[name]={"offset":pos,"type":letter,"count":n}
        pos+=width
    if pos!=int(h["NAXIS1"]):
        raise ValueError(f"Field width mismatch: {pos} != {h['NAXIS1']}")
    return offset,int(h["NAXIS1"]),int(h["NAXIS2"]),fields


def scalar(buf,spec):
    if spec is None or spec["count"]!=1 or spec["type"] not in SCALARS:
        return None
    return struct.unpack_from(SCALARS[spec["type"]],buf,spec["offset"])[0]


def sky_sep_arcsec(ra,dec,other_ra,other_dec):
    mean_dec=math.radians((dec+other_dec)/2)
    dra=(ra-other_ra)*math.cos(mean_dec)*3600
    ddec=(dec-other_dec)*3600
    return math.hypot(dra,ddec)


def selection_scenarios(rows):
    result={}
    bands=("F115W","F150W","F277W","F444W")
    for error_kind in ("e","new_e"):
        for mode in ("signed_non_detection_and_positive_ratio",
                     "absolute_non_detection_and_positive_ratio",
                     "signed_non_detection_and_flux_threshold",
                     "signed_non_detection_and_2sigma_upper_limit"):
            passed=0; eligible=0;indeterminate=0
            for row in rows:
                f=[float(row.get("f_"+b,float("nan"))) for b in bands]
                e=[float(row.get(error_kind+"_"+b,float("nan"))) for b in bands]
                if not all(math.isfinite(x) for x in f+e) or min(e)<=0:
                    indeterminate+=1;continue
                eligible+=1
                snr=[fi/ei for fi,ei in zip(f,e)]
                nondetect=(all(abs(snr[i])<2 for i in (0,1,2))
                           if mode.startswith("absolute") else
                           all(snr[i]<2 for i in (0,1,2)))
                if mode.endswith("positive_ratio"):
                    color=f[2]>0 and f[3]/f[2]>30
                elif mode.endswith("flux_threshold"):
                    color=f[2]<f[3]/30
                else:
                    color=f[2]+2*e[2]<f[3]/30
                if snr[3]>5 and nondetect and color:
                    passed+=1
            result[error_kind+"__"+mode]={
                "eligible_complete_numerical_rows":eligible,
                "passed_among_eligible":passed,"indeterminate_rows":indeterminate,
                "note":"Descriptive preselected-roster sensitivity ONLY, not recovered full-catalog selection, no usable area."
            }
    return result


def run(catalog,roster,output,tolerance_arcsec=0.5):
    catalog=Path(catalog);roster=Path(roster);output=Path(output)
    offset,row_bytes,nrows,fields=table_specs(catalog)
    for k in ["ID","RA","Dec"]+["f_"+b for b in ("F115W","F150W","F277W","F444W")]:
        if k not in fields or fields[k]["count"]!=1:
            raise ValueError("Missing/invalid expected field: "+k)
    with roster.open(newline="") as f:
        candidates=list(csv.DictReader(f))
    if len(candidates)!=125 or len(set(x["candidate_id"] for x in candidates))!=125:
        raise ValueError("Expected 125 unique frozen candidate IDs")
    checksum=hashlib.sha256(roster.read_bytes()).hexdigest()
    if checksum!="d43ce32e33b898f6bf57551affe41f209e0b455a4e4e63b4e99a8eea03236bbd":
        raise ValueError("Candidate roster checksum differs from frozen source")
    output.parent.mkdir(parents=True,exist_ok=True)
    rows=[];invalid=[]
    fd=os.open(catalog,os.O_RDONLY)
    try:
        for s in candidates:
            target=int(s["candidate_id"])
            if target<1 or target>nrows:
                invalid.append({"candidate_id":target,"problem":"outside table row range"})
                continue
            blob=os.pread(fd,row_bytes,offset+(target-1)*row_bytes)
            if len(blob)!=row_bytes:raise IOError("Truncated catalog table row")
            actual_id=scalar(blob,fields["ID"])
            actual_ra=scalar(blob,fields["RA"])
            actual_dec=scalar(blob,fields["Dec"])
            separation=sky_sep_arcsec(actual_ra,actual_dec,
                                      float(s["ra_deg"]),float(s["dec_deg"]))
            if actual_id!=target or not math.isfinite(separation) or separation>tolerance_arcsec:
                invalid.append({"candidate_id":target,"at_row_index_id":actual_id,
                                "angular_separation_arcsec":separation,
                                "problem":"ID-row or coordinate mismatch; catalogs may differ"})
                continue
            r={"candidate_id":target,"crossmatch_separation_arcsec":separation}
            for key in KEYS:
                r[key]=scalar(blob,fields.get(key))
            rows.append(r)
    finally:
        os.close(fd)
    report={
        "status":("ALL_125_DIRECT_IDS_AND_COORDINATES_VALIDATED"
                  if len(rows)==125 else "PARTIAL_OR_FAILED_ID_COORD_MATCH_DO_NOT_USE_AS_DAVE_PHOTOMETRY"),
        "catalog_basename":catalog.name,
        "catalog_rows":nrows,
        "catalog_row_bytes":row_bytes,
        "candidate_count":len(candidates),
        "matched_source_count":len(rows),
        "mismatched_source_count":len(invalid),
        "first_12_mismatch_diagnostics":invalid[:12],
        "max_match_separation_arcsec":max([r["crossmatch_separation_arcsec"] for r in rows],default=None),
        "field_status":{k:k in fields for k in KEYS},
        "roster_sha256":checksum,
        "selection_scenarios":selection_scenarios(rows) if len(rows)==125 else None,
        "warnings":[
            "Direct table offset ID-1 is ONLY a fast-path hypothesis and was explicitly checked against both actual IDs and sky coordinates.",
            "If fewer than 125 match, do not compare numerical selection scenarios: run an independent sky crossmatch/scan.",
            "Both initial e_* and notebook empirically inflated new_e_* are present; do not silently choose one as Dave's definition.",
            "Selection scenario counts on preselected sources cannot measure selection function/completeness.",
            "No 6.6GB table scan; at most 125 individual ~2KB rows were read."
        ]
    }
    output.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    csv_file=output.parent/(output.stem+"_candidate_rows.csv")
    if rows:
        with csv_file.open("w",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]))
            writer.writeheader();writer.writerows(rows)
    print("Matched",len(rows),"of 125 candidate ID+position pairs")
    print("Wrote",output,"and",csv_file if rows else "no candidate CSV")
    if invalid: print("Warning: ID/position mismatch; do NOT treat row offsets as Dave's source photometry.")
    return report


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--catalog",required=True)
    p.add_argument("--roster",default="pm_footprint/dave125_positions.csv")
    p.add_argument("--out",default="pm_footprint/results/local_v13_catalog_idcheck.json")
    p.add_argument("--max-separation-arcsec",type=float,default=0.5)
    a=p.parse_args()
    run(a.catalog,a.roster,a.out,a.max_separation_arcsec)


if __name__=="__main__":main()
