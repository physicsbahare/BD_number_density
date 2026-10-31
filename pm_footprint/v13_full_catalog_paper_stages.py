"""One-pass bounded-memory full v13 BINTABLE paper-selection audit (6.57GB local).

Reproduce 2026-10-29 manuscript Sec 2.2 stage counts 2997 / 335 / 125
under corrected new_e_* and, separately, raw e_* error scenarios, with
and without use_phot==1. Does not require loading the full catalog,
SCI/ERR images, or remote FITS. Preserves original files unchanged.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path

import numpy as np

from .local_v13_catalog_idcheck import table_specs

BANDS=("F115W","F150W","F277W","F444W")
NAMES=["ID","RA","Dec","use_phot","mag_F444W"]+[
    key+band for key in ("f_","e_","new_e_") for band in BANDS]


def extract_columns(buffer, offset_bytes, nrows, row_bytes, specs):
    columns={}
    for name in NAMES:
        if name not in specs:
            raise KeyError("Column missing from original v13 table: "+name)
        spec=specs[name]
        if spec["count"]!=1 or spec["type"] not in ("D","E","K","J","I","B"):
            raise ValueError("Unsupported selected FITS field: "+name)
        fmt={"D":">f8","E":">f4","K":">i8","J":">i4","I":">i2","B":"u1"}[spec["type"]]
        columns[name]=np.ndarray(shape=(nrows,),dtype=np.dtype(fmt),buffer=buffer,
                                 offset=spec["offset"],strides=(row_bytes,))
    return columns


def selected_stages(c,prefix,require_use_phot):
    f={b:c["f_"+b].astype("float64",copy=False) for b in BANDS}
    e={b:c[prefix+b].astype("float64",copy=False) for b in BANDS}
    m=c["mag_F444W"].astype("float64",copy=False)
    valid=np.isfinite(m)
    for b in BANDS:
        valid &= np.isfinite(f[b]) & np.isfinite(e[b]) & (e[b]>0)
    if require_use_phot:
        valid &= c["use_phot"]==1
    with np.errstate(divide="ignore",invalid="ignore"):
        stage_snr=valid & (f["F444W"]/e["F444W"]>5)
        for b in BANDS[:3]:
            stage_snr &= (f[b]/e[b]<2)
        effective_f277=np.where(f["F277W"]>0,f["F277W"],2*e["F277W"])
        stage_ratio=stage_snr & (f["F444W"]/effective_f277>30)
    stage_mag=stage_ratio & (m<26.2)
    return {
        "snr":stage_snr,"ratio":stage_ratio,"mag":stage_mag,
        "count_valid":int(valid.sum())
    }


def scan(catalog,roster,chunk_rows=20000):
    if chunk_rows<100 or chunk_rows>100000:
        raise ValueError("chunk_rows must be between 100 and 100000")
    catalog,roster=Path(catalog),Path(roster)
    offset,row_bytes,total,specs=table_specs(catalog)
    with roster.open(newline="") as stream:
        original=list(csv.DictReader(stream))
    frozen={int(r["candidate_id"]) for r in original}
    if len(original)!=125 or len(frozen)!=125:
        raise ValueError("Original frozen roster must have 125 unique IDs")
    if hashlib.sha256(roster.read_bytes()).hexdigest()!=(
            "d43ce32e33b898f6bf57551affe41f209e0b455a4e4e63b4e99a8eea03236bbd"):
        raise ValueError("Frozen candidate roster changed")
    keys=[f"{kind}_{phot}" for kind in ("new_e","e") for phot in ("all","use_phot")]
    accum={k:{"snr":0,"ratio":0,"mag":0,"valid":0,
              "last_stage_ids":set()} for k in keys}
    fd=os.open(catalog,os.O_RDONLY)
    try:
        for start in range(0,total,chunk_rows):
            count=min(chunk_rows,total-start)
            raw=os.pread(fd,count*row_bytes,offset+start*row_bytes)
            if len(raw)!=count*row_bytes:
                raise IOError("FITS catalog read truncated")
            cols=extract_columns(raw,offset,count,row_bytes,specs)
            for prefix in ("new_e_","e_"):
                for required,name in ((False,"all"),(True,"use_phot")):
                    key=prefix.rstrip("_")+"_"+name
                    stages=selected_stages(cols,prefix,required)
                    target=accum[key]
                    target["valid"]+=stages["count_valid"]
                    for stage in ("snr","ratio","mag"):
                        target[stage]+=int(stages[stage].sum())
                    target["last_stage_ids"].update(
                        map(int,cols["ID"][stages["mag"]]))
            if (start//chunk_rows)%25==0 or start+count>=total:
                print(f"read {start+count:,}/{total:,} v13 rows",flush=True)
    finally:
        os.close(fd)
    out={}
    for k,d in accum.items():
        final=d["last_stage_ids"]
        out[k]={
            "complete_numeric_rows":d["valid"],
            "stage_counts":[d["snr"],d["ratio"],d["mag"]],
            "final_ID_count":len(final),
            "frozen_125_IDs_recovered":len(final & frozen),
            "not_in_frozen_125_count":len(final-frozen),
            "missing_frozen_125_count":len(frozen-final),
            "unexpected_final_IDs_first_30":sorted(final-frozen)[:30],
            "missing_frozen_IDs_first_30":sorted(frozen-final)[:30]
        }
    return {
        "status":"FULL_LOCAL_V13_CATALOG_PAPER_SELECTION_STAGE_AUDIT",
        "input_catalog_basename":catalog.name,
        "catalog_bytes":catalog.stat().st_size,"catalog_rows":total,
        "FITS_row_bytes":row_bytes,
        "paper_reported_stage_counts":[2997,335,125],
        "scenarios":out,
        "caveats":[
            "Counts apply paper's signed SNR<2 and F277W positive-flux / nonpositive-2sigma split; compare paper's original software before inferring cause of disagreement.",
            "Use_phot was created from F444W 0.32-arcsec SNR and artifact flags before new_e; counts with and without it are both retained.",
            "Source-level selection through 3.28M catalogue is not photometric completeness or spatial F444W=26.2 recovery area.",
            "This scan uses original source rows only; masks, blending and cross-filter PM detectability remain independent tasks.",
            "One sequential 6.57GB read on user's local drive; no 50GB image input or MAST download."
        ]
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--catalog",required=True)
    p.add_argument("--roster",default="pm_footprint/dave125_positions.csv")
    p.add_argument("--chunk-rows",type=int,default=20000)
    p.add_argument("--out",type=Path,default=Path("pm_footprint/results/v13_full_catalog_paper_stages.json"))
    a=p.parse_args()
    r=scan(a.catalog,a.roster,a.chunk_rows)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(r,indent=2)+"\n")
    print(json.dumps({k:v["stage_counts"] for k,v in r["scenarios"].items()},indent=2))


if __name__=="__main__":
    main()
