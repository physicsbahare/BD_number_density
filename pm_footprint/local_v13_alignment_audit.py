"""Check v13 four-band pixel alignment WITHOUT trusting inconsistent blue WCS.

Zhen-Kai notebook performs forced photometry on the common F444W
segmentation pixel coordinates. Uploaded blue SCI/ERR FITS headers currently
map the same sky points ~67565 pixels from F444W WHT! This diagnostic
samples native ERR at BOTH sky-WCS and *shared F444W pixel* positions.

Shared-index results are a CONDITIONAL co-registration hypothesis, NOT
proof that the images are aligned. Never label WCS-out-of-bounds as missing
observations. Read-only, 3x3 pixels/125 positions, no full image load.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import numpy as np
from .local_v13_photometry_probe import primary_header, read_native_patch, summarize_patch

FILTERS=("F115W","F150W","F277W","F444W")


def audit(root,roster,out):
    root,roster,out=Path(root),Path(roster),Path(out)
    _,_,wht,wcs=primary_header(root/"F444W_psf_matched_F444W_WHT_v13.fits")
    if wcs is None: raise ValueError("Reference F444W WHT celestial WCS missing")
    with roster.open(newline="") as f: sources=list(csv.DictReader(f))
    if len(sources)!=125 or len(set(x["candidate_id"] for x in sources))!=125:
        raise ValueError("Expected 125 unique frozen sources")
    ra=np.array([float(x["ra_deg"]) for x in sources])
    dec=np.array([float(x["dec_deg"]) for x in sources])
    xr,yr=wcs.all_world2pix(ra,dec,0)
    rows=[{"candidate_id":x["candidate_id"]} for x in sources]
    report={"status":"HEADER_WCS_DISCREPANCY_AND_SHARED_PIXEL_HYPOTHESIS_NOT_DEPTH",
            "reference":"F444W_WHT WCS source-center coordinates",
            "per_filter":{}}
    for filt in FILTERS:
        p=root/(f"{filt}_psf_matched_F444W_ERR_v13.fits")
        h,off,info,w=primary_header(p)
        if w is None:raise ValueError(f"{filt} has no celestial WCS")
        x,y=w.all_world2pix(ra,dec,0)
        mismatch=np.hypot(x-xr,y-yr)
        fd=os.open(p,os.O_RDONLY)
        vals_ref=[];vals_header=[]
        try:
            for k in range(len(sources)):
                a,_=read_native_patch(fd,off,h,xr[k],yr[k])
                b,_=read_native_patch(fd,off,h,x[k],y[k])
                vals_ref.append(summarize_patch(a)["center_err"])
                vals_header.append(summarize_patch(b)["center_err"])
        finally:
            os.close(fd)
        positive=lambda v: sum(z is not None and z>0 for z in v)
        finite=lambda v: sum(z is not None and np.isfinite(z) for z in v)
        report["per_filter"][filt]={
            "header_to_reference_max_offset_pixels":float(mismatch.max()),
            "header_to_reference_median_offset_pixels":float(np.median(mismatch)),
            "same_index_valid_centers":finite(vals_ref),
            "same_index_positive_err_centers":positive(vals_ref),
            "header_wcs_valid_centers":finite(vals_header),
            "header_wcs_positive_err_centers":positive(vals_header),
            "pixel_err_units":info["keywords"]["BUNIT"],
            "interpretation":"Same-index test assumes image arrays co-registered despite FITS WCS; must check bright-source SCI morphology before trusting."
        }
        for row,v,hv in zip(rows,vals_ref,vals_header):
            row[f"{filt}_err_same_F444W_pixel"]=v
            row[f"{filt}_err_header_WCS_pixel"]=hv
        print(filt,report["per_filter"][filt],flush=True)
    report["limitations"]=[
        "WCS discrepancy may be a metadata problem or actual science-data misregistration; presence of nonzero ERR at same pixel does NOT demonstrate alignment.",
        "No PSF/aperture completeness, flag, image-source detection or calibrated depth inferred.",
        "F444W WHT-positive source points and blue science image source positions require an astrometric/morphological check.",
        "ERR 3x3 center is not an aperture-level or correlation-corrected noise estimate."
    ]
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    with (out.parent/(out.stem+"_centers.csv")).open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)
    return report


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--root",required=True)
    p.add_argument("--roster",default="pm_footprint/dave125_positions.csv")
    p.add_argument("--out",default="pm_footprint/results/local_v13_alignment_audit.json")
    a=p.parse_args()
    audit(a.root,a.roster,a.out)


if __name__=="__main__":main()
