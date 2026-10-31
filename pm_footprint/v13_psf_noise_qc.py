"""Audit Zhen-Kai v13 small PSF, EE table and notebook noise correction.

Requires only ~167kB PSF FITS, 5kB EE table, ~355kB original notebook.
No notebook execution, no 6.6GB catalog or 50GB mosaic read. This is a
calibration-consistency diagnostic, NOT a spatial completeness measurement.
"""
import argparse
import ast
import hashlib
import json
import re
from pathlib import Path
import numpy as np

BANDS=("F115W","F150W","F277W","F444W")
APERTURES=(0.32,0.48,0.70,1.00,1.40)
RADII=(0.16,0.24,0.35,0.5,0.7,1.0)


def read_small_psf(path):
    with open(path,"rb") as f:
        h={}
        while True:
            block=f.read(2880)
            if len(block)!=2880:raise ValueError("Truncated FITS header")
            stop=False
            for k in range(0,2880,80):
                line=block[k:k+80].decode("ascii")
                key=line[:8].strip()
                if key=="END":stop=True;break
                if line[8:10]=="= ":
                    value=line[10:].split("/",1)[0].strip().strip("'")
                    try:value=float(value) if "." in value else int(value)
                    except ValueError:pass
                    h[key]=value
            if stop:break
        if h.get("BITPIX")!=-32 or h.get("NAXIS")!=2:
            raise ValueError("Expected primary float32 2D PSF")
        nx,ny=int(h["NAXIS1"]),int(h["NAXIS2"])
        if nx*ny>1000000:raise ValueError("Refuse large image")
        blob=f.read(4*nx*ny)
        if len(blob)!=4*nx*ny:raise ValueError("Truncated PSF")
    return np.frombuffer(blob,dtype=">f4").astype(float).reshape(ny,nx)


def enclosed_fraction(psf,radius_arcsec,pixel_arcsec=.03,subsamples=8):
    if any(v%2!=1 for v in psf.shape):raise ValueError("Expected odd PSF shape")
    y,x=np.indices(psf.shape,dtype=float)
    cy,cx=(np.array(psf.shape)-1)/2
    coverage=np.zeros_like(psf)
    for oy in (np.arange(subsamples)+.5)/subsamples-.5:
        for ox in (np.arange(subsamples)+.5)/subsamples-.5:
            coverage+=((x+ox-cx)**2+(y+oy-cy)**2)<(radius_arcsec/pixel_arcsec)**2
    return float((psf*coverage).sum()/(psf.sum()*subsamples**2))


def notebook_factors(path):
    nb=json.loads(Path(path).read_text())
    source="\n".join("".join(c.get("source",[])) for c in nb["cells"]
                     if c.get("cell_type")=="code")
    if "pixel_scale=0.03" not in source or "norm_radius=1" not in source:
        raise ValueError("Unverified PSF scale/normalization")
    m=re.search(r"corr_factors\s*=\s*(\{[\s\S]*?\})",source)
    if not m:raise ValueError("No empirical correction dictionary")
    factors=ast.literal_eval(m.group(1))
    if "sp[f'new_e_{band}'] = sp[f'e_{band}'] * corr" not in source:
        raise ValueError("No verified new_e = e*correction equation")
    if any(len(factors[b])!=5 for b in BANDS):
        raise ValueError("Wrong correction factor lengths")
    return factors


def run(psf_path,ee_path,notebook_path,out):
    psf=read_small_psf(psf_path)
    table=np.genfromtxt(ee_path,names=True)
    radius=table["aper_radius"]
    if not np.all(np.diff(radius)>0):raise ValueError("Nonmonotonic EE radius")
    factors=notebook_factors(notebook_path)
    psf_ref=enclosed_fraction(psf,1.)
    tab_ref=float(np.interp(1.,radius,table["F444W"]))
    comparisons=[]
    for r in RADII:
        p=enclosed_fraction(psf,r)
        tab=float(np.interp(r,radius,table["F444W"])/tab_ref)
        comparisons.append({
            "radius_arcsec":r,"psf_fraction_total":p,
            "psf_relative_to_1arcsec":p/psf_ref,
            "table_relative_to_1arcsec":tab,
            "relative_difference_percent":100*(p/psf_ref/tab-1)})
    f444=3631e9*10**(-.4*26.2)
    sigma277=f444/30/2
    factor_info={b:[{
        "aperture_diameter_arcsec":d,
        "factor":float(c),
        "nominal_depth_overstatement_mag_if_new_e_used":float(2.5*np.log10(c))
    } for d,c in zip(APERTURES,factors[b])] for b in BANDS}
    requirements=[{
        "aperture_diameter_arcsec":d,
        "max_raw_e_F444W_nJy_for_corrected_SNR5":f444/(5*c444),
        "max_raw_e_F277W_nJy_for_corrected_zero_flux_2sigma_ratio30":sigma277/c277
    } for d,c444,c277 in zip(APERTURES,factors["F444W"],factors["F277W"])]
    result={
        "status":"SOURCE_DERIVED_SMALL_PSF_AND_ERROR_FACTOR_QC_NOT_SELECTION_AREA",
        "sha256":{Path(p).name:hashlib.sha256(Path(p).read_bytes()).hexdigest()
                  for p in (psf_path,ee_path,notebook_path)},
        "psf_shape":list(psf.shape),"psf_sum":float(psf.sum()),
        "psf_max_pixel_yx":list(map(int,np.unravel_index(np.argmax(psf),psf.shape))),
        "pixel_scale_arcsec":.03,"normalization_radius_arcsec":1.,
        "curve_of_growth_F444W":comparisons,
        "empirical_error_factor_by_aperture":factor_info,
        "F444W_26p2_AB_flux_nJy":f444,
        "F277W_max_flux_nJy_for_ratio30":f444/30,
        "F277W_max_corrected_sigma_nJy_if_2sigma_zero_flux_upper_limit":sigma277,
        "conditional_raw_error_thresholds":requirements,
        "limitations":[
            "Actual Dave selection e_* vs new_e_* and negative F277W flux handling remain unknown.",
            "These are conditional raw-error thresholds, not measured depths or selection completeness.",
            "PSF 8x8 subpixel integration need not reproduce exactly the original notebook aperture implementation.",
            "Notebook's new_e_* are created after use_phot; use_phot is not a substitute for a corrected SNR>5 selection.",
            "No WHT/epoch/star-mask spatial integration or F356W detection here."
        ]}
    Path(out).parent.mkdir(parents=True,exist_ok=True)
    Path(out).write_text(json.dumps(result,indent=2)+"\n")
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--psf",required=True)
    p.add_argument("--ee",required=True)
    p.add_argument("--notebook",required=True)
    p.add_argument("--out",default="pm_footprint/results/v13_psf_noise_qc.json")
    a=p.parse_args()
    r=run(a.psf,a.ee,a.notebook,a.out)
    print(json.dumps({"F444W_0p16arcsec":r["curve_of_growth_F444W"][0],
                      "conditional_error_thresholds":r["conditional_raw_error_thresholds"]},indent=2))


if __name__=="__main__":
    main()
