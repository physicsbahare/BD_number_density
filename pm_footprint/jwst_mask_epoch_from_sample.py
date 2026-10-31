"""Masked-sky intersection with frozen MAST epoch polygons; NO WHT download.

First generate a weighted sample from the user's local JWST mask archive:
 python -m pm_footprint.jwst_mask_union_repro --archive FILE.tar.gz \
   --output mask_union.json --sample-out mask_points.csv --samples 10000
Then:
 python -m pm_footprint.jwst_mask_epoch_from_sample \
   --points mask_points.csv --mask-summary mask_union.json \
   --output epoch_mask_intersection.json

Outputs are archive-geometry intersections, NOT dropout or PM completeness.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import numpy as np
import pandas as pd
from shapely import STRtree, points
from .exact_geometry import polygon_region

BASELINES=(30,180,365,730)
RED=("F277W","F356W","F410M","F444W")


def epoch_flags(starts,ends,is_f444,is_primary,baselines=BASELINES):
    """Check F444W plus any red-band independent epoch, using interval gap."""
    starts=np.asarray(starts,float);ends=np.asarray(ends,float)
    f=np.asarray(is_f444,bool);primary=np.asarray(is_primary,bool)
    survey=bool(primary.any())
    flags={int(g):{"red":False,"same_f444w":False} for g in baselines}
    if not survey or not f.any():return survey,flags
    sep=np.maximum(starts[f,None]-ends[None,:],starts[None,:]-ends[f,None])
    for gap in baselines:
        flags[int(gap)]["red"]=bool((sep>=gap).any())
        flags[int(gap)]["same_f444w"]=bool((sep[:,f]>=gap).any())
    return survey,flags


def wht_fraction_bounds(S,R,W,WR,b,c):
    """Sharp set bounds if WHT W is a subset of archive S; R is red area.

    b=|star-mask ∩ S|, c=|star-mask ∩ R|, R subset S.
    x=|star-mask ∩ W ∩ R|; y=|star-mask ∩ W ∩ (S\\R)|.
    Since S\\W is known in R and its complement, bound x,y separately.
    """
    if not (0<=c<=b<=S and 0<=WR<=R<=S and 0<=WR<=W<=S):
        raise ValueError("Inconsistent nested areas")
    delta_r=R-WR;delta_not=(S-W)-delta_r
    if delta_not< -1e-8:raise ValueError("WHT red not nested in archive")
    xlo,xhi=max(0,c-delta_r),c
    ylo,yhi=max(0,b-c-delta_not),b-c
    if W-xhi-ylo<=0 or W-xlo-yhi<=0:
        raise ValueError("No positive unmasked WHT denominator")
    return {"fraction_lower":(WR-xhi)/(W-xhi-ylo),
            "fraction_upper":(WR-xlo)/(W-xlo-yhi),
            "star_mask_in_wht_red_deg2_bounds":[xlo,xhi],
            "star_mask_in_wht_nonred_deg2_bounds":[ylo,yhi],
            "assumption":"WHT-positive grid is nested inside the exact archival F444W survey polygon; star-mask areas are Monte Carlo estimates"}


def evaluate(points_csv,mask_summary,metadata_csv,exact_json,wht_json):
    pts=pd.read_csv(points_csv)
    if list(pts.columns)!=["x_lambert_rad","y_lambert_rad"]:
        raise ValueError("Expected Lambert equal-area sample columns in radians")
    xy=points(pts.x_lambert_rad.to_numpy(),pts.y_lambert_rad.to_numpy())
    n=len(pts)
    if n<100:raise ValueError("Need >=100 weighted mask samples")
    area=mask_summary.get("grid_results",{}).get("0.6",{}).get("union_upper_bound_deg2")
    if area is None:
        area=mask_summary["grids"]["0.6arcsec"]["union_masked_area_approx_upper_deg2"]
    area=float(area)
    raw=pd.read_csv(metadata_csv)
    frame=raw.loc[(raw.obs_collection.astype(str).str.upper()=="JWST")
        &(raw.instrument_name.astype(str).str.upper()=="NIRCAM/IMAGE")
        &(raw.dataproduct_type.astype(str).str.lower()=="image")
        &(raw.intentType.astype(str).str.lower()=="science")
        &(raw.dataRights.astype(str).str.upper()=="PUBLIC")
        &raw.filters.astype(str).str.upper().str.contains("|".join(RED))].copy()
    frame["t_min"]=pd.to_numeric(frame.t_min,errors="coerce")
    frame["t_max"]=pd.to_numeric(frame.t_max,errors="coerce")
    frame=frame[np.isfinite(frame.t_min)&np.isfinite(frame.t_max)&(frame.t_max>=frame.t_min)]
    poly=[polygon_region(s) for s in frame.s_region]
    ii,jj=STRtree(poly).query(xy,predicate="intersects")
    matches=[[] for _ in range(n)]
    for i,j in zip(ii,jj):matches[int(i)].append(int(j))
    f444=frame.filters.astype(str).str.upper().str.contains("F444W").to_numpy()
    primary=f444&frame.proposal_id.astype(str).str.replace(r"\.0$","",regex=True).eq("1727").to_numpy()
    starts,ends=frame.t_min.to_numpy(),frame.t_max.to_numpy()
    survey=np.zeros(n,bool);red={g:np.zeros(n,bool) for g in BASELINES}
    same={g:np.zeros(n,bool) for g in BASELINES}
    for i,match in enumerate(matches):
        if not match:continue
        q=np.asarray(match,dtype=int)
        survey[i],flags=epoch_flags(starts[q],ends[q],f444[q],primary[q])
        for gap in BASELINES:
            red[gap][i]=flags[gap]["red"]
            same[gap][i]=flags[gap]["same_f444w"]
    b=area*survey.mean();c=area*red[365].mean()
    S=exact_json["archive_f444w_area_deg2"]
    R=exact_json["baseline_results"]["365"]["area_deg2"]["f444w_plus_red"]
    W=wht_json["area_deg2"]["0"]["archive_f444w"]
    WR=wht_json["area_deg2"]["0"]["f444w_plus_red_365d"]
    out={"status":"PROVISIONAL_STAR_MASK_ARCHIVE_INTERSECTION_NOT_SELECTION_MATCHED",
        "method":"weighted Monte Carlo on 0.6arcsec mask grid; exact MAST CAOM polygon point membership; conservative interval-separated epochs",
        "n_samples":n,"n_valid_red_imaging_rows":len(frame),
        "mask_union_upper_raster_area_deg2":area,
        "mask_in_archive":{"count":int(survey.sum()),"area_deg2":b},
        "epoch_intersections":{},
        "archive_only_mask_corrected_365d":{"denominator_deg2":S-b,
            "red_numerator_deg2":R-c,"red_fraction":(R-c)/(S-b)},
        "conditional_wht_mask_corrected_365d":wht_fraction_bounds(S,R,W,WR,b,c),
        "caveats":["WHT pixel occupancy was not saved; bounds are not a measured WHT-star-mask spatial intersection.",
            "0.6arcsec raster and finite weighted samples introduce systematic/statistical uncertainty.",
            "Value 1 interpreted as excluded, pending official mask semantics verification.",
            "No F444W=26.2 dropout selection completeness or F356W astrometric SNR modeled."]}
    for gap in BASELINES:
        k=int(red[gap].sum());p=k/n
        out["epoch_intersections"][str(gap)]={"red_count":k,"same_f444w_count":int(same[gap].sum()),
            "red_mask_area_deg2":area*p,"same_f444w_mask_area_deg2":area*same[gap].mean(),
            "red_mc_95pct_halfwidth_deg2":1.96*area*math.sqrt(p*(1-p)/n)}
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--points",required=True)
    p.add_argument("--mask-summary",required=True)
    p.add_argument("--metadata",default="pm_footprint/results/mast_observations.csv")
    p.add_argument("--exact",default="pm_footprint/results/exact/summary.json")
    p.add_argument("--wht",default="pm_footprint/results/wht60/summary.json")
    p.add_argument("--output",default="pm_footprint/results/jwst_starmask/epoch_intersection.json")
    a=p.parse_args()
    res=evaluate(a.points,json.loads(Path(a.mask_summary).read_text()),a.metadata,
        json.loads(Path(a.exact).read_text()),json.loads(Path(a.wht).read_text()))
    Path(a.output).parent.mkdir(parents=True,exist_ok=True)
    Path(a.output).write_text(json.dumps(res,indent=2)+"\n")
    print(json.dumps({"archive_only":res["archive_only_mask_corrected_365d"],
                      "conditional_wht":res["conditional_wht_mask_corrected_365d"]},indent=2))

if __name__=="__main__":main()
