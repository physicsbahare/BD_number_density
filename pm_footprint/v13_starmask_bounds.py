"""Conditional star-mask bounds from *measured* sampled v13 WHT∩epoch polygons.

Inputs are committed compact JSON (v13 spatial epoch overlay, archive
polygon areas and previous user JWST starmask Monte Carlo). Unlike naive
subtraction, this uses set inclusion in red/nonred regions independently.
No external queries or FITS downloads. Not a direct WHT∩mask pixel overlay.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path


def validate_nested(survey, red, wht, red_wht, mask_survey, mask_red):
    values = (survey, red, wht, red_wht, mask_survey, mask_red)
    if not all(isinstance(v, (float, int)) and 0 <= v for v in values):
        raise ValueError("All areas must be nonnegative")
    tol = 1e-8
    if not (red <= survey + tol and wht <= survey + tol and
            red_wht <= red + tol and red_wht <= wht + tol and
            0 <= mask_red <= mask_survey + tol and
            mask_survey <= survey + tol and
            mask_survey-mask_red <= survey-red+tol):
        raise ValueError("Areas violate S/R/W or M nested-set conditions")
    missing_red = red - red_wht
    missing_nonred = (survey - red) - (wht - red_wht)
    if missing_red < -tol or missing_nonred < -tol:
        raise ValueError("Sampled WHT red/nonred not contained in archive geometry")
    return max(0.0, missing_red), max(0.0, missing_nonred)


def bounds(survey, red, wht, red_wht, mask_survey, mask_red):
    """Sharp algebraic bounds conditional on the six supplied areas.

    S = archive; R subset S = independent red epoch; W subset S = v13
    WHT-positive sample centers inside S; M = code-1 starmask union.
    Unknown x=|M ∩ R ∩ W|, y=|M ∩ (W\\R)|.

    This is conditional on exact areas, not statistical confidence bounds.
    """
    dr, dn = validate_nested(survey, red, wht, red_wht,
                             mask_survey, mask_red)
    xlo, xhi = max(0.0, mask_red-dr), min(mask_red, red_wht)
    ylo = max(0.0, mask_survey-mask_red-dn)
    yhi = min(mask_survey-mask_red, wht-red_wht)
    if xlo > xhi+1e-10 or ylo > yhi+1e-10:
        raise ValueError("Incompatible mask/footprint allocations")
    if wht-xhi-yhi <= 0:
        raise ValueError("No positive unmasked WHT area")
    # Fraction decreases with excluded red x, increases with excluded nonred y.
    lower = (red_wht-xhi)/(wht-xhi-ylo)
    upper = (red_wht-xlo)/(wht-xlo-yhi)
    if not (0 <= lower <= upper <= 1+1e-10):
        raise AssertionError("Invalid corrected fraction bounds")
    return {
        "red_missing_from_v13_wht_deg2":dr,
        "nonred_missing_from_v13_wht_deg2":dn,
        "masked_red_wht_area_bounds_deg2":[xlo,xhi],
        "masked_nonred_wht_area_bounds_deg2":[ylo,yhi],
        "masked_total_wht_area_bounds_deg2":[xlo+ylo,xhi+yhi],
        "unmasked_red_wht_area_bounds_deg2":[red_wht-xhi,red_wht-xlo],
        "unmasked_total_wht_area_bounds_deg2":[wht-xhi-yhi,wht-xlo-ylo],
        "geometric_red_fraction_before_mask":red_wht/wht,
        "geometric_red_fraction_after_mask_bounds":[lower,upper]
    }


def analyze(overlay, mask_mc):
    s = overlay["exact_archive_deg2"]
    r = overlay["exact_red_epoch_deg2"]
    m = mask_mc["measured_mask_intersections_deg2"]
    b = m["mask_within_archive_f444w"]
    c = m["mask_within_archive_red_365d"]
    db = m["survey_mc_95pct_halfwidth_deg2"]
    dc = m["red_365d_mc_95pct_halfwidth_deg2"]
    results = []
    envelope = []
    for phase in overlay["phase_results"]:
        g = phase["area_intersections_deg2"]
        w = g["archive"]   # WHT sample clipped to original archive polygon
        wr = g["red_365d"]
        nominal = bounds(s,r,w,wr,b,c)
        sensitivities = []
        for bb,cc in itertools.product((b-db,b+db),(c-dc,c+dc)):
            sensitivities.append(bounds(s,r,w,wr,bb,cc))
        local_env = [
            min(z["geometric_red_fraction_after_mask_bounds"][0] for z in sensitivities),
            max(z["geometric_red_fraction_after_mask_bounds"][1] for z in sensitivities)
        ]
        envelope.extend([local_env[0],local_env[1]])
        results.append({
            "phase":phase["phase"],
            "wht_area_within_archive_deg2":w,
            "positive_wht_outside_archive_deg2":phase["wht_positive_area_deg2"]-w,
            "red_365d_wht_intersection_deg2":wr,
            "nominal_mc_mask_area_conditional_bounds":nominal,
            "marginal_mc_halfwidth_corner_sensitivity_fraction":local_env
        })
    phase = [p["nominal_mc_mask_area_conditional_bounds"][
        "geometric_red_fraction_after_mask_bounds"] for p in results]
    return {
        "status":"CONDITIONAL_BOUNDS_ON_SAMPLED_V13_WHT_AND_PROVISIONAL_MC_STARMASK",
        "source_interpretation":"S/R from exact JWST CAOM geometry; W/WR from user's two 3.6-arcsec-spaced v13 WHT-positive sample grids; M areas estimated by previous 1000-draw MC 0.6arcsec JWST mask sky union.",
        "mask_code_excluded":1,
        "mask_in_archive_area_deg2_mc_nominal":b,
        "mask_in_archive_red365_area_deg2_mc_nominal":c,
        "phase_results":results,
        "nominal_mc_mask_fraction_bounds_across_phases":[
            min(z[0] for z in phase),max(z[1] for z in phase)],
        "illustrative_mc_marginal_95pct_corner_sensitivity_envelope_across_phases":[
            min(envelope),max(envelope)],
        "caveats":[
            "These are SHARP SET bounds conditional on supplied numerical area values, not a measured v13 WHT∩star-mask spatial union.",
            "The mask areas come from 1000 MC points and coarse union raster; its uncertainty and systematic WCS/raster errors are NOT included in nominal set bounds.",
            "Corner sensitivity combines separate marginal 95% MC halfwidths; it is NOT a joint 95% confidence interval (and cannot account for raster/WCS systematics).",
            "Two 3.6arcsec WHT sampling phases do not bound native-pixel area uncertainty; much thinner mask edges and holes may be missed.",
            "Clipped WHT area excludes a small number of positive sample centers outside the exact archival survey polygon.",
            "Mask code 1 interpreted as excluded; documentation/semantics remain to be confirmed.",
            "None of this establishes F444W=26.2 dropout-selection completeness, F356W detectability or astrometric PM sensitivity."
        ]
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--overlay",default="pm_footprint/results/local_v13_epoch_overlay_20261031.json")
    p.add_argument("--mask-mc",default="pm_footprint/results/jwst_starmask/epoch_mc_20261031.json")
    p.add_argument("--out",default="pm_footprint/results/jwst_starmask/v13_conditional_mask_bounds_20261031.json")
    a=p.parse_args()
    inputs = {"overlay":Path(a.overlay),"mask_mc":Path(a.mask_mc)}
    d = analyze(json.loads(inputs["overlay"].read_text()),
                json.loads(inputs["mask_mc"].read_text()))
    d["input_sha256"]={k:hashlib.sha256(path.read_bytes()).hexdigest()
                       for k,path in inputs.items()}
    output=Path(a.out)
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
    print(json.dumps({
        "conditional_fraction_bounds":d["nominal_mc_mask_fraction_bounds_across_phases"],
        "illustrative_mc_corner_sensitivity":d["illustrative_mc_marginal_95pct_corner_sensitivity_envelope_across_phases"]
    },indent=2))


if __name__=="__main__":
    main()
