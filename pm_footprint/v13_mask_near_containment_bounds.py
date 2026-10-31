"""Set-theoretic v13 JWST-mask × WHT × red-epoch bounds without rereading WHT.

Inputs: saved two-phase v13 WHT/archive overlap, archive-only star-mask
Monte Carlo and provisional mask union raster. NOT a direct spatial
intersection, calibrated dropout area or PM completeness.

Run: python -m pm_footprint.v13_mask_near_containment_bounds
"""
import argparse
import json
from pathlib import Path

BASE = Path("pm_footprint/results")


def bound_phase(archive_area, wht_area, wht_archive_area, wht_red_area,
                mask_archive_area, mask_archive_red_area, mask_union_upper_area):
    """Feasible fraction extrema from three unobserved disjoint areas.

    D=archive\\WHT, O=WHT\\archive; a=mask∩red∩D,
    b=mask∩(archive\\red)∩D, c=mask∩O.
    corrected_red=R-mask_archive_red+a
    corrected_wht=W-mask_archive+a+b-c.
    a+b<=|D|; c<=min(|O|, mask_union-mask_archive).
    """
    vals = (archive_area, wht_area, wht_archive_area, wht_red_area,
            mask_archive_area, mask_archive_red_area, mask_union_upper_area)
    if any(not isinstance(x, (int, float)) or x < 0 for x in vals):
        raise ValueError("All areas must be nonnegative numbers")
    if not (wht_red_area <= wht_archive_area <= min(archive_area, wht_area)):
        raise ValueError("Invalid WHT/archive/red nesting")
    if not (mask_archive_red_area <= mask_archive_area <= mask_union_upper_area):
        raise ValueError("Invalid mask nesting")
    missing = archive_area - wht_archive_area
    outside = wht_area - wht_archive_area
    amax = min(missing, mask_archive_red_area)
    bmax = min(missing, mask_archive_area - mask_archive_red_area)
    cmax = min(outside, mask_union_upper_area - mask_archive_area)
    lower = (wht_red_area - mask_archive_red_area) / (wht_area - mask_archive_area + bmax)
    upper = (wht_red_area - mask_archive_red_area + amax) / (wht_area - mask_archive_area + amax - cmax)
    if not 0 <= lower <= upper <= 1:
        raise ValueError("Infeasible corrected fraction")
    return {
        "archive_not_in_wht_deg2": missing,
        "wht_outside_archive_deg2": outside,
        "max_mask_wht_outside_archive_deg2": cmax,
        "conditional_mask_wht_red_intersection_deg2": [mask_archive_red_area-amax, mask_archive_red_area],
        "conditional_mask_wht_intersection_deg2": [mask_archive_area-min(missing, mask_archive_area), mask_archive_area+cmax],
        "conditional_corrected_red_fraction": [lower, upper],
        "uncorrected_red_fraction": wht_red_area / wht_area,
    }


def analyze(overlay, mask_mc, mask_grid):
    archive = float(overlay["exact_archive_deg2"])
    m = mask_mc["measured_mask_intersections_deg2"]
    ma, mr = float(m["mask_within_archive_f444w"]), float(m["mask_within_archive_red_365d"])
    da, dr = float(m["survey_mc_95pct_halfwidth_deg2"]), float(m["red_365d_mc_95pct_halfwidth_deg2"])
    mu = float(mask_grid["grids"]["0.6arcsec"]["union_masked_area_approx_upper_deg2"])
    phases, sensitivity = [], []
    for p in overlay["phase_results"]:
        w = float(p["wht_positive_area_deg2"])
        wa = float(p["area_intersections_deg2"]["archive"])
        r = float(p["area_intersections_deg2"]["red_365d"])
        v = bound_phase(archive, w, wa, r, ma, mr, mu)
        phases.append({"phase": p["phase"], **v})
        for ma0 in (ma-da, ma+da):
            for mr0 in (mr-dr, mr+dr):
                vv = bound_phase(archive, w, wa, r, ma0, mr0, mu)
                sensitivity.extend(vv["conditional_corrected_red_fraction"])
    return {
        "status": "CONDITIONAL_V13_NEAR_CONTAINMENT_BOUNDS_NOT_DIRECT_WHT_MASK_MEASUREMENT",
        "source_files": [
            "pm_footprint/results/local_v13_epoch_overlay_20261031.json",
            "pm_footprint/results/jwst_starmask/epoch_mc_20261031.json",
            "pm_footprint/results/jwst_starmask/union_area_grid_qc.json",
        ],
        "phase_bounds": phases,
        "nominal_corrected_red_fraction_envelope": [
            min(x["conditional_corrected_red_fraction"][0] for x in phases),
            max(x["conditional_corrected_red_fraction"][1] for x in phases),
        ],
        "separate_MC_halfwidth_corner_sensitivity_envelope_NOT_95CI": [min(sensitivity), max(sensitivity)],
        "mask_archive_MC_95_halfwidth_deg2": da,
        "mask_red_MC_95_halfwidth_deg2": dr,
        "limitations": [
            "Not a direct spatial mask-WHT intersection; assumes mask code 1 means excluded.",
            "0.6arcsec union upper raster approximation is not a rigorous geometric bound.",
            "Mask MC 1000 draws and WHT 3.6arcsec phase sampling have unmodeled systematics.",
            "MC endpoint perturbations are NOT a joint 95pct confidence interval.",
            "Not 26.2AB dropout selection area, F356W detection or PM completeness.",
        ],
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--overlay", default=str(BASE/"local_v13_epoch_overlay_20261031.json"))
    p.add_argument("--mask-mc", default=str(BASE/"jwst_starmask/epoch_mc_20261031.json"))
    p.add_argument("--mask-grid", default=str(BASE/"jwst_starmask/union_area_grid_qc.json"))
    p.add_argument("--out", default=str(BASE/"jwst_starmask/v13_near_containment_bounds_20261101.json"))
    a = p.parse_args()
    result = analyze(*(json.loads(Path(x).read_text()) for x in (a.overlay, a.mask_mc, a.mask_grid)))
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result["nominal_corrected_red_fraction_envelope"]))


if __name__ == "__main__":
    main()
