"""Propagate finite mask-sample uncertainty into v13 red-epoch area bounds.

Reuses the 1000 archived weighted mask-grid draws (575 red, 995 survey).
No large WHT downloads, no MAST query, no binary mask reprocessing.
Wilson intervals for each of two binomial proportions have 97.5% marginal
coverage; Bonferroni gives approximately >=95% simultaneous coverage
under iid sampling. Raster/WCS/WHT errors are NOT included.

Run: python -m pm_footprint.v13_epoch_mc_uncertainty
"""
import argparse
import json
import math
from pathlib import Path
from statistics import NormalDist
from pm_footprint.v13_mask_epoch_stratified_bounds import phase_bounds

BASE = Path("pm_footprint/results")


def wilson(k, n, confidence):
    if type(k) is not int or type(n) is not int or n < 1 or not 0 <= k <= n:
        raise ValueError("Require integer 0 <= k <= n, n > 0")
    if not 0 < confidence < 1:
        raise ValueError("Invalid confidence")
    z = NormalDist().inv_cdf((1 + confidence) / 2)
    p = k / n
    den = 1 + z*z/n
    mid = (p + z*z/(2*n))/den
    half = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n))/den
    return [max(0., mid-half), min(1., mid+half)]


def analyze(overlay, mc, mask_summary, raster, alpha=0.05):
    n = int(mc["provenance"]["sample_draws_used"])
    kr = int(mc["sample_membership_counts"]["red365"])
    ks = int(mc["sample_membership_counts"]["survey"])
    if not 0 <= kr <= ks <= n:
        raise ValueError("Red sample must be nested in survey sample")
    marginal = 1 - alpha/2
    pr = wilson(kr, n, marginal)
    ps = wilson(ks, n, marginal)
    area = float(raster["grids"]["0.6arcsec"]["union_masked_area_approx_upper_deg2"])
    cap = float(mask_summary["sum_per_tile_mask_code1_area_deg2_NOT_UNION"])
    A = float(overlay["exact_archive_deg2"])
    AR = float(overlay["exact_red_epoch_deg2"])
    phases = []
    for phase in overlay["phase_results"]:
        W = float(phase["wht_positive_area_deg2"])
        WA = float(phase["area_intersections_deg2"]["archive"])
        WR = float(phase["area_intersections_deg2"]["red_365d"])
        nominal = phase_bounds(A, AR, W, WA, WR, area*ks/n, area*kr/n, cap)
        corners = []
        for sp in ps:
            for rp in pr:
                if rp > sp:
                    raise ValueError("Infeasible confidence rectangle")
                b = phase_bounds(A, AR, W, WA, WR, area*sp, area*rp, cap)
                corners.append(b["conditional_corrected_red_fraction"])
        # Both extrema decrease with red probability and increase with
        # survey probability (piecewise monotone min-cap expressions).
        envelope = [min(x[0] for x in corners), max(x[1] for x in corners)]
        if cap >= min(WR, W-WR):
            raise ValueError("Mask cap too large for simple all-placement bounds")
        # Independent of MC, but conditional on the nominal native-pixel cap.
        placement = [(WR-cap)/(W-cap), WR/(W-cap)]
        phases.append({
            "phase": phase["phase"], "unmasked_fraction": WR/W,
            "nominal_conditional_fraction": nominal["conditional_corrected_red_fraction"],
            "mc_wilson_bonferroni_fraction_envelope": envelope,
            "all_mask_placement_cap_fraction_envelope": placement})
    return {
        "status": "MC_UNCERTAINTY_AUDIT_NOT_DIRECT_WHT_MASK_INTERSECTION",
        "n_mask_draws": n, "red_count": kr, "survey_count": ks,
        "mask_raster_area_deg2": area, "native_tile_sum_cap_deg2": cap,
        "marginal_wilson_confidence": marginal,
        "bonferroni_nominal_joint_confidence": 1-alpha,
        "red_membership_wilson": pr, "survey_membership_wilson": ps,
        "phase_results": phases,
        "mc_fraction_envelope_all_phases": [
            min(x["mc_wilson_bonferroni_fraction_envelope"][0] for x in phases),
            max(x["mc_wilson_bonferroni_fraction_envelope"][1] for x in phases)],
        "placement_cap_fraction_envelope_all_phases": [
            min(x["all_mask_placement_cap_fraction_envelope"][0] for x in phases),
            max(x["all_mask_placement_cap_fraction_envelope"][1] for x in phases)],
        "limitations": [
            "Wilson/Bonferroni covers iid sample noise only, approximately; not a total scientific confidence interval.",
            "0.6arcsec raster upper is not a rigorous exact-sky upper bound.",
            "Nominal 30mas tile-sum cap is not an independently calibrated WCS pixel-area bound.",
            "WHT 3.6arcsec grid phase, area deposition losses and WCS systematics are excluded.",
            "Not a direct mask-WHT intersection, F444W=26.2 dropout area or F356W PM completeness."]}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--overlay", default=str(BASE/"local_v13_epoch_overlay_20261031.json"))
    p.add_argument("--mc", default=str(BASE/"jwst_starmask/epoch_mc_20261031.json"))
    p.add_argument("--mask", default=str(BASE/"jwst_starmask/summary.json"))
    p.add_argument("--raster", default=str(BASE/"jwst_starmask/union_area_grid_qc.json"))
    p.add_argument("--out", default=str(BASE/"jwst_starmask/v13_epoch_mc_uncertainty_20261101.json"))
    a = p.parse_args()
    sources = [json.loads(Path(x).read_text()) for x in (a.overlay,a.mc,a.mask,a.raster)]
    result = analyze(*sources)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"MC":result["mc_fraction_envelope_all_phases"],
                      "all_placement":result["placement_cap_fraction_envelope_all_phases"]}))


if __name__ == "__main__":
    main()
