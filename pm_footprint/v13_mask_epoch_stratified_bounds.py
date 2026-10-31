"""Epoch-stratified conditional JWST mask x v13 WHT bounds.

Unlike v13_mask_near_containment_bounds.py (2026-10-31), this version
separates missing WHT sky INSIDE the red-epoch archive polygon from missing
WHT sky OUTSIDE it. Both can be recovered from the frozen exact archive
area and the two-phase v13 WHT overlay. No large FITS access.

This is a conditional set-theoretic bound given *estimated* mask/archive
Monte Carlo intersections, not a direct mask-WHT measurement or a 95% CI.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

BASE = Path("pm_footprint/results")


def phase_bounds(archive, archive_red, wht, wht_archive, wht_red,
                 mask_archive, mask_archive_red, mask_total_cap):
    """Extrema of (WHT_red - mask_WHT_red)/(WHT - mask_WHT).

    Define D_red = archive_red \\ WHT, D_nonred =
    (archive \\ archive_red) \\ WHT, O = WHT \\ archive.
    Let a,b,c be the masked portions of D_red,D_nonred,O.
    Then f=(wht_red-mask_archive_red+a)/
           (wht-mask_archive+a+b-c).
    f increases with a,c and decreases with b for feasible 0<f<1.
    """
    v = (archive, archive_red, wht, wht_archive, wht_red,
         mask_archive, mask_archive_red, mask_total_cap)
    if any(not math.isfinite(x) or x < 0 for x in v):
        raise ValueError("Nonfinite or negative area")
    eps = 1e-10
    if not (archive_red <= archive + eps
            and wht_archive <= min(wht, archive) + eps
            and wht_red <= min(wht_archive, archive_red) + eps
            and mask_archive_red <= min(mask_archive, archive_red) + eps
            and mask_archive <= min(archive, mask_total_cap) + eps
            and mask_archive - mask_archive_red <= archive - archive_red + eps):
        raise ValueError("Inconsistent nested regions")
    d_red = archive_red - wht_red
    d_nonred = (archive - archive_red) - (wht_archive - wht_red)
    outside = wht - wht_archive
    if min(d_red, d_nonred, outside) < -eps:
        raise ValueError("Negative archive/epoch/WHT set difference")
    amax = min(max(0., d_red), mask_archive_red)
    bmax = min(max(0., d_nonred), mask_archive - mask_archive_red)
    cmax = min(max(0., outside), mask_total_cap - mask_archive)
    # The extrema use (a,b,c)=(0,bmax,0) and (amax,0,cmax).
    # The masked sky not in WHT can be placed in the ample unobserved
    # archive complement under the stated conditional area assumptions.
    lo = (wht_red - mask_archive_red) / (wht - mask_archive + bmax)
    hi = (wht_red - mask_archive_red + amax) / (
        wht - mask_archive + amax - cmax)
    if not (0 <= lo <= hi <= 1):
        raise ValueError("No feasible positive unmasked area")
    return {
        "archive_red_not_in_wht_deg2": d_red,
        "archive_nonred_not_in_wht_deg2": d_nonred,
        "wht_outside_archive_deg2": outside,
        "max_mask_red_outside_wht_deg2": amax,
        "max_mask_nonred_outside_wht_deg2": bmax,
        "max_mask_outside_archive_in_wht_deg2": cmax,
        "conditional_mask_in_wht_red_deg2": [
            mask_archive_red - amax, mask_archive_red],
        "conditional_corrected_red_fraction": [lo, hi],
        "uncorrected_fraction": wht_red / wht,
    }


def analyze(overlay, mc, mask_summary, raster):
    archive = float(overlay["exact_archive_deg2"])
    archive_red = float(overlay["exact_red_epoch_deg2"])
    measured = mc["measured_mask_intersections_deg2"]
    mask_archive = float(measured["mask_within_archive_f444w"])
    mask_red = float(measured["mask_within_archive_red_365d"])
    native_sum = float(mask_summary["sum_per_tile_mask_code1_area_deg2_NOT_UNION"])
    raster_cap = float(
        raster["grids"]["0.6arcsec"]["union_masked_area_approx_upper_deg2"])
    cases = []
    for phase in overlay["phase_results"]:
        wht = float(phase["wht_positive_area_deg2"])
        intersection = phase["area_intersections_deg2"]
        wht_archive = float(intersection["archive"])
        wht_red = float(intersection["red_365d"])
        cases.append({
            "phase": phase["phase"],
            "nominal_native_tile_sum_cap": phase_bounds(
                archive, archive_red, wht, wht_archive, wht_red,
                mask_archive, mask_red, native_sum),
            "raster_union_cap_sensitivity": phase_bounds(
                archive, archive_red, wht, wht_archive, wht_red,
                mask_archive, mask_red, raster_cap),
        })
    def envelope(key):
        bounds = [x[key]["conditional_corrected_red_fraction"] for x in cases]
        return [min(x[0] for x in bounds), max(x[1] for x in bounds)]
    return {
        "status": "EPOCH_STRATIFIED_CONDITIONAL_SET_BOUNDS_NOT_DIRECT_MASK_WHT",
        "source_files": [
            "pm_footprint/results/local_v13_epoch_overlay_20261031.json",
            "pm_footprint/results/jwst_starmask/epoch_mc_20261031.json",
            "pm_footprint/results/jwst_starmask/summary.json",
            "pm_footprint/results/jwst_starmask/union_area_grid_qc.json"],
        "correction_to_prior": (
            "Prior v13_mask_near_containment_bounds.py used total "
            "archive-minus-WHT area independently as a red and nonred "
            "allocation cap. Exact red archive geometry allows separate, "
            "strictly smaller epoch-stratified caps."),
        "mask_total_nominal_native_tile_sum_cap_deg2": native_sum,
        "mask_total_0p6arcsec_raster_upper_approx_deg2": raster_cap,
        "cases": cases,
        "nominal_native_sum_cap_fraction_envelope": envelope(
            "nominal_native_tile_sum_cap"),
        "raster_cap_sensitivity_fraction_envelope": envelope(
            "raster_union_cap_sensitivity"),
        "limitations": [
            "Mask/archive areas use 1000 weighted MC draws; neither envelope is a statistical confidence interval.",
            "Mask code 1 interpreted as excluded; official semantics not independently verified.",
            "Native pixel tile-sum cap uses nominal 30mas pixel area and overlaps are double-counted; this is a conservative area cap, not unique-sky union.",
            "0.6arcsec raster union upper is only an approximation, not a guaranteed cap.",
            "WHT is a two-phase 3.6arcsec sampled occupancy map with boundary aliasing.",
            "Not a direct spatial mask-WHT intersection, not F444W=26.2 selection area, and not F356W astrometric completeness."
        ],
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--overlay", default=str(BASE / "local_v13_epoch_overlay_20261031.json"))
    p.add_argument("--mc", default=str(BASE / "jwst_starmask/epoch_mc_20261031.json"))
    p.add_argument("--mask", default=str(BASE / "jwst_starmask/summary.json"))
    p.add_argument("--raster", default=str(BASE / "jwst_starmask/union_area_grid_qc.json"))
    p.add_argument("--out", default=str(BASE / "jwst_starmask/v13_epoch_stratified_bounds_20261101.json"))
    a = p.parse_args()
    sources = [json.loads(Path(x).read_text()) for x in
               (a.overlay, a.mc, a.mask, a.raster)]
    result = analyze(*sources)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "native_sum_cap": result["nominal_native_sum_cap_fraction_envelope"],
        "raster_sensitivity": result["raster_cap_sensitivity_fraction_envelope"]
    }))


if __name__ == "__main__":
    main()
