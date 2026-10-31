"""Conservative test: can the JWST star mask alone reconcile red-epoch source densities?

Uses frozen WHT geometric areas and a generous UPPER bound on the unique
code-1 mask area (sum over all 20 native FITS tiles, including duplicates).
No WHT download, mask rerun, or Monte Carlo is required.
"""
import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def reconcile(red_area, total_area, n_red, n_nonred, mask_upper):
    """Extrema over red/nonred excluded areas x,y>=0 with x+y<=mask_upper.

    Density ratio D=(n_red/(R-x))/(n_nonred/(T-R-y)).
    It increases with x and decreases with y. Extremes are (0,M),(M,0).
    This tests a hypothetical uniform-density explanation, NOT uniformity.
    """
    R, T, M = map(float, (red_area, total_area, mask_upper))
    O = T-R
    if not (0 < R < T and 0 <= M < min(R, O)):
        raise ValueError("Require nested positive areas and M smaller than both components")
    if not (isinstance(n_red, int) and isinstance(n_nonred, int)
            and n_red > 0 and n_nonred > 0):
        raise ValueError("Require positive integer counts")
    q = n_red/n_nonred
    raw = q*O/R
    lower = q*(O-M)/R
    upper = q*O/(R-M)
    required = R-q*O  # red-only masking needed to make D(x,0)=1
    return {
        "n_red": n_red, "n_nonred": n_nonred,
        "observed_red_count_fraction": n_red/(n_red+n_nonred),
        "red_area_deg2": R, "nonred_area_deg2": O,
        "unmasked_density_ratio_red_over_nonred": raw,
        "masked_density_ratio_absolute_bounds": [lower, upper],
        "red_area_fraction_absolute_bounds": [(R-M)/(T-M), R/(T-M)],
        "red_only_mask_area_required_to_equalize_densities_deg2":
            max(0.0, required) if raw < 1 else None,
        "required_red_only_mask_over_available_upper_bound":
            max(0.0, required)/M if raw < 1 and M > 0 else None,
        "can_star_mask_alone_equalize_densities": lower <= 1 <= upper,
        "minimum_uniform_expected_red_count_after_mask":
            (n_red+n_nonred)*(R-M)/(T-M),
    }


def analyze(wht, mask, spatial):
    M = float(mask["sum_per_tile_mask_code1_area_deg2_NOT_UNION"])
    if mask["n_mask_fits"] != 20 or mask["candidate_centers_in_code1"]:
        raise ValueError("Frozen JWST mask provenance/center check changed")
    scenarios = {}
    for threshold in ("0", "0.1", "0.5"):
        a, c = wht["area_deg2"][threshold], wht["candidate_counts"][threshold]
        T, N = float(a["archive_f444w"]), int(c["archive_f444w"])
        scenarios[threshold] = {}
        for days in (30, 180, 365, 730):
            key = f"f444w_plus_red_{days}d"
            nr = int(c[key])
            scenarios[threshold][str(days)] = reconcile(float(a[key]), T, nr, N-nr, M)
    g = spatial["positional_groups_0p5arcsec"]
    grouped = reconcile(
        float(wht["area_deg2"]["0"]["f444w_plus_red_365d"]),
        float(wht["area_deg2"]["0"]["archive_f444w"]),
        int(g["n_inside"]), int(g["n_outside"]), M)
    return {
        "status": "DETERMINISTIC_MASK_INSUFFICIENCY_FOR_UNIFORM_COUNT_RECONCILIATION",
        "mask_code_interpretation": "1=excluded, conditional; confirm official binary semantics",
        "mask_native_tile_sum_area_deg2_UPPER_BOUND": M,
        "candidate_centers_in_mask_code1": 0,
        "wht_threshold_baseline_days": scenarios,
        "positional_groups_0p5arcsec_365d_at_wht0": grouped,
        "interpretation": "Even if all native code-1 mask area were unique, inside WHT and allocated adversarially, masking alone cannot reconcile observed source counts with uniform surface density. This is NOT evidence of an astrophysical underdensity; depth, source selection, clustering and deblending remain unmodeled.",
        "limitations": [
            "Conditional on 60mas WHT gridded areas (~0.96% per-tile deposition loss).",
            "Mask M is a per-tile native-pixel SUM, not a direct unique-sky or WHT overlap.",
            "Source-center unmasked does not establish unmasked PSF/aperture.",
            "Preselected F444W-only candidates are not a spatially complete population.",
            "0.5arcsec proximity groups are not verified distinct astrophysical objects.",
            "No 26.2AB selection-matched area, F356W SNR, or astrometric completeness.",
        ],
        "reproduce": "python -m pm_footprint.mask_count_reconciliation --out pm_footprint/results/mask_count_reconciliation_20261101.json",
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--wht", type=Path, default=ROOT/"pm_footprint/results/wht60/summary.json")
    p.add_argument("--mask", type=Path, default=ROOT/"pm_footprint/results/jwst_starmask/summary.json")
    p.add_argument("--counts", type=Path, default=ROOT/"pm_footprint/results/candidate_spatial_significance.json")
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    result = analyze(json.loads(a.wht.read_text()), json.loads(a.mask.read_text()),
                     json.loads(a.counts.read_text()))
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result["wht_threshold_baseline_days"]["0"]["365"], indent=2))


if __name__ == "__main__":
    main()
