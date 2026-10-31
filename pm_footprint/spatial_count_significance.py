"""Exploratory binomial diagnostic for Dave125 positions; NOT PM completeness.
Run: python -m pm_footprint.spatial_count_significance
"""
import json
from pathlib import Path
from scipy.stats import binomtest

AREA_TOTAL = 0.528823684998808
AREA_RED365 = 0.31462445722188803

def conditional_test(n_inside, n_total, area_red, area_total):
    if not (0 < area_red < area_total and 0 <= n_inside <= n_total and n_total > 0):
        raise ValueError("Invalid nested areas or counts")
    p0 = area_red / area_total
    test = binomtest(n_inside, n_total, p=p0, alternative="two-sided")
    ci = test.proportion_ci(confidence_level=0.95, method="exact")
    rate_in = n_inside / area_red
    rate_out = (n_total - n_inside) / (area_total - area_red)
    def rr(q):
        return float("inf") if q == 1 else (q / (1-q)) * ((area_total-area_red)/area_red)
    return {
        "n_inside": n_inside, "n_outside": n_total-n_inside, "n_total": n_total,
        "null_inside_probability": p0, "observed_inside_fraction": n_inside/n_total,
        "density_inside_per_deg2": rate_in, "density_outside_per_deg2": rate_out,
        "inside_outside_density_ratio": rate_in/rate_out if rate_out else None,
        "conditional_exact_two_sided_pvalue": test.pvalue,
        "conditional_exact_95pct_inside_fraction_ci": [ci.low, ci.high],
        "transformed_95pct_density_ratio_ci": [rr(ci.low), rr(ci.high)]
    }

def main():
    result = {
        "status": "CONDITIONAL_UNIFORM_SURFACE_DENSITY_DIAGNOSTIC_NOT_COMPLETENESS",
        "area_total_deg2": AREA_TOTAL, "red365_area_deg2": AREA_RED365,
        "individual_catalog_entries": conditional_test(62, 125, AREA_RED365, AREA_TOTAL),
        "positional_groups_0p5arcsec": conditional_test(60, 115, AREA_RED365, AREA_TOTAL),
        "assumptions": [
            "Independent uniform source positions conditional on N and fixed areas.",
            "0.5arcsec positional grouping is NOT validated deblending.",
            "No correction for spatial clustering, selection, depth or WHT rasterization.",
            "Exploratory conditional null test, not PM detection probability or completeness."
        ]
    }
    path = Path("pm_footprint/results/candidate_spatial_significance.json")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))

if __name__ == "__main__":
    main()
