"""Conditional Gaussian F277W selection probabilities for F444W-only dropouts.

Diagnostic only: NOT a calibrated COSMOS-Web effective area. Assumes
nominal F277W 5-sigma depths and fixed, noise-free F444W magnitude.
"""
import argparse
import json
import math
from pathlib import Path

def normal_cdf(x):
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))

def ab_ujy(m):
    return 10 ** ((23.9 - m) / 2.5)

def evaluate(m444=26.2, m5_depths=(27.71, 28.28, 28.898),
             true_ratios=(float("inf"), 100, 30, 10)):
    f444 = ab_ujy(m444)
    threshold = f444 / 30
    m2_needed = m444 + 2.5 * math.log10(30)
    result = {
        "status": "CONDITIONAL_NOISE_DIAGNOSTIC_NOT_MEASURED_COMPLETENESS",
        "m444_ab": m444,
        "f444_ujy": f444,
        "ratio_threshold": 30,
        "f277_measured_upper_threshold_ujy": threshold,
        "f277_2sigma_limit_needed_ab": m2_needed,
        "f277_5sigma_depth_needed_ab": m2_needed - 2.5 * math.log10(5/2),
        "depth_scenarios": [],
    }
    for m5 in m5_depths:
        sigma = ab_ujy(m5) / 5
        row = {
            "f277_5sigma_depth_ab": m5,
            "f277_1sigma_ujy": sigma,
            "f444_over_30_in_sigma": threshold / sigma,
            "true_ratio_scenarios": [],
        }
        for ratio in true_ratios:
            mu = 0 if math.isinf(ratio) else f444 / ratio
            upper = min(threshold, 2 * sigma)
            # If selection computes F444W / measured F277W and requires >30,
            # a negative F277W measurement FAILS the signed ratio inequality.
            # The F277W SNR<2 criterion also caps the upper measured flux.
            positive_ratio = normal_cdf((upper-mu)/sigma) - normal_cdf(-mu/sigma)
            # Alternative *noise-threshold* policy allowing negative measured
            # fluxes as passing the upper-flux threshold, NOT a secure
            # two-sigma color lower limit.
            threshold_pass = normal_cdf((threshold-mu)/sigma)
            row["true_ratio_scenarios"].append({
                "true_f444_f277_ratio": "infinite" if math.isinf(ratio) else ratio,
                "true_f277_ujy": mu,
                "p_positive_measured_f277_and_ratio_gt30_and_f277_snr_lt2": positive_ratio,
                "p_measured_f277_below_f444_over30_allowing_negative": threshold_pass,
                "p_measured_f277_negative": normal_cdf(-mu/sigma),
            })
        result["depth_scenarios"].append(row)
    result["limitations"] = [
        "Gaussian F277W flux noise only; fixed noise-free F444W, no correlated noise, aperture or extraction effects.",
        "Nominal planned depth scenarios, not calibrated local COSMOS2025 DR1 F277W noise maps.",
        "Actual Dave code for negative flux and upper-limit treatment must be inspected; these are alternative hypothetical policies.",
        "A 2sigma upper-limit color proof is different from a measured-flux ratio selection.",
        "These conditional probabilities are not a selection-weighted survey area or proper-motion completeness.",
    ]
    return result

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    Path(args.output).write_text(json.dumps(evaluate(), indent=2) + "\n")

if __name__ == "__main__":
    main()
