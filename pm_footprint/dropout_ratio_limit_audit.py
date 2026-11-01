"""Audit whether F277W nondetection *limits* certify a F444W/F277W > 30 dropout.

Uses the frozen 125-row v13 photometry fixture only; does NOT measure sky area,
COSMOS2025 selection completeness, or astrometric completeness. Fluxes are nJy.
No network, MAST, WHT image or star-mask archive required.

Run: python -m pm_footprint.dropout_ratio_limit_audit
"""
import argparse
import csv
import json
import math
from pathlib import Path

AB_NJY_ZEROPOINT = 31.4
RATIO = 30.0
CONFIRMED_PM = (387132, 1651050)


def flux_njy(mag):
    return 10.0 ** ((AB_NJY_ZEROPOINT - mag) / 2.5)


def required_5sigma_f277_depth(m444=26.2, ratio=RATIO, limit_sigma=2.0):
    """AB 5sigma depth needed for ratio lower bound > ratio from limit_sigma*err."""
    if not math.isfinite(m444) or ratio <= 0 or limit_sigma <= 0:
        raise ValueError("Invalid depth inputs")
    return m444 + 2.5 * math.log10(ratio * limit_sigma / 5.0)


def median(values):
    a = sorted(values)
    n = len(a)
    if not n:
        raise ValueError("Empty group")
    return (a[(n - 1) // 2] + a[n // 2]) / 2.0


def audit(rows):
    rows = list(rows)
    ids = [int(row["candidate_id"]) for row in rows]
    if len(ids) != 125 or len(set(ids)) != 125 or not set(CONFIRMED_PM).issubset(ids):
        raise ValueError("Expected frozen 125 unique IDs including both PM objects")
    samples = []
    for row in rows:
        cid = int(row["candidate_id"])
        f444 = float(row["f_F444W"])
        e444 = float(row["new_e_F444W"])
        f277 = float(row["f_F277W"])
        e277 = float(row["new_e_F277W"])
        m444 = float(row["mag_F444W"])
        epoch = int(row["red365"])
        if (not all(map(math.isfinite, (f444, e444, f277, e277, m444)))
                or f444 <= 0 or e444 <= 0 or e277 <= 0 or epoch not in (0, 1)):
            raise ValueError(f"Invalid photometry for {cid}")
        if abs(m444 - (AB_NJY_ZEROPOINT - 2.5 * math.log10(f444))) > 1e-5:
            raise ValueError(f"AB flux/magnitude mismatch for {cid}; expected nJy")
        snr277 = f277 / e277
        samples.append({
            "id": cid, "mag": m444, "red365": epoch,
            "snr277": snr277, "f277": f277, "e277": e277,
            "ratio_from_zero_centered_2sigma_limit": f444 / (2.0 * e277),
            "ratio_from_positive_flux_plus_2sigma": f444 / (max(0.0, f277) + 2.0 * e277),
            "f277_5sigma_depth": AB_NJY_ZEROPOINT - 2.5 * math.log10(5.0 * e277),
        })

    def group(seq):
        return {
            "n": len(seq),
            "f277_signed_snr_lt2": sum(s["snr277"] < 2 for s in seq),
            "f277_abs_snr_lt2": sum(abs(s["snr277"]) < 2 for s in seq),
            "f277_negative_flux": sum(s["f277"] < 0 for s in seq),
            "n_ratio_gt30_from_zero_centered_2sigma": sum(
                s["ratio_from_zero_centered_2sigma_limit"] > RATIO for s in seq),
            "n_ratio_gt30_from_positive_flux_plus_2sigma": sum(
                s["ratio_from_positive_flux_plus_2sigma"] > RATIO for s in seq),
            "median_f277_sigma_njy": median([s["e277"] for s in seq]),
            "median_f277_5sigma_depth_ab": median([s["f277_5sigma_depth"] for s in seq]),
            "median_ratio_from_zero_centered_2sigma": median(
                [s["ratio_from_zero_centered_2sigma_limit"] for s in seq]),
        }

    groups = {
        "all": group(samples),
        "red365": group([s for s in samples if s["red365"] == 1]),
        "nonred365": group([s for s in samples if s["red365"] == 0]),
        "f444w_bright_lt25": group([s for s in samples if s["mag"] < 25]),
        "f444w_25_to_26p2": group([s for s in samples if 25 <= s["mag"] <= 26.2]),
    }
    pm = {}
    for s in samples:
        if s["id"] in CONFIRMED_PM:
            pm[str(s["id"])] = {
                "f444w_ab": s["mag"],
                "red365": s["red365"],
                "f277w_sigma_njy": s["e277"],
                "ratio_lower_bound_zero_centered_2sigma": s["ratio_from_zero_centered_2sigma_limit"],
                "ratio_bound_positive_flux_plus_2sigma": s["ratio_from_positive_flux_plus_2sigma"],
                "f277w_5sigma_depth_ab": s["f277_5sigma_depth"],
            }
    return {
        "status": "V13_CANDIDATE_LEVEL_LIMIT_AUDIT_NOT_SELECTION_AREA_OR_COSMOS2025_RESELECTION",
        "photometry_source": "pm_footprint/results/v13_epoch_candidate_photometry_fixture_20261031.csv",
        "flux_unit": "nJy (validated against AB magnitudes)",
        "f444w_reference_ab": 26.2,
        "f444w_reference_flux_njy": flux_njy(26.2),
        "required_f277w_sigma_njy_for_ratio_gt30_at_26p2": flux_njy(26.2) / (2.0 * RATIO),
        "required_f277w_5sigma_depth_ab_for_ratio_gt30_at_26p2": required_5sigma_f277_depth(),
        "groups": groups,
        "confirmed_pm_objects": pm,
        "interpretation": [
            "Signed F277W SNR<2 and measured F444W/F277W>30 are NOT equivalent to a ratio lower bound >30 using a 2sigma upper limit.",
            "The zero-centered 2sigma convention uses 2*ERR; the alternative max(0,measured flux)+2*ERR is a deliberately more conservative diagnostic, not an exact posterior bound.",
            "F277W errors here come from the independent v13 fixture, not necessarily Dave's original COSMOS2025 aperture errors or official spatially calibrated depths.",
            "The ratio test is not a proposed veto of any of Dave's 125 sources; it evaluates strength of their color-limit evidence.",
            "A 5sigma depth requirement of about 28.90 AB at F444W=26.2 is conditional on a 2sigma zero-centered F277W limit and ratio 30.",
            "Candidate-level errors cannot establish any sky-area completeness or F356W astrometric sensitivity.",
        ],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fixture", default="pm_footprint/results/v13_epoch_candidate_photometry_fixture_20261031.csv")
    parser.add_argument("--out", default="pm_footprint/results/dropout_ratio_limit_audit_20261101.json")
    args = parser.parse_args()
    with open(args.fixture, newline="") as f:
        result = audit(csv.DictReader(f))
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"required_f277w_5sigma_depth": result["required_f277w_5sigma_depth_ab_for_ratio_gt30_at_26p2"],
                      "groups": result["groups"], "confirmed_pm": result["confirmed_pm_objects"]}, indent=2))


if __name__ == "__main__":
    main()
