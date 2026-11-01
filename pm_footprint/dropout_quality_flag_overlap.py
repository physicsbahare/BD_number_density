"""Audit *potential* conflict between F444W-only selection and COSMOS2025 warn_flag=4.

Source for COSMOS2025 rule: https://cosmos2025.iap.fr/catalog.html
The rule is delta_mag < 0.2 in ONE NIRCam band, other NIRCam+HST+MIRI
bands delta_mag > 0.5 or undetected. We test only the F444W half of this
rule using the frozen *v13* candidate photometry fixture and corrected
aperture errors. These are NOT actual COSMOS2025 flag values.

Run: python -m pm_footprint.dropout_quality_flag_overlap
"""
import argparse
import csv
import json
import math
from pathlib import Path

K = 2.5 / math.log(10.0)
BANDS = ("F115W", "F150W", "F277W")
MOVERS = (387132, 1651050)


def audit(rows):
    """Return counts without reclassifying any candidate as an artifact."""
    rows = list(rows)
    ids = [int(row["candidate_id"]) for row in rows]
    if len(rows) != 125 or len(set(ids)) != 125:
        raise ValueError("Expected frozen 125 distinct Dave IDs")
    if not set(MOVERS).issubset(ids):
        raise ValueError("Both confirmed moving-source IDs must be present")
    threshold = K / 0.2
    bright, marginal, strong_negative = [], [], []
    movers = {}
    by_epoch = {"red365": {"n": 0, "f444_mag_error_lt_0p2": 0},
                "nonred365": {"n": 0, "f444_mag_error_lt_0p2": 0}}
    for row in rows:
        cid = int(row["candidate_id"])
        flux = float(row["f_F444W"])
        err = float(row["new_e_F444W"])
        if not math.isfinite(flux) or not math.isfinite(err) or flux <= 0 or err <= 0:
            raise ValueError(f"Invalid F444W flux/error for {cid}")
        snr = flux / err
        dm = K / snr
        group = "red365" if int(row["red365"]) == 1 else "nonred365"
        by_epoch[group]["n"] += 1
        if dm < 0.2:
            bright.append(cid)
            by_epoch[group]["f444_mag_error_lt_0p2"] += 1
        elif snr > 5:
            marginal.append(cid)
        for band in BANDS:
            blue_flux = float(row[f"f_{band}"])
            blue_err = float(row[f"new_e_{band}"])
            if not math.isfinite(blue_flux) or not math.isfinite(blue_err) or blue_err <= 0:
                raise ValueError(f"Invalid {band} flux/error for {cid}")
            signed_snr = blue_flux / blue_err
            if signed_snr <= -2:
                strong_negative.append({"candidate_id": cid, "band": band,
                                        "signed_snr": signed_snr})
        if cid in MOVERS:
            movers[str(cid)] = {"f444w_snr": snr,
                                "approx_f444w_mag_error": dm,
                                "passes_f444w_half_of_warn4": dm < 0.2}
    return {
        "status": "POTENTIAL_FLAG4_SELECTION_OVERLAP_NOT_ACTUAL_COSMOS2025_FLAGS",
        "source_catalog": "Zhen-Kai v13 frozen Dave125 photometry fixture",
        "official_rule": "COSMOS2025 warn_flag=4: delta_mag<0.2 in only one NIRCam band; all other NIRCam+HST+MIRI delta_mag>0.5 or undetected",
        "f444w_snr_threshold_for_approx_mag_error_lt_0p2": threshold,
        "n_candidates": len(rows),
        "n_passing_f444w_half_of_warn4": len(bright),
        "n_between_snr5_and_magerror0p2": len(marginal),
        "f444w_half_exception_ids": marginal,
        "epoch_counts": by_epoch,
        "confirmed_movers": movers,
        "negative_blue_snr_below_minus2": strong_negative,
        "n_unique_candidates_with_negative_blue_snr_below_minus2": len({x["candidate_id"] for x in strong_negative}),
        "limitations": [
            "The 124/125 result checks only the F444W magnitude-error condition; it does NOT determine actual warn_flag=4.",
            "Official COSMOS2025 photometry, other NIRCam+HST+MIRI bands, and its actual flags differ from this v13 fixture.",
            "Signed SNR<2 permits large negative deviations; flag these for background QA, not automatic deletion.",
            "Do not veto all warn_flag=4 objects in a deliberate F444W-only brown-dwarf selection.",
            "flag_star measures source-pixel overlap with JWST mask, not just source-center mask membership.",
            "Neither quality flags nor this audit measure geometric mask-WHT intersection or PM completeness."
        ]
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--fixture", default="pm_footprint/results/v13_epoch_candidate_photometry_fixture_20261031.csv")
    p.add_argument("--out", default="pm_footprint/results/dropout_quality_flag_overlap_20261101.json")
    a = p.parse_args()
    with open(a.fixture, newline="") as f:
        result = audit(csv.DictReader(f))
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"n": result["n_candidates"],
                      "f444_condition": result["n_passing_f444w_half_of_warn4"],
                      "negative_blue_outliers": result["negative_blue_snr_below_minus2"]},
                     indent=2))


if __name__ == "__main__":
    main()
