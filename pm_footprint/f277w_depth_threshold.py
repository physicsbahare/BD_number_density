"""Aperture-noise-equivalent F277W depth for Dave's hybrid flux-ratio selection.

Only the frozen 125-source fixture is used. Not a sky-area completeness estimate.
"""
import argparse
import csv
import json
import math
from pathlib import Path


def m5_from_sigma(sigma_njy):
    if not math.isfinite(sigma_njy) or sigma_njy <= 0:
        raise ValueError("positive finite uncertainty required")
    return 31.4 - 2.5 * math.log10(5 * sigma_njy)


def required_m5(m444):
    if not math.isfinite(m444):
        raise ValueError("finite magnitude required")
    return m444 + 2.5 * math.log10(12.0)


def median(values):
    v = sorted(values)
    n = len(v)
    return v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2


def summarize(rows):
    if len(rows) != 125 or len({int(x["candidate_id"]) for x in rows}) != 125:
        raise ValueError("Expected frozen 125 unique candidates")
    details = []
    for row in rows:
        f444 = float(row["f_F444W"])
        sigma = float(row["new_e_F277W"])
        red = int(row["red365"])
        if not (math.isfinite(f444) and f444 > 0 and red in (0, 1)):
            raise ValueError("Invalid flux or epoch flag")
        m444 = float(row["mag_F444W"])
        m5 = m5_from_sigma(sigma)
        ratio = f444 / (2 * sigma)
        robust = ratio > 30
        if robust != (m5 > required_m5(m444)):
            raise AssertionError("Magnitude/flux threshold equivalence failed")
        details.append(dict(id=int(row["candidate_id"]), red365=red,
                            m5_F277W_aperture_equivalent=m5,
                            ratio_using_2sigma=ratio,
                            conservative_ratio_pass=robust,
                            observed_nonpositive_F277W=float(row["f_F277W"]) <= 0))
    def stats(data):
        ds = [x["m5_F277W_aperture_equivalent"] for x in data]
        return dict(n=len(data), median_m5_F277W_aperture_equivalent=median(ds),
                    range_m5_F277W_aperture_equivalent=[min(ds), max(ds)],
                    n_m5_at_least_28p9=sum(x >= 28.9 for x in ds),
                    n_conservative_ratio_pass=sum(x["conservative_ratio_pass"] for x in data),
                    n_nonpositive_observed_F277W=sum(x["observed_nonpositive_F277W"] for x in data))
    return dict(status="FROZEN_SOURCE_APERTURE_DEPTH_DIAGNOSTIC_NOT_COMPLETENESS",
                m444_reference_AB=26.2,
                required_F277W_5sigma_AB_at_m444_26p2=required_m5(26.2),
                nominal_F277W_5sigma_AB=28.9,
                nominal_margin_mag=28.9-required_m5(26.2),
                derivation="F444/(2sigma277)>30 iff m5_277 > m444 + 2.5log10(60/5) = m444 + 2.5log10(12)",
                all=stats(details),
                red365=stats([x for x in details if x["red365"]]),
                no_red365=stats([x for x in details if not x["red365"]]),
                pm_signature_candidates=[x for x in details if x["id"] in (387132, 1651050)],
                limitations=[
                    "5sigma equivalent uses each selected source new_e_F277W aperture error, not an independent local depth map or injection recovery.",
                    "The 125 preselected sources cannot estimate the sky area passing a selection threshold.",
                    "Positive measured F277W flux may pass the hybrid ratio when its 2sigma upper limit cannot certify ratio>30.",
                    "Nominal 28.9 is a near-threshold reference; spatial noise, aperture, mask, and F444W errors remain unmodeled.",
                    "Neither geometric multi-epoch overlap nor F356W astrometric completeness is measured here."
                ])


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--fixture", type=Path, default=Path("pm_footprint/results/v13_epoch_candidate_photometry_fixture_20261031.csv"))
    p.add_argument("--out", type=Path, default=Path("pm_footprint/results/f277w_depth_threshold_20261101.json"))
    a = p.parse_args()
    with a.fixture.open(newline="") as f:
        result = summarize(list(csv.DictReader(f)))
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
