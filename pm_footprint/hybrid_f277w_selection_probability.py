"""Analytic Gaussian probability of passing Dave's hybrid F277W-only cuts.

Conditioned on fixed F444W flux and known F277W Gaussian sigma. This is NOT
measured survey completeness, nor a PM-detection probability.
"""
import argparse
import csv
import json
import math
from pathlib import Path


def normal_cdf(z):
    return 0.5 * math.erfc(-z / math.sqrt(2.0))


def hybrid_pass_probability(f444, sigma277, true_f277=0.0):
    """Signed F277W SNR<2 AND f444/(f277 if positive else 2*sigma)>30."""
    if not all(map(math.isfinite, (f444, sigma277, true_f277))):
        raise ValueError("Inputs must be finite")
    if f444 <= 0 or sigma277 <= 0:
        raise ValueError("Flux and uncertainty must be positive")
    threshold = f444 / 30.0
    z0 = -true_f277 / sigma277
    zupper = (min(2.0*sigma277, threshold)-true_f277) / sigma277
    positive = normal_cdf(zupper)-normal_cdf(z0)
    negative = normal_cdf(z0) if 2.0*sigma277 < threshold else 0.0
    return max(0.0, min(1.0, positive+negative))


def median(values):
    v = sorted(values)
    n = len(v)
    return v[n//2] if n%2 else (v[n//2-1]+v[n//2])/2.0


def summarize(rows):
    if len(rows)!=125 or len({int(r["candidate_id"]) for r in rows})!=125:
        raise ValueError("Expected frozen 125 unique candidate rows")
    details = []
    for r in rows:
        f444 = float(r["f_F444W"])
        f277 = float(r["f_F277W"])
        sigma = float(r["new_e_F277W"])
        red = int(r["red365"])
        if red not in (0,1):
            raise ValueError("red365 must be 0 or 1")
        critical = f444/60.0
        details.append({
            "candidate_id": int(r["candidate_id"]),
            "red365": red,
            "observed_nonpositive_F277W": f277<=0,
            "negative_branch_allowed": sigma<critical,
            "sigma_over_critical": sigma/critical,
            "sigma277_njy": sigma,
            "critical_sigma277_njy": critical,
            "zero_true_f277_gaussian_pass_probability":
                hybrid_pass_probability(f444,sigma,0.0)
        })
    def stats(group):
        sub = [x for x in details if group is None or x["red365"]==group]
        return {
            "n":len(sub),
            "negative_branch_allowed_n":sum(x["negative_branch_allowed"] for x in sub),
            "negative_branch_forbidden_n":sum(not x["negative_branch_allowed"] for x in sub),
            "observed_nonpositive_F277W_n":sum(x["observed_nonpositive_F277W"] for x in sub),
            "zero_true_f277_gaussian_probability_median":
                median([x["zero_true_f277_gaussian_pass_probability"] for x in sub])
        }
    return {
        "status":"ANALYTIC_CONDITIONAL_GAUSSIAN_SELECTION_DIAGNOSTIC_NOT_COMPLETENESS",
        "paper_rule":"F277W signed SNR<2 and f444/(f277 if f277>0 else 2*sigma277)>30",
        "f444_26p2_brighter_side_limit_njy":10**((31.4-26.2)/2.5),
        "sigma277_critical_at_26p2_njy":10**((31.4-26.2)/2.5)/60,
        "zero_true_flux_discontinuity":{
            "probability_from_below":normal_cdf(2),
            "probability_at_or_from_above":normal_cdf(2)-0.5,
            "jump_probability_points":50.0
        },
        "all_candidates":stats(None),
        "red365":stats(1),
        "no_red365":stats(0),
        "within_10percent_of_critical_sigma":[
            x for x in details if abs(x["sigma_over_critical"]-1)<0.1],
        "pm_signature_candidates":[
            x for x in details if x["candidate_id"] in (387132,1651050)],
        "limitations":[
            "F444W is fixed and F277W errors are assumed Gaussian and known.",
            "The m444<26.2 cut is strict: 26.2 is the brighter-side limiting flux.",
            "Only F277W cuts are modeled, not the other bands, masks, or PM detection.",
            "The 125 preselected sources are not an unbiased survey-area sample.",
            "Probabilities are counterfactual conditional diagnostics, NOT completeness."
        ]
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--fixture",type=Path,default=Path(
        "pm_footprint/results/v13_epoch_candidate_photometry_fixture_20261031.csv"))
    p.add_argument("--out",type=Path,default=Path(
        "pm_footprint/results/hybrid_f277w_selection_probability_20261101.json"))
    a=p.parse_args()
    with a.fixture.open(newline="") as f:
        rows=list(csv.DictReader(f))
    result=summarize(rows)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("Negative branch allowed:",result["all_candidates"]["negative_branch_allowed_n"],"/125")


if __name__=="__main__":
    main()
