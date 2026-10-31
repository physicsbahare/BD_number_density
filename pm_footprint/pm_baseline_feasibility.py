"""Conditional 3-sigma proper-motion feasibility using F444W/F356W time baselines.

This is a *time-baseline-only* diagnostic. Positional errors are assumptions,
not measured F356W detections, and outputs are NOT PM completeness or area.

Run from repo root:
 python -m pm_footprint.pm_baseline_feasibility \
  --input pm_footprint/results/f356w_baselines/candidate_f356w_baselines.csv \
  --output pm_footprint/results/f356w_baselines/pm_feasibility_scenarios.json
"""
import argparse
import csv
import json
import math
from pathlib import Path

YEAR_DAYS = 365.25
SIGMA_LEVEL = 3.0
SIGMAS_MAS = (10, 15, 20, 30)
MUS_MASYR = (20, 30, 50, 80, 100)


def mu_threshold(gap_days, sigma_per_epoch_mas, significance=SIGMA_LEVEL):
    """Minimum |mu| at significance for equal independent epoch errors."""
    if gap_days <= 0 or sigma_per_epoch_mas <= 0 or significance <= 0:
        raise ValueError("gap, sigma, significance must be positive")
    return significance * math.sqrt(2) * sigma_per_epoch_mas * YEAR_DAYS / gap_days


def sigma_requirement(gap_days, mu_masyr, significance=SIGMA_LEVEL):
    """Largest allowed equal per-epoch 1-D error for target |mu|."""
    if gap_days <= 0 or mu_masyr <= 0 or significance <= 0:
        raise ValueError("gap, mu, significance must be positive")
    return mu_masyr * gap_days / (YEAR_DAYS * significance * math.sqrt(2))


def analyze(rows):
    if len(rows) != 125 or len({r['candidate_id'] for r in rows}) != 125:
        raise ValueError("Requires exactly 125 distinct frozen candidate IDs")
    eligible = []
    for r in rows:
        g = r.get('best_f356w_gap_days', '').strip()
        flag = r.get('f356w_ge365', '').strip().lower()
        if flag not in ('true', 'false'):
            raise ValueError("Invalid f356w_ge365 flag")
        if flag == 'true':
            if not g or not math.isfinite(float(g)) or float(g) < 365:
                raise ValueError("Inconsistent >=365d F356W baseline")
            eligible.append((int(r['candidate_id']), float(g)))
    matrix = {}
    for sigma in SIGMAS_MAS:
        matrix[str(sigma)] = {str(mu): sum(
            mu >= mu_threshold(g, sigma) for _, g in eligible)
            for mu in MUS_MASYR}
    ranked = [{
        'candidate_id': cid,
        'max_available_conservative_gap_days': round(g, 6),
        'min_mu_3sigma_if_sigma20_masyr': round(mu_threshold(g, 20), 3),
        'max_sigma_per_epoch_3sigma_if_mu50_mas': round(sigma_requirement(g, 50), 3)
    } for cid, g in sorted(eligible, key=lambda v: (-v[1], v[0]))]
    return {
        'status': 'CONDITIONAL_BASELINE_ONLY_FEASIBILITY_NOT_MEASURED_PM_COMPLETENESS',
        'source': 'pm_footprint/results/f356w_baselines/candidate_f356w_baselines.csv',
        'n_candidates': len(rows),
        'n_geometric_f356w_ge365': len(eligible),
        'significance_sigma': SIGMA_LEVEL,
        'assumed_independent_equal_epoch_positional_errors_mas': list(SIGMAS_MAS),
        'assumed_true_proper_motions_mas_per_year': list(MUS_MASYR),
        'conditional_baseline_sufficient_counts_by_sigma_then_mu': matrix,
        'example_median_gap_days': 752.1011143171309,
        'example_median_gap_required_equal_sigma_per_epoch_for_mu50_3sigma_mas': sigma_requirement(752.1011143171309, 50),
        'candidate_priority_by_longest_baseline': ranked,
        'limitations': [
            'No measured source flux or F356W detection: all positional uncertainties and true proper motions are hypothetical.',
            'Only the maximum available gap is considered; the pair with maximum gap may have worse SNR or registration.',
            'Cross-filter F444W/F356W chromatic astrometric biases, morphology, confusion, correlations and registration systematics are omitted.',
            '3sigma is a scalar displacement approximation with equal independent 1D epoch errors, not a calibrated two-dimensional PM test.',
            'Counts are conditional among 62 geometrically covered candidates, not 125-source PM completeness and NOT survey effective area.'
        ]
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input', default='pm_footprint/results/f356w_baselines/candidate_f356w_baselines.csv')
    p.add_argument('--output', default='pm_footprint/results/f356w_baselines/pm_feasibility_scenarios.json')
    a = p.parse_args()
    with open(a.input, newline='') as f:
        result = analyze(list(csv.DictReader(f)))
    Path(a.output).parent.mkdir(parents=True, exist_ok=True)
    Path(a.output).write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'n': result['n_geometric_f356w_ge365'], 'matrix': result['conditional_baseline_sufficient_counts_by_sigma_then_mu']}, indent=2))


if __name__ == '__main__':
    main()
