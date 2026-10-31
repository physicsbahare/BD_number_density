"""Exploratory candidate-density contrast across the 365-day red-epoch footprint.

Uses only saved, small GitHub numerical products. This is NOT a completeness
measurement: the binomial null assumes independent uniformly selected sources
on the WHT>0 footprint, which the dropout selection does not establish.
"""
import argparse
import csv
import json
import math
from pathlib import Path


def exact_binomial_two_sided(k, n, p):
    """Probability-ordered exact two-sided test (like scipy.stats.binomtest)."""
    if not (0 <= k <= n and 0 < p < 1):
        raise ValueError("Invalid binomial parameters")
    def logpmf(i):
        return (math.lgamma(n + 1) - math.lgamma(i + 1)
                - math.lgamma(n - i + 1) + i * math.log(p)
                + (n - i) * math.log1p(-p))
    observed = logpmf(k)
    return min(1., sum(math.exp(logpmf(i)) for i in range(n + 1)
                       if logpmf(i) <= observed + 1e-10))


def positional_groups(rows, radius_arcsec=0.5):
    parent = list(range(len(rows)))
    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for i, a in enumerate(rows):
        for j in range(i + 1, len(rows)):
            b = rows[j]
            dec = math.radians((float(a['dec_deg']) + float(b['dec_deg'])) / 2)
            d = 3600 * math.hypot((float(a['ra_deg']) - float(b['ra_deg'])) * math.cos(dec),
                                  float(a['dec_deg']) - float(b['dec_deg']))
            if d < radius_arcsec:
                parent[root(j)] = root(i)
    groups = {}
    for i, row in enumerate(rows):
        groups.setdefault(root(i), []).append(row)
    return list(groups.values())


def diagnostic(rows, area_total, area_red):
    if not (0 < area_red < area_total):
        raise ValueError("Invalid areas")
    ids = [r['candidate_id'] for r in rows]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate candidate IDs")
    red = lambda r: str(r['in_f444w_plus_red_365d']).lower() == 'true'
    if any(str(r['in_archive_f444w']).lower() != 'true' for r in rows):
        raise ValueError("Some candidate centers are outside archive WHT")
    groups = positional_groups(rows)
    mixed = [g for g in groups if any(red(r) for r in g) and not all(red(r) for r in g)]
    p = area_red / area_total
    def stat(k, n):
        area_other = area_total - area_red
        return {'red': k, 'outside_red': n-k, 'total': n,
                'red_fraction': k/n, 'red_density_deg2': k/area_red,
                'outside_density_deg2': (n-k)/area_other,
                'red_to_outside_density_ratio': (k/area_red)/((n-k)/area_other),
                'uniform_independent_binomial_two_sided_p': exact_binomial_two_sided(k,n,p)}
    k = sum(red(r) for r in rows)
    kg = sum(any(red(r) for r in g) for g in groups)
    return {'area_total_deg2': area_total, 'area_red_deg2': area_red,
            'geometric_red_area_fraction': p,
            'entries': stat(k, len(rows)),
            'proximity_groups_0p5arcsec': stat(kg, len(groups)),
            'mixed_coverage_groups': [[r['candidate_id'] for r in g] for g in mixed],
            'status': 'EXPLORATORY_NOT_SELECTION_COMPLETENESS',
            'caveats': ['Uniform independent binomial null is not established: spatial selection depth and clustering vary.',
                        '0.5 arcsec proximity groups are not verified unique astrophysical objects.',
                        'Area estimates inherit ~0.96% WHT raster deposition loss.',
                        'Cross-filter F356W coverage is geometric only; no detection/PM completeness.']}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--membership', type=Path, default=Path('pm_footprint/results/wht60/dave125_wht_membership.csv'))
    parser.add_argument('--summary', type=Path, default=Path('pm_footprint/results/wht60/summary.json'))
    parser.add_argument('--out', type=Path, default=Path('pm_footprint/results/candidate_spatial_diagnostic.json'))
    args = parser.parse_args()
    with args.membership.open(newline='') as f:
        rows = list(csv.DictReader(f))
    summary = json.loads(args.summary.read_text())
    areas = summary['area_deg2']['0']
    report = diagnostic(rows, areas['archive_f444w'], areas['f444w_plus_red_365d'])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'entries': report['entries'], 'groups': report['proximity_groups_0p5arcsec']}, indent=2))


if __name__ == '__main__':
    main()
