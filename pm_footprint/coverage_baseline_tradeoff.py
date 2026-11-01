"""Audit WHT-threshold and time-baseline robustness using existing SMALL JSON.

No WHT image, JWST mask binary, network, or MAST query required.
This reports geometric opportunities, never PM detection completeness.
"""
import argparse
import json
from pathlib import Path


def evaluate(wht, exact):
    if wht['candidate_catalog_sha256'] != exact['candidate_catalog_sha256']:
        raise ValueError('Candidate roster SHA mismatch')
    thresholds = ('0', '0.1', '0.5')
    rows = []
    for t in thresholds:
        a = wht['area_deg2'][t]
        c = wht['candidate_counts'][t]
        base = a['archive_f444w']
        if not 0 < base < 1.0:
            raise ValueError('Implausible survey area')
        for days in (365, 730):
            red = a[f'f444w_plus_red_{days}d']
            same = a[f'same_f444w_{days}d']
            if not 0 <= same <= red <= base:
                raise ValueError('Non-nested epoch geometry')
            rows.append({
                'wht_relative_median_threshold': float(t),
                'baseline_days': days,
                'wht_area_deg2': base,
                'red_epoch_area_deg2': red,
                'red_epoch_fraction': red/base,
                'same_f444w_area_deg2': same,
                'candidate_red_count': c[f'f444w_plus_red_{days}d'],
                'candidate_same_f444w_count': c[f'same_f444w_{days}d'],
            })
        if c['f444w_plus_red_365d'] != 62 or c['f444w_plus_red_730d'] != 42:
            raise ValueError('Frozen candidate counts unexpectedly changed')
        if c['same_f444w_365d'] or c['same_f444w_730d']:
            raise ValueError('Expected zero same-band candidates at >=365d')
    baseline = next(r for r in rows if r['wht_relative_median_threshold']==0 and r['baseline_days']==365)
    strict = next(r for r in rows if r['wht_relative_median_threshold']==.5 and r['baseline_days']==365)
    two_year = next(r for r in rows if r['wht_relative_median_threshold']==0 and r['baseline_days']==730)
    ids365 = set(exact['baseline_results']['365']['candidate_ids']['f444w_plus_red'])
    ids730 = set(exact['baseline_results']['730']['candidate_ids']['f444w_plus_red'])
    if len(ids365) != 62 or len(ids730) != 42 or not ids730.issubset(ids365):
        raise ValueError('Exact candidate baseline membership inconsistent')
    return {
        'status':'GEOMETRIC_BASELINE_AND_WHT_THRESHOLD_ROBUSTNESS_NOT_PM_COMPLETENESS',
        'inputs': ['pm_footprint/results/wht60/summary.json', 'pm_footprint/results/exact/summary.json'],
        'threshold_baseline_table': rows,
        'wht_50pct_median_cut_vs_positive': {
            'survey_area_reduction_pct': 100*(1-strict['wht_area_deg2']/baseline['wht_area_deg2']),
            'red365_area_reduction_pct': 100*(1-strict['red_epoch_area_deg2']/baseline['red_epoch_area_deg2']),
            'red365_fraction_change_percentage_points':100*(strict['red_epoch_fraction']-baseline['red_epoch_fraction']),
            'candidate_count_change':strict['candidate_red_count']-baseline['candidate_red_count'],
        },
        '730d_vs_365d_positive_wht': {
            'red_area_retained_fraction': two_year['red_epoch_area_deg2']/baseline['red_epoch_area_deg2'],
            'red_area_lost_deg2': baseline['red_epoch_area_deg2']-two_year['red_epoch_area_deg2'],
            'red_candidates_retained_fraction': 42/62,
            'red_candidates_lost':20,
            'candidate_density_365_per_deg2':62/baseline['red_epoch_area_deg2'],
            'candidate_density_730_per_deg2':42/two_year['red_epoch_area_deg2'],
            'candidate_density_lost_annulus_per_deg2':20/(baseline['red_epoch_area_deg2']-two_year['red_epoch_area_deg2']),
            'density_ratio_730_to_365':(42/two_year['red_epoch_area_deg2'])/(62/baseline['red_epoch_area_deg2']),
            'confirmed_pm_ids_at_365d': exact['baseline_results']['365']['confirmed_pm_ids_in_red_mask'],
            'confirmed_pm_ids_at_730d': exact['baseline_results']['730']['confirmed_pm_ids_in_red_mask'],
        },
        'limitations':[
            'Relative WHT>0.5 median is a per-tile weight cut, NOT a calibrated F444W=26.2 selection-depth cut.',
            'WHT 0.6 arcsec raster deposition loses about 0.96 percent of per-tile area; this sensitivity check does not remove that systematic.',
            'Multi-epoch red imaging is overwhelmingly cross-filter (F356W), not proof of F356W detection or astrometric precision.',
            'Mask WCS provenance requires revalidation; no JWST mask correction is included.',
            'No source quality, local F115W/F150W/F277W errors, F356W SNR, or proper-motion completeness included.',
        ],
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--wht', default='pm_footprint/results/wht60/summary.json')
    p.add_argument('--exact', default='pm_footprint/results/exact/summary.json')
    p.add_argument('--out', default='pm_footprint/results/coverage_baseline_tradeoff_20261101.json')
    args=p.parse_args()
    result=evaluate(json.loads(Path(args.wht).read_text()),json.loads(Path(args.exact).read_text()))
    Path(args.out).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result['730d_vs_365d_positive_wht'],indent=2))


if __name__=='__main__':
    main()
