"""Package executed validation diagnostics; never produce a finite adopted count.

Run from the repository: python -m src.ali_validation_report
All expensive inputs are completed products. No new photometry is performed.
"""
import json
import shutil
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def run():
    repo = Path(__file__).resolve().parents[1]
    results = Path('/home/bahareh/results')
    out = repo / 'data/ali_validation_2026-11-28'
    out.mkdir(parents=True, exist_ok=True)
    for directory, names in [
        ('ali_domain_validation', ['domain_status.json', 'native_zero_source_selection.json',
                                   'native_radial_shell_counts.csv', 'relative_floor_convergence.json']),
        ('ali_detection_validation', ['detection_footprint_validation.json']),
    ]:
        for name in names:
            shutil.copyfile(results / directory / name, out / name)
    shutil.copyfile(results / 'ali_validation_review_20261128.md', out / 'scientific_review.md')
    d = pd.read_csv(out / 'native_radial_shell_counts.csv')
    absolute = d[d.r_high_pc.notna()]
    cumulative, bins = [], []
    for radius in sorted(absolute.r_high_pc.unique()):
        subset = absolute[absolute.r_high_pc <= radius]
        # A bin's density error is fully correlated between its radial shells.
        sigma = np.sqrt(np.sum(subset.groupby('teff_lo_K').density_error.sum()**2))
        row = dict(r_max_pc=radius, formal_measured_density_sigma=sigma,
                   lower_limit_bin_count=float(subset.loc[subset.density_status == 'lower_limit', 'total'].sum()),
                   quadrature_error_upper_bound=float(subset.quadrature_error.sum()))
        for key in ['thin', 'thick', 'halo', 'total', 'noise_floor_total', 'excess_above_no_source_floor']:
            row[key] = float(subset[key].sum())
        cumulative.append(row)
        for (lo, hi), b in subset.groupby(['teff_lo_K', 'teff_hi_K']):
            bins.append(dict(r_max_pc=radius, teff_lo_K=lo, teff_hi_K=hi,
                             total=b.total.sum(), thin=b.thin.sum(), thick=b.thick.sum(), halo=b.halo.sum(),
                             density_status=b.density_status.iloc[0],
                             formal_density_sigma=None if b.density_error.isna().all() else b.density_error.sum()))
    c = pd.DataFrame(cumulative)
    c['conditional_surface_density_arcmin2'] = c.total / (0.5324308306666868 * 3600)
    c.to_csv(out / 'cutoff_diagnostics.csv', index=False)
    pd.DataFrame(bins).to_csv(out / 'cutoff_bin_diagnostics.csv', index=False)
    r = d[d.r_high_pc.isna()].groupby(['relative_r_low', 'relative_r_high'])[
        ['thin', 'thick', 'halo', 'total', 'noise_floor_total', 'excess_above_no_source_floor']].sum().reset_index()
    r.to_csv(out / 'nominal_boundary_scatter.csv', index=False)
    footprint = json.loads((out / 'detection_footprint_validation.json').read_text())
    ledger = dict(
        historical_area_deg2=0.496306,
        historical_mask_provenance='Unrecovered: scalar benchmark alone cannot identify the spatial mask.',
        original_four_positive_wht_centers_deg2=footprint['original_positive_wht_area_deg2'],
        native_four_aperture_valid_deg2=footprint['native_four_area_deg2'],
        native_five_aperture_valid_deg2=footprint['native_five_area_deg2'],
        aperture_support_loss_deg2=footprint['original_positive_wht_area_deg2']-footprint['native_four_area_deg2'],
        f814w_coverage_loss_deg2=footprint['native_four_area_deg2']-footprint['native_five_area_deg2'],
        erroneous_detection_positive_original_four_deg2=footprint['erroneous_chi_positive_original_four_area_deg2'],
        adopted_geometric_area_deg2=footprint['native_five_area_deg2'],
        adoption_scope='Valid forced-aperture photometric measurements; excludes no-data, not all catalog masks or detection completeness.')
    (out / 'footprint_ledger.json').write_text(json.dumps(ledger, indent=2)+'\n')
    status = dict(
        base_commit='5c8cd562ef79da3ef8dc87b1b146e5c1650fd69b',
        final_count_accepted=False, final_five_band_count=None,
        final_surface_density_arcmin2=None, final_count_uncertainty=None,
        formal_all_distance_forced_photometry_integral='divergent',
        geometric_five_band_area_deg2=footprint['native_five_area_deg2'],
        native_catalog_snr_adopted=True,
        positive_auto_flux_joint_kernel_validated=False,
        historical_area_reproduced=False,
        required_dependency='Validated joint source-associated detection and photometric-selection kernel with a sufficiently vanishing zero-flux tail, or an independently justified finite population domain. Neither is supplied by the existing products.',
        prohibited_substitutions=['arbitrary radial cutoff', '50-percent boundary as hard support', 'noise-floor-subtracted excess as selected-source count'],
        no_repeated_processing=['F814W injections', 'depth reconstruction', 'atmosphere calculations', 'Galactic population parameter fitting'])
    (out / 'result_status.json').write_text(json.dumps(status, indent=2)+'\n')
    plt.rcParams.update({'font.size': 10, 'svg.fonttype': 'none', 'svg.hashsalt': 'ali-validation-20261128'})
    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    for key, label in [('total', 'Forced-photometry integral'), ('noise_floor_total', 'Zero-source floor'),
                       ('excess_above_no_source_floor', 'Floor-subtracted excess (different estimand)')]:
        ax.plot(c.r_max_pc/1000, c[key], 'o-', label=label)
    ax.set(xscale='log', xlabel='Diagnostic outer radius (kpc)', ylabel='Expected integral',
           title='No outer radius is adopted: unchanged halo has a divergent tail')
    ax.legend(fontsize=8); ax.grid(alpha=.2); fig.tight_layout()
    fig.savefig(out / 'radial_domain.svg', metadata={'Date': None}); plt.close(fig)
    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    labels = [f'{low:g}–{high:g}' for low, high in zip(r.relative_r_low, r.relative_r_high)]
    ax.bar(labels, r.noise_floor_total, label='Zero-source floor')
    ax.bar(labels, r.excess_above_no_source_floor, bottom=r.noise_floor_total, label='Flux-dependent excess')
    ax.set(xlabel='Distance / local nominal native-SNR F444W radius', ylabel='Shell integral',
           title='Noise scattering extends beyond the nominal boundary')
    ax.legend(fontsize=8); fig.tight_layout()
    fig.savefig(out / 'nominal_boundary_scatter.svg', metadata={'Date': None}); plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 4.1))
    names = ['Historical (mask unknown)', 'Four-band positive WHT', 'Four-band aperture valid',
             'Five-band aperture valid', 'Erroneous χ-positive requirement']
    values = [0.496306, footprint['original_positive_wht_area_deg2'], footprint['native_four_area_deg2'],
              footprint['native_five_area_deg2'], footprint['erroneous_chi_positive_original_four_area_deg2']]
    ax.barh(names, values, color=['gray', 'C0', 'C0', 'C2', 'C3'])
    for i, value in enumerate(values): ax.text(value+.0003, i, f'{value:.6f}', va='center', fontsize=8)
    ax.set(xlim=(.48, .542), xlabel='Area (deg²)', title='Geometric coverage and an invalid intensity mask')
    ax.invert_yaxis(); fig.tight_layout()
    fig.savefig(out / 'footprint_accounting.svg', metadata={'Date': None}); plt.close(fig)
    print(json.dumps(status, indent=2))


if __name__ == '__main__':
    run()
