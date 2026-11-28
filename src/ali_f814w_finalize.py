"""Collect executed runs without promoting a failed historical validation."""
import argparse,json,shutil
from pathlib import Path
import pandas as pd

RUNS=['ali_preserved_nircam_counts','ali_preserved_convergence1024','ali_preserved_shot31','ali_preserved_noise_low','ali_preserved_noise_high','ali_preserved_threshold_low','ali_preserved_threshold_high','ali_preserved_flux_low','ali_preserved_flux_high','ali_preserved_correlation','ali_final_counts','ali_catalog_depth_counts','ali_counts_extended20000','ali_counts_extended50000','ali_legacy_baseline','ali_grid_legacy_diagnostic']


def run(a):
    a.out.mkdir(parents=True,exist_ok=True);runs={};rows=[]
    for name in RUNS:
        src=a.results/name;d=json.loads((src/'count_summary.json').read_text());runs[name]=d
        for key,v in d['results'].items():rows.append(dict(run=name,selection=key,area_deg2=d['four_area_deg2'] if key=='four_original' else d['five_area_deg2'],**v))
    pd.DataFrame(rows).to_csv(a.out/'executed_run_comparison.csv',index=False)
    (a.out/'executed_runs.json').write_text(json.dumps(runs,indent=2)+'\n')
    primary=runs['ali_preserved_nircam_counts'];r=primary['results'];n=r['five_empirical']['total'];area=primary['five_area_deg2'];n4=r['four_matched']['total']
    for key in ['ali_preserved_nircam_counts','ali_preserved_convergence1024']:
        result=runs[key]
        assert 0<result['five_area_deg2']<=result['four_area_deg2']
        table=pd.read_csv(a.results/key/'temperature_component_counts.csv')
        assert len(table)==52 and set(table.teff_lo_K)==set(range(300,2250,150))
        assert (abs(table.total-table[['thin','thick','halo']].sum(axis=1))<1e-10).all()
        pivot=table.pivot(index='teff_lo_K',columns='selection',values='total')
        assert (pivot.five_empirical<=pivot.four_matched+1e-8).all()
        assert (pivot.five_gaussian<=pivot.four_matched+1e-8).all()
    for filename in ['temperature_component_counts.csv','count_summary.json']:
        shutil.copyfile(a.results/'ali_preserved_nircam_counts'/filename,a.out/('conditional_'+filename))
    unc={}
    for nodes in [256,512]:
        p=a.results/f'ali_preserved_cdf_uncertainty{nodes}.json';d=json.loads(p.read_text());unc[str(nodes)]=d;shutil.copyfile(p,a.out/p.name)
    get=lambda name:runs[name]['results']['five_empirical']['total']
    convergence={key:runs['ali_preserved_convergence1024']['results'][key]['total']-r[key]['total'] for key in r}
    summary=dict(scientific_status='HISTORICAL_BASELINE_AND_EXACT_NATIVE_FIVE_BAND_VALIDATION_UNRESOLVED',
        final_accepted_count=None,conditional_expected_count=n,five_band_area_deg2=area,
        conditional_surface_density_arcmin2=n/(3600*area),components={k:r['five_empirical'][k] for k in ['thin','thick','halo']},
        four_band_reconstructed_full_footprint_count=r['four_original']['total'],four_band_reconstructed_full_area_deg2=primary['four_area_deg2'],four_band_matched_count=n4,
        f814w_alone_removed_count=n4-n,f814w_alone_fraction_removed=(n4-n)/n4,
        empirical_minus_gaussian_count=n-r['five_gaussian']['total'],difference_from_provisional_25_8331=n-25.8331,
        formal_measured_local_density_standard_error=r['five_empirical']['measured_density_uncertainty'],
        incomplete_local_density_bin_contribution=r['five_empirical']['incomplete_density_bin_contribution'],
        conditional_f814w_cdf_sampling_standard_error=unc['512']['conditional_cdf_sampling_standard_error'],
        linear_minus_step_ecdf_count=n-unc['512']['expected_count_ecdf_step_quadrature'],
        convergence_1024_minus_512=convergence,source_noise_31_minus_15_node_count=get('ali_preserved_shot31')-n,
        sensitivity_count_ranges=dict(noise_scale_0_95_to_1_05=[get('ali_preserved_noise_low'),get('ali_preserved_noise_high')],
            catalog_snr_threshold_shift_minus_plus_0_5=[get('ali_preserved_threshold_low'),get('ali_preserved_threshold_high')],
            f814w_flux_scale_0_95_to_1_05=[get('ali_preserved_flux_low'),get('ali_preserved_flux_high')]),
        nircam_correlation_sensitivity_count=get('ali_preserved_correlation'),
        limitations=['The reconstructed full four-band footprint is 0.533333 deg2; the archived 0.496306 deg2 mask is unrecovered.',
        'The exact historical driver/mask/selection configuration producing 26.500648 has not been identified. Missing intermediate caches were regenerated; they are no longer a blocker.',
        'The preserved-convention NIRCam calculation uses empirical-noise SNR, while stored catalog SNR uses another denominator. It is not an exact native five-band selection calculation.',
        'The default radial integration retains the existing intrinsic five-sigma F444W depth bound. Native-SNR and extended-domain diagnostic counts differ materially.',
        'ACS estimator residuals, band covariance, absolute aperture response and unpropagated population systematics prevent an unconditional combined error bar.'])
    (a.out/'result_and_validation_status.json').write_text(json.dumps(summary,indent=2)+'\n')
    bins=pd.read_csv(a.out/'conditional_temperature_component_counts.csv');b=bins[bins.selection=='five_empirical']
    table='| Temperature (K) | Thin | Thick | Halo | Total | Local-density status |\n|---|---:|---:|---:|---:|---|\n'
    for _,x in b.iterrows():table+=f"| {int(x.teff_lo_K)}–{int(x.teff_hi_K)} | {x.thin:.6f} | {x.thick:.6f} | {x.halo:.6f} | {x.total:.6f} | {x.density_status} |\n"
    report=f'''# Ali F814W real-image calculation: executed results, 2026-11-28

**The numerical reconstruction and targeted F814W tests were executed. The historical-baseline and exact-native-selection acceptance checks remain unresolved. Do not use the conditional result below as an accepted final Section 4.6 prediction.**

The conditional result preserving the reconstructed 0.15 arcsec NIRCam convention is **{n:.6f} objects** over **{area:.9f} deg²**, or **{n/(area*3600):.8f} arcmin⁻²**. The Galactic density model, thirteen 150 K bins, local densities and atmosphere normalization are preserved. All 26 archived two-band atmosphere reaches are reproduced to their recorded rounding. No survey image was downloaded.

## Footprints and the F814W effect

| Calculation | Area (deg²) | Expected count |
|---|---:|---:|
| Archived four-band benchmark, not reproduced | 0.496306 | approximately 26.500648 |
| Reconstructed four-band positive-weight footprint | {primary['four_area_deg2']:.9f} | {r['four_original']['total']:.6f} |
| Four bands on common five-band footprint | {area:.9f} | {n4:.6f} |
| Five bands, Gaussian F814W probability | {area:.9f} | {r['five_gaussian']['total']:.6f} |
| Five bands, empirical F814W probability | {area:.9f} | {n:.6f} |

The area/footprint restriction changes the four-band count by {r['four_original']['total']-n4:.6f}. F814W alone removes {n4-n:.6f} objects ({100*(n4-n)/n4:.3f}%) on the matched footprint. Replacing the Gaussian optical probability by the empirical probability changes the count by {n-r['five_gaussian']['total']:.6f}. The conditional count differs from the provisional 25.8331 by {n-25.8331:.6f} ({100*(n/25.8331-1):.3f}%); this comparison includes a changed footprint and is not an isolated F814W correction to the archived run.

`four_original` in the raw execution tables means the reconstructed native positive-weight footprint. It does **not** identify the missing archived 0.496306 deg² spatial mask. Scalar area-normalizing 28.886864 to the archived area gives 26.881364, 1.44% above the benchmark, but cannot reproduce spatial depth/mask correlations and is not accepted as a baseline pass.

## Temperature and Galactic components

{table}
Thin/thick/halo totals are {r['five_empirical']['thin']:.6f}, {r['five_empirical']['thick']:.6f}, and {r['five_empirical']['halo']:.6f}. Full signed-cut, matched-footprint and density-error tables are in `conditional_temperature_component_counts.csv`.

## Observations and targeted validation

The twenty local 30 mas ACS SCI/ERR/WHT tiles are in `/mnt/E_202/JWST/COSMOS_Web/HST_F814W_v0.5.2026jul/tiles`; the four NIRCam bands are in `/mnt/E_202/JWST/COSMOS-Web/NIRCam_mosaics`. NIRCam MJy/sr values use PIXAR_SR to obtain nJy/pixel; ACS SCI and ERR are nJy. WHT units and the relation to ERR were inspected. Reciprocal total ERR variance is not assumed. Eight empirical WHT strata per tile calibrate blank-aperture noise and background, including correlated noise.

ACS native SNR is reconstructed from a fractional 0.10 arcsec-radius noise-equalized aperture, with a paired-source fitted scale/offset and static registration. It is not AUTO flux/error. Held-out ACS median absolute residuals are 0.412–0.563 SNR, requiring a systematic sensitivity rather than an exact-estimator claim. The full calibration and input inventory are saved beside this report. Missing ACS coverage is excluded, not counted as a nondetection. All nondetection inequalities are signed.

Targeted tests injected 851,643 sources at nine flux levels into 94,627 distinct real background scenes, including near-zero, subthreshold, near-threshold and brighter controls. Another population-flux test uses the adopted thirteen-bin atmosphere fluxes and local NIRCam depth at true F444W SNR 5, 8 and 15. These are forced-aperture tests, not detection, blending, morphology, quality-flag or visual-selection completeness.

The observed blank-source pass fraction is 0.832173 versus Gaussian 0.845457. This discrepancy motivates the empirical CDF. The CDF uses 47,318 held-out scenes; its independent injection check uses the other 47,309. Source photon noise is included in recovery and convolved with the CDF. The 2σ injected-source test gives 0.178740 recovered versus 0.173521 predicted: approximately 3 binomial standard errors before accounting for CDF-estimation uncertainty, so agreement is approximate rather than exact. Other tested levels and depth strata are tabulated. The 118 unmasked stellar controls give observed/model core aperture-ratio p16/median/p84 = 0.97569/1.01208/1.04531; this validates a core ratio, not the absolute wing correction.

## Uncertainties and convergence

Available independent-bin local-density errors give ±{r['five_empirical']['measured_density_uncertainty']:.6f} objects. Lower-limit bins contribute {r['five_empirical']['incomplete_density_bin_contribution']:.6f} and have no propagated formal density error; this is not zero physical uncertainty. Density covariance is unavailable. Finite-CDF sampling gives a conditional standard error of {unc['512']['conditional_cdf_sampling_standard_error']:.6f}, with shared scene dependence across temperature bins retained. The step-CDF uncertainty estimator differs from the interpolated-CDF main count by {summary['linear_minus_step_ecdf_count']:.6f}; record this approximation separately. Its radial 256/512-node central values differ by {unc['512']['expected_count_ecdf_step_quadrature']-unc['256']['expected_count_ecdf_step_quadrature']:.9f}.

Doubling spatial leaves changes the five-band count by {convergence['five_empirical']:.6f}; 15→31 source-noise quadrature nodes changes it by {summary['source_noise_31_minus_15_node_count']:.9f}. Adaptive quadrature's reported bound is {r['five_empirical']['quadrature_error_bound']:.3g}, which does not include spatial compression, calibration, footprint or model uncertainty.

F814W empirical noise ±5% gives {get('ali_preserved_noise_low'):.6f} / {get('ali_preserved_noise_high'):.6f}; a ±0.5 native-SNR threshold perturbation gives {get('ali_preserved_threshold_low'):.6f} / {get('ali_preserved_threshold_high'):.6f}; aperture-flux response ±5% gives {get('ali_preserved_flux_low'):.6f} / {get('ali_preserved_flux_high'):.6f}. These are sensitivities, not formal confidence intervals. A Gaussian F277W/F444W correlation sensitivity at the measured r015 value 0.224419 gives {get('ali_preserved_correlation'):.6f}; the other NIRCam probabilities remain in the recovered framework. Absolute PSF response, source-noise gain approximation, catalog-estimator residuals, unresolved population-model systematics and possible background correlations remain outside a combined error bar.

One extra native-pixel footprint erosion bounds phase-support area sensitivity at 0.00052336 deg² (0.0983%). Rounded native WCS translations and optional SCI/ERR masks in older validity caches are additional documented footprint approximations. A 3 arcsec-center audit of existing JWST star masks removes about 0.0125784 deg², leaving 0.5198524 deg²; it does not recover the archived area. No star-mask or catalog quality cut is silently appended to the requested five photometric criteria.

## Required scientific acceptance check

The missing input is the **exact historical Ali spatial driver plus footprint/mask and SNR/distance-domain configuration** that produced 26.500648 over 0.496306 deg². Its filename is not recorded in the recovered repository, so no fictitious exact path is supplied. The old intermediate depth files are no longer the obstacle: they have been regenerated.

The preserved NIRCam convention uses empirical-noise SNR, whereas actual catalog SNR has a different denominator. A diagnostic using the native estimators in all bands gives 12.360844 on the retained empirical five-sigma depth domain, versus 29.488820 if that intrinsic boundary is instead set by native catalog SNR. Extending the same native-Gaussian model to 20/50 kpc gives 62.187447/84.149728. These are domain diagnostics, not alternative adopted predictions. They establish that choosing an undocumented domain or replacing the NIRCam selection convention would change the scientific estimand. No Galactic parameters are changed to resolve this.

Required action is to recover the historical driver/configuration or have its author reconcile these footprint and SNR/domain definitions against the documented calculations. Until then, an exact five-band photometric prediction and an accepted publication paragraph cannot be certified. This is a scientific validation failure, not a missing-cache execution failure.

## Section 4.6 wording: conditional draft only

> We evaluated the F814W nondetection probability using the local ACS mosaics, empirically calibrated blank-aperture noise and spatial PSF models, while retaining the recovered Kirkpatrick et al. (2024) temperature-bin densities, Sonora Elf Owl flux normalization and Galactic components. Within the reconstructed NIRCam selection convention and its intrinsic five-sigma F444W distance domain, the common five-band footprint is {area:.6f} deg² and the predicted count is {n:.2f}, corresponding to {n/(area*3600):.5f} arcmin⁻²; the four-band prediction on the same footprint is {n4:.2f}. Targeted real-background injections motivate an empirical F814W probability with source photon noise. This conditional calculation concerns photometric criteria, not the final visually vetted sample. The historical footprint and native-SNR conventions remain to be reconciled before adopting the prediction.

This draft is honest English method/result text, but is **not** a certified final paragraph for the paper. Reproduction commands are in `docs/ALI_F814W_REPRODUCE_20261128.md`; machine-readable status is `result_and_validation_status.json`.
'''
    (a.out/'REPORT.md').write_text(report)
    print(json.dumps({k:v for k,v in summary.items() if k not in ['limitations','sensitivity_count_ranges']},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--results',type=Path,default=Path('/home/bahareh/results'));p.add_argument('--out',type=Path,required=True);run(p.parse_args())
