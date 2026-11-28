# Ali F814W real-image calculation: executed results, 2026-11-28

**The numerical reconstruction and targeted F814W tests were executed. The historical-baseline and exact-native-selection acceptance checks remain unresolved. Do not use the conditional result below as an accepted final Section 4.6 prediction.**

The conditional result preserving the reconstructed 0.15 arcsec NIRCam convention is **23.843714 objects** over **0.532430831 deg²**, or **0.01243965 arcmin⁻²**. The Galactic density model, thirteen 150 K bins, local densities and atmosphere normalization are preserved. All 26 archived two-band atmosphere reaches are reproduced to their recorded rounding. No survey image was downloaded.

## Footprints and the F814W effect

| Calculation | Area (deg²) | Expected count |
|---|---:|---:|
| Archived four-band benchmark, not reproduced | 0.496306 | approximately 26.500648 |
| Reconstructed four-band positive-weight footprint | 0.533333203 | 28.886864 |
| Four bands on common five-band footprint | 0.532430831 | 28.805101 |
| Five bands, Gaussian F814W probability | 0.532430831 | 24.165091 |
| Five bands, empirical F814W probability | 0.532430831 | 23.843714 |

The area/footprint restriction changes the four-band count by 0.081763. F814W alone removes 4.961387 objects (17.224%) on the matched footprint. Replacing the Gaussian optical probability by the empirical probability changes the count by -0.321377. The conditional count differs from the provisional 25.8331 by -1.989386 (-7.701%); this comparison includes a changed footprint and is not an isolated F814W correction to the archived run.

`four_original` in the raw execution tables means the reconstructed native positive-weight footprint. It does **not** identify the missing archived 0.496306 deg² spatial mask. Scalar area-normalizing 28.886864 to the archived area gives 26.881364, 1.44% above the benchmark, but cannot reproduce spatial depth/mask correlations and is not accepted as a baseline pass.

## Temperature and Galactic components

| Temperature (K) | Thin | Thick | Halo | Total | Local-density status |
|---|---:|---:|---:|---:|---|
| 300–450 | 5.494112 | 0.395805 | 0.069411 | 5.959328 | lower_limit |
| 450–600 | 8.155405 | 1.781872 | 0.419529 | 10.356806 | measurement |
| 600–750 | 2.345123 | 1.816249 | 0.589559 | 4.750931 | measurement |
| 750–900 | 0.248197 | 0.896872 | 0.425346 | 1.570415 | measurement |
| 900–1050 | 0.019929 | 0.380190 | 0.271646 | 0.671765 | measurement |
| 1050–1200 | 0.000861 | 0.132091 | 0.159049 | 0.292001 | measurement |
| 1200–1350 | 0.000056 | 0.067180 | 0.136461 | 0.203697 | measurement |
| 1350–1500 | 0.000000 | 0.005953 | 0.022999 | 0.028953 | measurement |
| 1500–1650 | 0.000000 | 0.000988 | 0.007831 | 0.008819 | measurement |
| 1650–1800 | 0.000000 | 0.000053 | 0.000882 | 0.000935 | measurement |
| 1800–1950 | 0.000000 | 0.000001 | 0.000060 | 0.000061 | measurement |
| 1950–2100 | 0.000000 | 0.000000 | 0.000003 | 0.000003 | measurement |
| 2100–2250 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | lower_limit |

Thin/thick/halo totals are 16.263684, 5.477255, and 2.102775. Full signed-cut, matched-footprint and density-error tables are in `conditional_temperature_component_counts.csv`.

## Observations and targeted validation

The twenty local 30 mas ACS SCI/ERR/WHT tiles are in `/mnt/E_202/JWST/COSMOS_Web/HST_F814W_v0.5.2026jul/tiles`; the four NIRCam bands are in `/mnt/E_202/JWST/COSMOS-Web/NIRCam_mosaics`. NIRCam MJy/sr values use PIXAR_SR to obtain nJy/pixel; ACS SCI and ERR are nJy. WHT units and the relation to ERR were inspected. Reciprocal total ERR variance is not assumed. Eight empirical WHT strata per tile calibrate blank-aperture noise and background, including correlated noise.

ACS native SNR is reconstructed from a fractional 0.10 arcsec-radius noise-equalized aperture, with a paired-source fitted scale/offset and static registration. It is not AUTO flux/error. Held-out ACS median absolute residuals are 0.412–0.563 SNR, requiring a systematic sensitivity rather than an exact-estimator claim. The full calibration and input inventory are saved beside this report. Missing ACS coverage is excluded, not counted as a nondetection. All nondetection inequalities are signed.

Targeted tests injected 851,643 sources at nine flux levels into 94,627 distinct real background scenes, including near-zero, subthreshold, near-threshold and brighter controls. Another population-flux test uses the adopted thirteen-bin atmosphere fluxes and local NIRCam depth at true F444W SNR 5, 8 and 15. These are forced-aperture tests, not detection, blending, morphology, quality-flag or visual-selection completeness.

The observed blank-source pass fraction is 0.832173 versus Gaussian 0.845457. This discrepancy motivates the empirical CDF. The CDF uses 47,318 held-out scenes; its independent injection check uses the other 47,309. Source photon noise is included in recovery and convolved with the CDF. The 2σ injected-source test gives 0.178740 recovered versus 0.173521 predicted: approximately 3 binomial standard errors before accounting for CDF-estimation uncertainty, so agreement is approximate rather than exact. Other tested levels and depth strata are tabulated. The 118 unmasked stellar controls give observed/model core aperture-ratio p16/median/p84 = 0.97569/1.01208/1.04531; this validates a core ratio, not the absolute wing correction.

## Uncertainties and convergence

Available independent-bin local-density errors give ±1.768642 objects. Lower-limit bins contribute 5.959328 and have no propagated formal density error; this is not zero physical uncertainty. Density covariance is unavailable. Finite-CDF sampling gives a conditional standard error of 0.049560, with shared scene dependence across temperature bins retained. The step-CDF uncertainty estimator differs from the interpolated-CDF main count by 0.024391; record this approximation separately. Its radial 256/512-node central values differ by -0.000000693.

Doubling spatial leaves changes the five-band count by 0.014101; 15→31 source-noise quadrature nodes changes it by -0.000006004. Adaptive quadrature's reported bound is 1.73e-06, which does not include spatial compression, calibration, footprint or model uncertainty.

F814W empirical noise ±5% gives 24.173685 / 23.540840; a ±0.5 native-SNR threshold perturbation gives 21.872094 / 25.355248; aperture-flux response ±5% gives 23.854183 / 23.833204. These are sensitivities, not formal confidence intervals. A Gaussian F277W/F444W correlation sensitivity at the measured r015 value 0.224419 gives 23.748729; the other NIRCam probabilities remain in the recovered framework. Absolute PSF response, source-noise gain approximation, catalog-estimator residuals, unresolved population-model systematics and possible background correlations remain outside a combined error bar.

One extra native-pixel footprint erosion bounds phase-support area sensitivity at 0.00052336 deg² (0.0983%). Rounded native WCS translations and optional SCI/ERR masks in older validity caches are additional documented footprint approximations. A 3 arcsec-center audit of existing JWST star masks removes about 0.0125784 deg², leaving 0.5198524 deg²; it does not recover the archived area. No star-mask or catalog quality cut is silently appended to the requested five photometric criteria.

## Required scientific acceptance check

The missing input is the **exact historical Ali spatial driver plus footprint/mask and SNR/distance-domain configuration** that produced 26.500648 over 0.496306 deg². Its filename is not recorded in the recovered repository, so no fictitious exact path is supplied. The old intermediate depth files are no longer the obstacle: they have been regenerated.

The preserved NIRCam convention uses empirical-noise SNR, whereas actual catalog SNR has a different denominator. A diagnostic using the native estimators in all bands gives 12.360844 on the retained empirical five-sigma depth domain, versus 29.488820 if that intrinsic boundary is instead set by native catalog SNR. Extending the same native-Gaussian model to 20/50 kpc gives 62.187447/84.149728. These are domain diagnostics, not alternative adopted predictions. They establish that choosing an undocumented domain or replacing the NIRCam selection convention would change the scientific estimand. No Galactic parameters are changed to resolve this.

Required action is to recover the historical driver/configuration or have its author reconcile these footprint and SNR/domain definitions against the documented calculations. Until then, an exact five-band photometric prediction and an accepted publication paragraph cannot be certified. This is a scientific validation failure, not a missing-cache execution failure.

## Section 4.6 wording: conditional draft only

> We evaluated the F814W nondetection probability using the local ACS mosaics, empirically calibrated blank-aperture noise and spatial PSF models, while retaining the recovered Kirkpatrick et al. (2024) temperature-bin densities, Sonora Elf Owl flux normalization and Galactic components. Within the reconstructed NIRCam selection convention and its intrinsic five-sigma F444W distance domain, the common five-band footprint is 0.532431 deg² and the predicted count is 23.84, corresponding to 0.01244 arcmin⁻²; the four-band prediction on the same footprint is 28.81. Targeted real-background injections motivate an empirical F814W probability with source photon noise. This conditional calculation concerns photometric criteria, not the final visually vetted sample. The historical footprint and native-SNR conventions remain to be reconciled before adopting the prediction.

This draft is honest English method/result text, but is **not** a certified final paragraph for the paper. Reproduction commands are in `docs/ALI_F814W_REPRODUCE_20261128.md`; machine-readable status is `result_and_validation_status.json`.
