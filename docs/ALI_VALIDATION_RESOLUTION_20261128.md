# Ali: footprint, native SNR and radial-domain validation

Subsequent scope clarification: the user authorized a distinct depth-limited
population estimand. Its finite result is reported in
[ALI_DEPTH_LIMITED_20261128.md](ALI_DEPTH_LIMITED_20261128.md). This report's
all-distance divergence result remains unchanged; its acceptance ruling applies
to the unrestricted estimand considered here.

Continuation of `5c8cd562ef79da3ef8dc87b1b146e5c1650fd69b`, executed 2026-11-28.

**No finite final five-band population prediction is accepted.** The appropriate
geometric five-band measurement footprint is **0.532430831 deg²**. Adopting the
native catalog SNR rather than treating it as a correlated-noise-corrected
aperture significance exposes a nonzero zero-source selection floor. With the
unchanged, untruncated halo, the all-distance forced-photometry population
integral diverges. A finite catalog-source prediction requires an additional
validated source-associated detection/selection kernel or an independently
justified finite population domain. Neither has been identified in the existing
products. No arbitrary distance limit or population-model change was adopted.

All completed ACS calibration, injection–recovery, NIRCam depth reconstruction,
atmosphere fluxes and Galactic parameters were reused. None of those processing
steps was repeated. The new calculations inspect existing detection images,
compare cached native masks and integrate radial diagnostics with cached
probabilities. Numerical results, figures and independent review are in
[`data/ali_validation_2026-11-28`](../data/ali_validation_2026-11-28).

## 1. Footprint decision and limits of historical attribution

| Definition | Area (deg²) |
| --- | ---: |
| Historical scalar; exact spatial mask/driver unrecovered | 0.496306000 |
| Reconstructed original four-band positive-WHT pixel centers | 0.533333203 |
| Four bands with native aperture support | 0.532647461 |
| Five bands with native aperture support | **0.532430831** |
| Incorrectly requiring positive χ detection intensity, original four-band centers | 0.499244069 |

The physical, measured reductions are 0.000685742 deg² for aperture support and
0.000216630 deg² for additional ACS coverage. Together they remove 0.000902372
deg², much less than the 0.037027203 deg² historical-to-reconstructed difference.
The previously executed JWST star-mask center audit removes another 0.012578391
deg² from the five-band footprint, leaving 0.519852440 deg²; it also does not
recover the historical area. Applying extra catalog quality/bright-star flags
would change the five photometric criteria specified for this calculation.

The combined detection image sets negative band contributions to zero. Its zero
pixels can therefore be real, covered background pixels. The streamed local
detection archive test shows that requiring χ>0 removes 0.034089134 deg² from
the reconstructed original four-band footprint. This identifies a concrete
methodological way to lose a comparable area, **not the proven cause of the
historical footprint**: 0.499244069 is not 0.496306. No recovered historical
spatial mask supports attributing the difference to this operation. Area scaling
cannot reconstruct depth–mask correlations and was not used to tune the count.

For exactly the listed forced-photometry cuts, adopt valid WHT/measurement and
aperture support, including actual ACS intersection. Missing ACS measurements
are excluded rather than interpreted as non-detections. This is a geometric
measurement footprint, not a claim of complete catalog detection or source
masking. The existing one-pixel erosion audit gives a 0.00052336 deg² (0.0983%)
area sensitivity; rounded native translations and aperture-phase support remain
approximations. The catalog's source/segmentation masks cannot in general be
represented by a single intensity-based footprint cut.

Archive execution note: all 20 FITS members and packed positivity masks were
processed before a trailing non-FITS member caused a header-read error. The
completed bitmap and all 20 grid caches were verified and reused; the code now
skips non-FITS members. Unpersisted per-tile header values remain null in the
saved diagnostic, rather than being reconstructed or fabricated. The A1 pilot
header is retained. The 30 GB archive was not processed a second time.

## 2. Native catalog SNR is not the empirical aperture z-score

The primary catalog methodology distinguishes the pixel-noise-normalized SNR
maps from correlated-noise aperture uncertainty calibration; detection uses a
positive multi-band statistic and hot/cold extraction. Its F277W comparison
completeness is not a F444W-only point-source kernel. These facts supply no
hard intrinsic F444W flux boundary. See [COSMOS2025, §§3.3–3.4 and
3.6](https://arxiv.org/html/2506.03243v1#S3.SS3).

The existing paired-source calibration reconstructs the native 0.2 arcsec
diameter SNR estimator. Write it as

\[
S_b=(\widehat F_b-b_b)/s_{\mathrm{cat},b},
\]

where b is the fitted catalog offset in flux units. This offset is incorporated consistently in the mean modeled
measurement. In the count implementation, `signal + mu` is the mean measurement
before applying thresholds `k * denominator`; `denominator` is the native
catalog scale, and `noise` is the empirically measured standard deviation in
the same flux units. Thus, in that implementation's convention,

\[
P(S_{444}>5\mid F)=\Phi\left[
 (F_{444,\mathrm{ap}}+\mu_{444}-5s_{\mathrm{cat},444})/
 \sigma_{\mathrm{emp},444}\right].
\]

This is not `Phi(F / sigma_emp - 5)`. The signed upper cuts use the corresponding
lower-tail probability at `2 * s_cat`; ACS uses the already validated empirical
CDF with source photon noise. Positive F444W is automatically satisfied by its
positive five-SNR threshold in the reconstructed estimator. No absolute-SNR
cut or AUTO/ERR substitution is made. Across all tiles the F444W fitted
offset/scale lies between −0.27092 and −0.17416, so the raw flux implied by
native SNR>5 is positive.
This guarantees positivity of the reconstructed aperture measurement only.
The recovered Ali workflow specifies positive **AUTO** flux
(`IMAGE_INJECTION_ALI_DAVE_20261109.md`); that is a separate catalog measurement.
Its joint distribution with native SNR is not supplied by these forced-aperture
products. AUTO positivity must therefore also remain in the outstanding joint
catalog-selection dependency. No claim is made that native SNR alone implements
every catalog photometric condition.

Independent A1 header check: the pixel-noise normalization predicts the native
F444W aperture scale 145.1412, versus the existing paired-source fit 145.5680
(0.29% difference), in noise-equalized flux units. The predicted offset is
−30.8967 versus fitted −29.9657. Empirical aperture noise is 2.3312 times the
catalog scale in this tile. Correlated pixel noise explains why catalog SNR
must not be treated as a unit-variance aperture z-score. This single-header
cross-check does not replace the completed all-tile paired-source calibration.
The previously measured ACS estimator residual (median absolute residual
0.412–0.563 native SNR) still limits exact catalog reconstruction.

**Convention clarification:** the paired-source raw relation is
`S = (F_measured * sqrt(WHT) - B) / A`. The count model therefore stores
`s_cat=A/sqrt(WHT)` and `mu=mu_blank-B/sqrt(WHT)`. The probability expression
above uses this net `mu`; it must not subtract B a second time.

## 3. Nominal depth boundary and sources beyond it

For the Gaussian F444W estimator, the local 50% boundary is

\[
r_{50}=\left[\frac{F_{25}\,d_{25,444}^2\,\mathrm{EE}_{444}}
 {5s_{\mathrm{cat},444}-\mu_{444}}\right]^{1/2}.
\]

All 770,595 covered map cells have positive denominators (minimum 3.2861 nJy);
the defensive positivity clamp is inactive. Neither five times the empirical
noise nor five times the native denominator alone is the same boundary. More
fundamentally, the 50% radius is not a support boundary: measured-flux selection
has a soft tail. The original `count_to_distance` is an integral over an
explicitly supplied finite domain; it does not establish that domain from the
catalog rule.

Executed five-band shell sums using the native estimator in all bands:

| Distance / local r50 | Total | Zero-source floor | Flux-dependent excess |
| --- | ---: | ---: | ---: |
| 0.001–1 | 28.706317 | 3.768855 | 24.937462 |
| 1–1.5 | 15.071112 | 1.950848 | 13.120265 |
| 1.5–2 | 7.071083 | 1.906634 | 5.164448 |
| 2–5 | 16.262648 | 9.884280 | 6.378368 |
| 5–10 | 13.470587 | 12.152686 | 1.317901 |

Already 34.43% of the integral to 1.5 r50 is beyond r50; to 2 r50 this fraction
is 43.55%. The excess decomposition shows that this is not solely a zero-source
floor: 13.120265 flux-dependent excess units occur from 1 to 1.5 r50. These are
diagnostic integrals, not accepted finite survey counts or completeness
measurements. The floor-subtracted quantity is a different estimand.

## 4. Why the unrestricted prediction does not converge

Using only held-out backgrounds already measured during calibration, **526 of
84,203** common five-band forced apertures satisfy all five signed native cuts:
0.00624681, with a nominal binomial Wilson 95% interval
[0.00573670, 0.00680197]. Four-band blank passes are 629. These are reconstructed
forced-aperture estimators, not detected source entries. Spatial dependence and
calibration uncertainty are outside the binomial interval. The area-weighted
model zero-source probability is 0.00717404. Its weighting differs from the
pooled blank sample, so the two values are not a direct goodness-of-fit test.
Their nonzero nature is the relevant convergence result.

As r tends to infinity, true source flux tends to zero, the forced-measurement
selection probability approaches p0>0, and the unchanged halo density behaves
as r^(−2.77). Therefore

\[
\frac{dN_{\rm halo}}{dr}\ \propto\ r^2\rho(r)p_0
 \ \sim\ r^{-0.77},\qquad N_{\rm halo}(<R)\ \sim\ R^{0.23}.
\]

This is an analytical divergence, not failed numerical quadrature. The disks
converge, but their convergence does not remove the halo tail. The real finite
catalog does not contain infinitely many sources: the invalid extrapolation
assigns a chance noise selection separately to arbitrarily many invisible
modeled stars, without a finite detected-peak/source-association process.

The additional native-SNR calculations reuse the empirical ACS kernel:

| Diagnostic Rmax | Thin | Thick | Halo | Total | Formal measured-density σ |
| --- | ---: | ---: | ---: | ---: | ---: |
| 20 kpc | 12.1760 | 13.7935 | 35.3530 | 61.3225 | 3.2243 |
| 50 kpc | 12.1760 | 13.7935 | 56.9926 | 82.9622 | 4.1692 |
| 100 kpc | 12.1760 | 13.7935 | 76.1284 | 102.0979 | 5.0709 |
| 200 kpc | 12.1760 | 13.7935 | 98.8685 | 124.8380 | 6.1605 |
| 1 Mpc | 12.1760 | 13.7935 | 169.3608 | 195.3303 | 9.5601 |

No Rmax in this table is adopted. Per-bin and component contributions at every
diagnostic cutoff are saved in `cutoff_bin_diagnostics.csv`. Its formal density
errors sum fully correlated radial-shell contributions within each temperature
bin, then combine bins independently. Lower-limit bins have no assigned
uncertainty: their contribution is separately recorded (11.3537 objects at
20 kpc, 32.1830 at 1 Mpc). These errors are not full population uncertainties.

The 1 Mpc integral contains a 143.7747 zero-source floor. Subtracting this floor
gives an apparently converging excess of 51.5556, but **this is not an adopted
number of UCDs passing the cuts**. Subtraction changes the estimand and can
produce signed responses; convergence alone does not validate it.

A joint *source-associated* detection and photometric-selection kernel must
vanish sufficiently rapidly as flux tends to zero. If its added detection
factor behaves as F^β, convergence of this halo requires β>0.115. A generic
recovery function retaining a random-association floor is insufficient. This
is a necessary asymptotic condition, not an assumed completeness law. No such
F444W-only source-associated kernel is present in the reviewed products.
Alternatively a physically justified finite halo domain would change the
specified untruncated model and is outside this continuation's restrictions.

## 5. Comparison with previous results and uncertainty status

The archived 25.8331 approximate-transfer result and reconstructed 23.8437
empirical-ACS result use restricted conventions/domains. They cannot be compared
as competing full native-SNR predictions. Within the *archived reconstructed
convention only*, four-band original/matched counts are 28.886864/28.805101 and
the empirical five-band count is 23.843714: ACS alone removes 4.961387
(17.2240%) on that matched footprint, while footprint reduction removes
0.081763 from the four-band count. That conditional five-band value is 1.989386
(7.7009%) below 25.8331. These numbers are retained as provenance, not tuned or
newly adopted. Native-SNR hard-domain diagnostics of 12.360844 and 29.488820
likewise do not establish a unique domain. The new native r50 test in §3 also
includes the calibrated offset in its boundary.

The new quadrature error estimates sum to <4.1e−5 objects through 1 Mpc, and the
relative-shell noise-floor 64/128-node check differs by <1.2e−14 per component.
Those checks support these diagnostics, not convergence of the unrestricted
integral. Existing ACS noise/threshold/PSF/CDF sensitivity and spatial-resolution
tests remain documented in the [prior report](../data/ali_f814w_2026-11-26/REPORT.md).
They were not rerun and cannot be transferred into a finite combined error bar
for a divergent integral. No finite final surface density or total uncertainty
is reported. Local-density incompleteness, catalog-estimator residuals, noise
correlations, source association and unpropagated Galactic/atmosphere
systematics remain distinct from numerical error.

## 6. Reproduction and exact inputs

Use the unchanged server environment and run from the repository:

```sh
cd /home/bahareh/projects/BD_number_density
PYTHON=/home/bahareh/miniforge3/envs/phd/bin/python
$PYTHON -m src.ali_detection_footprint_validation
$PYTHON -m src.ali_domain_validation --cells 256 --epsrel 3e-5
$PYTHON -m src.ali_domain_validation --cells 256 --relative-floor-only
$PYTHON -m src.ali_validation_report
$PYTHON -m pytest -q tests/test_ali_domain_validation.py tests/test_ali_probabilities.py
```

The footprint command reuses its completion marker when present. For this
execution only, recovery of already completed archive products used
`--reuse-completed-grid-cache`; this is not a substitute for processing an
absent cache on a new machine. No raw image is downloaded. The reporting script
needs the saved independent review text in the server results directory; it
only packages completed calculations and figures.

Inputs are the local detection archive
`/mnt/E_202/JWST/COSMOS-Web/auxiliary/detection_images.tar.gz`, cached native masks
`/home/bahareh/results/ali_native_coverage`, original WHT masks
`/home/bahareh/results/ali_original_coverage`, native calibration
`/home/bahareh/results/ali_calibration_final`, existing held-out scene stamps
`/home/bahareh/results/ali_fiveband_reconstruction`, native spatial maps
`/home/bahareh/results/ali_final_counts/spatial_maps.npz`, and the committed
`data/elfowl_fiveband_bin_reach.csv` and fixed density functions. Full original
science/weight/PSF/catalog paths are preserved in the prior committed
`data/ali_f814w_2026-11-26/survey_input_manifest.json`. No raw FITS or packed gigabyte
masks are uploaded.

The outstanding dependency is scientific, not a missing old depth cache. The
exact historical mask/driver still has no recorded filename. It is needed only
to attribute/reproduce the historical scalar; its recovery is not a reason to
force the current area or count. Independently, the joint source-associated
detection kernel (or justified population domain) is essential to certify a
finite native-selection prediction. A targeted ACS forced-aperture experiment
cannot supply that catalog-generation ingredient.

## 7. Section 4.6 wording supported by the executed validation

A paragraph claiming a finite final updated number would not be supported.
The following is publication-ready **validation/limitation wording**, not a
replacement numerical population measurement:

> We calibrated the F814W non-detection probability using the local ACS imaging
> and the completed real-background injection–recovery experiment, retaining
> the adopted Sonora Elf Owl photometry, temperature-bin space densities and
> Galactic components. The common region with valid five-band aperture
> measurements covers 0.532431 deg². Reconstructing the native COSMOS2025 SNR
> estimator reveals a nonzero zero-source probability of satisfying the signed
> photometric cuts. With the adopted untruncated halo density proportional to
> r⁻²·⁷⁷, the corresponding all-distance forced-photometry integral does not
> converge; its halo contribution grows asymptotically as R⁰·²³. We therefore
> do not adopt a finite expected count or surface density from an arbitrary
> distance cutoff. A validated source-associated catalog detection model is
> required before assigning a finite prediction to these photometric criteria.
> Such a prediction would not describe the completeness of the final visually
> vetted sample.

Independent scientific review, executed numerical tables and the machine-readable
acceptance status accompany this report. Final numerical acceptance remains
explicitly false.
