# Ali Section 4.6: adopted depth-limited foreground prediction

Continuation of d431c5e, 2026-11-28. This report adopts the finite population
definition explicitly authorized in the follow-up request. It supersedes the
previous acceptance status **for this depth-limited estimand only**. The earlier
all-distance divergence diagnosis remains valid for a different question.

## Definition and literature basis

The prediction is the expected number of photometric passes **among Galactic
UCDs intrinsically brighter than the local empirical five-sigma F444W
point-source depth**. For temperature bin T and sky position x,

\[
r_{\max}(T,x)=d_{25,444}(T)
 \sqrt{F_{25}\,\mathrm{EE}_{444}(x)/(5\sigma_{444,\mathrm{emp}}(x))},
\qquad F_{25}=10^{(31.4-25)/2.5}\ {\rm nJy},
\]

\[
N_{5,\mathrm{depth}}=\sum_T\int_{\Omega_5}d\Omega
 \int_0^{r_{\max}(T,x)}r^2\rho_T(r)
 P_{\mathrm{native},4}(T,r,x)P_{814,\mathrm{emp}}(T,r,x)\,dr.
\]

The two probabilities retain the existing band-independence framework and
ACS source-photon-noise treatment; a measured NIRCam covariance sensitivity is
reported separately. The boundary uses empirical background noise and the
existing PSF aperture response. It does not use the catalog's smaller native
SNR denominator as if it were correlated-noise-corrected noise. Net background
offsets enter measured-SNR probabilities, not true stellar luminosities.

Apparent-magnitude domains and sensitivity-derived effective volumes have
precedents in [Ryan & Reid (2016)](https://arxiv.org/pdf/1510.05019),
[Aganze et al. (2022)](https://escholarship.org/content/qt76c609m6/qt76c609m6.pdf),
[Chen et al. (2025)](https://arxiv.org/html/2503.18510v1), and
[FC-ENZO](https://arxiv.org/html/2607.02354v1). Their implementations differ:
the [focused literature review](ALI_DEPTH_LIMITED_LITERATURE_20261128.md)
records those differences, including Aganze's scatter extension and Ryan &
Reid's magnitude-count integral. They support this type of estimand, not a
unique universal five-sigma choice. No cited paper's distance, atmosphere model
or Galactic parameters were imported.

There is **no fixed 20/50 kpc cap**. The old execution JSON contains an unused
`rmax_pc=50000` argument: the `depth` branch integrates 0–1 times each local
luminosity/sensitivity reach and never reads that parameter. Positive-WHT
coverage and native aperture support give the adopted common area
0.532430830667 deg². The historical 0.496306 deg² spatial mask remains
unidentified and is not imposed. Missing ACS coverage is excluded.

## Adopted result and interpretation

The existing native-estimator, 512-leaf calculation on exactly this domain is
reused: **12.360844 objects**, or **0.006448852 arcmin⁻²**, over
**0.532430831 deg²**. Contributions are **8.335186 thin disk**, **2.972338 thick
disk**, and **1.053320 halo**. These are model expectations, not measured
population counts. All thirteen 150 K bins from 300 to 2250 K, their component
contributions and effective selected volumes N/rho_local are saved in
`data/ali_depth_limited_2026-11-28/five_band_bins_and_effective_volumes.csv`.

The implemented photometric probabilities use native SNR_F444W>5 and signed
SNR<2 in F115W, F150W, F277W and F814W. Positive reconstructed F444W aperture
flux follows from native SNR>5 in this calibration. Positive catalog AUTO flux
is a separate estimator whose joint distribution is not simulated. Thus this
is a **depth-limited photometric foreground approximation**, not an exact
prediction of all catalog entries passing every extraction/flux/quality rule.
No detection, deblending, morphology, bright-star masking, AUTO-positivity,
quality-flag or visual-vetting completeness is claimed.

This definition deliberately excludes sources intrinsically below the adopted
depth even if noise could scatter them into the cuts. That is a restriction of
the population, not an assertion that their selection probability is zero.
It does not reconstruct the unrestricted or exact catalog-selection count.

## Matched-footprint comparison and provenance

| Prediction | Area (deg²) | Expected count |
| --- | ---: | ---: |
| Four bands on original reconstructed positive-WHT footprint | 0.533333203 | 14.908490 |
| Four bands on five-band footprint | 0.532430831 | 14.861436 |
| Five bands, Gaussian ACS diagnostic | 0.532430831 | 12.531692 |
| **Five bands, adopted empirical ACS** | **0.532430831** | **12.360844** |

The original four-band area here is the already calculated pixel-center WHT
definition; the aperture-valid four-band area separately audited is
0.532647461 deg². These definitions are not silently conflated. The ACS-only
effect is evaluated on the identical five-band measurement footprint and the
same intrinsic F444W domain: **2.500592 objects removed**, or **16.8260%**.
The original-to-matched area/depth-mask change removes 0.047054 four-band objects.

The previous conditional 23.843714 used the recovered r=0.15 arcsec NIRCam
aperture calibration and empirical-noise SNR convention. The adopted calculation
uses native r=0.10 arcsec estimators and associated depth/PSF response, with
native SNR thresholds and calibrated net offsets. The change is −11.482870
(−48.1589%). It is **not** a new ACS-only suppression or a change to the Galactic
population. The old count is not tuned or reproduced as a target. Comparisons
to 25.8331 are also convention-dependent; that approximate-transfer result is
not a full native-estimator benchmark.

## Uncertainty and depth-boundary sensitivity

The central formal measured-density error is **0.841443 objects**, conditional
on the fixed model and measured-bin uncertainties. The lower-limit density
bins contribute **5.036239 objects**, almost entirely from 300–450 K
(the 2100–2250 K contribution is <3e−11). They have no propagated density uncertainty.
These incomplete bins prevent interpreting a symmetric error bar as a full
population confidence interval. The finite-CDF error and calibration/domain
sensitivities are tabulated in the executed result package.
The additional linear-weight propagation through all temperature bins gives
an approximate conditional ACS CDF sampling error of **0.025519 objects**.
Combining this with measured-density errors in quadrature gives **0.841830**,
still a partial uncertainty. Its step-CDF diagnostic count is 12.348251,
0.012593 (0.102%) below the adopted interpolated-CDF count. The uncertainty
calculation uses 256 radial nodes and preserves shared-CDF correlations between
temperature bins; it approximates different background scenes as independent.
The step/interpolated difference is a separate CDF implementation diagnostic,
not a measured calibration error or a replacement central result.

F814W empirical-noise scaling by 0.95/1.05 gives 12.532729/12.205160;
shifting its native-SNR cut to 1.5/2.5 gives 11.359221/13.125418.
These calibration stress tests are sensitivities, not assigned one-sigma errors.

New spatial grouping at 1024 versus 512 leaves per tile gives **12.369231**,
a change of **0.008387 (0.0679%)**. The central adaptive radial quadrature error
estimate is <4.3e−7 objects. These numerical errors are distinct from observational
calibration and model uncertainties.

Changing only the intrinsic F444W population boundary to four empirical sigma
(0.2423 mag deeper, distance ×1.1180) gives **16.085977**; six sigma
(0.1980 mag shallower, distance ×0.9129) gives **9.896737**. The measured SNR
selection remains unchanged. These are alternative population definitions,
**not confidence limits**. The alternative native mean-SNR=5 boundary is also
saved as a convention sensitivity: **28.722682** objects, not adopted in place
of the calibrated five-empirical-sigma depth. The native estimator's smaller
denominator admits an intrinsically fainter population than the five-empirical-
sigma definition. This large difference makes the estimand label essential;
it is not an uncertainty interval around 12.360844.

The five-sigma sensitivity-derived median reach varies across the existing
temperature bins from 0.735 kpc (300–450 K) to 11.331 kpc (2100–2250 K).
The spatial 16th/50th/84th percentiles for every bin and boundary are saved
alongside the sensitivities. These are sky-depth distributions, not distance
uncertainties for individual objects. The native-estimator covariance
sensitivity (F277W/F444W r=0.1671) gives 12.360266, a change of −0.000578;
the independent-band calculation is retained to preserve the existing framework.

Local-density incompleteness, absolute PSF response, correlated noise between
bands, native-estimator residuals, the chosen population boundary and Galactic/
atmosphere systematics remain outside a combined uncertainty. Completed F814W
injection validation supports its aperture probability; it does not validate
full catalog generation or eliminate these limitations.

## Section 4.6 publication paragraph

> We estimated the Galactic ultracool-dwarf foreground using a depth-limited
> effective-volume calculation, following the magnitude- and sensitivity-based
> approaches of [Ryan & Reid (2016)](https://arxiv.org/pdf/1510.05019),
> [Aganze et al. (2022)](https://escholarship.org/content/qt76c609m6/qt76c609m6.pdf)
> and [Chen et al. (2025)](https://arxiv.org/html/2503.18510v1), while retaining
> our adopted Kirkpatrick et al. (2024) temperature-bin densities, Sonora Elf Owl photometry and Galactic
> components. At each position, the intrinsic population was limited to the
> distance at which its F444W aperture flux equals five times the empirically
> calibrated background noise. Within this domain we applied the signed native
> SNR cuts, using the real-image-calibrated F814W probability validated by the
> completed injection–recovery experiment. Over the common 0.532431 deg²
> footprint, the predicted count is 12.36, corresponding to 0.00645 arcmin⁻²;
> the matched four-band prediction is 14.86. The formal propagated uncertainty
> from measured local densities is 0.84 objects, excluding unquantified
> uncertainties in lower-limit density bins and model systematics; the
> conditional finite-F814W-CDF sampling error is approximately 0.03 objects. Changing
> the intrinsic sensitivity boundary from four to six empirical sigma changes
> the prediction from 16.09 to 9.90; this is definition sensitivity rather than
> a confidence interval. The result describes photometric passes within the
> stated depth-limited population, not exact catalog-extraction completeness
> or the final visually vetted sample.

## Execution and reproduction

The central run is reused from `/home/bahareh/results/ali_final_counts` (saved
in the new compact package as `reused_central_execution.json` and its per-bin
table). No old depth caches, atmosphere calculations or injections were rerun.
New commands, from `/home/bahareh/projects/BD_number_density`:

```sh
PYTHON=/home/bahareh/miniforge3/envs/phd/bin/python
$PYTHON -m src.ali_depth_limited --out /home/bahareh/results/ali_depth_limited_sensitivity
$PYTHON -m src.ali_depth_limited --out /home/bahareh/results/ali_depth_limited_convergence --cells 1024 --scenarios convergence
$PYTHON -m src.ali_depth_limited --out /home/bahareh/results/ali_depth_limited_native_correlation --scenarios native_correlation
$PYTHON -m src.ali_selection_uncertainty --maps /home/bahareh/results/ali_final_counts/spatial_maps.npz --calibration /home/bahareh/results/ali_calibration_final --out /home/bahareh/results/ali_depth_limited_cdf_uncertainty.json --nodes 256
$PYTHON -m src.ali_depth_limited_report
$PYTHON -m pytest -q tests/test_ali_depth_boundary.py tests/test_ali_probabilities.py tests/test_ali_domain_validation.py
```

The sensitivity driver reads only existing `ali_final_counts/spatial_maps.npz`,
`ali_calibration_final/*.npz/*.json`, the committed five-band d25 table and fixed
density functions. Local raw data paths and prior calibration/injection commands
remain in `docs/ALI_F814W_REPRODUCE_20261128.md`; those steps are not needed here.
The initial execution additionally retained a legacy r015 correlation stress
test (`correlation`, r=0.224419); the current default uses the measured native
r010 value (`native_correlation`, r=0.1671). They are clearly labelled separately.
