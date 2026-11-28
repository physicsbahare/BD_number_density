# Ali native-domain and detection-footprint validation review, 2026-11-28

Scope: read-only review of `src/ali_domain_validation.py`, `/home/bahareh/results/ali_domain_validation/*.json`, `native_radial_shell_counts.csv`, and `/home/bahareh/results/ali_detection_validation/detection_footprint_validation.json`. No code edits, new depth reconstruction, injections, atmosphere changes or Galactic changes.

## Finding

No scientifically identified finite **all-distance source-associated photometric count** follows from the currently specified untruncated halo and forced-photometry selection kernel. A finite depth-limited count remains mathematically valid as an explicitly conditional estimand, but no existing ingredient reviewed here uniquely selects its outer boundary. Do not substitute an arbitrary cutoff, the 50% detection radius, or a floor-subtracted excess for the requested prediction.

This is stronger than numerical cutoff sensitivity: with the specified model the unsubtracted all-distance integral diverges. It does **not** mean an actual finite-pixel catalog contains infinitely many objects. The divergent integral assigns a nonzero probability of a chance noise selection separately to infinitely many essentially invisible model stars, rather than modeling a finite catalog of detected peaks and their source associations.

## Mathematical check

At large heliocentric distance along the fixed sightline, the halo density satisfies rho(r) proportional to r^-2.77. Modeled source flux and source photon variance vanish as r^-2. The probability function in `ali_domain_validation.py:67-71` therefore approaches the nonzero source-free selection probability p0. Consequently the halo count integrand approaches a positive constant times r^(2-2.77) = r^-0.77, and cumulative count grows as R^0.23. The thin and thick disks have convergent tails, but their convergence does not rescue the nonzero halo term. This conclusion does not depend on the particular diagnostic cutoff.

The empirical-CDF interpolation and Gaussian source-noise convolution are continuous in the zero-signal limit used here. Area-weighted modeled p0 = 0.0071740392. Real held-out scenes give 526/84,203 = 0.0062468083 joint passes. Those scenes establish a nonzero chance-pass phenomenon; they are reconstructed native-estimator forced-aperture measurements, not independent source-free entries extracted by the actual catalog pipeline. Their pooled sampling distribution differs from area weighting, so their discrepancy with modeled p0 is not a direct calibrated goodness-of-fit test. The ACS estimator's previously documented residual also still applies.

For a separate source-associated detection kernel D to regularize the same halo integral, it must vanish fast enough at large r. If D is asymptotically r^-alpha, convergence requires alpha > 0.23. If D is proportional to a positive power beta of source flux, beta > 0.115 suffices because flux scales as r^-2. Merely multiplying by an ordinary recovery curve with a nonzero chance-association floor would still diverge. This is a necessary asymptotic condition, not a proposal to invent or fit a new completeness curve.

## Numerical checks

Independently summed the saved absolute-shell rows. Cumulative five-band empirical predictions to 20/50/100/200/1000 kpc are 61.322510 / 82.962169 / 102.097881 / 124.837992 / 195.330299. At 1 Mpc the integrated zero-source floor is 143.774745 and the floor-subtracted excess is 51.555554. These match the supplied interpretation.

The excess appears numerically convergent because the leading positive constant tail is removed. It is a different, signed response/excess estimand. In general, individual-color contributions can have either sign, so it is not automatically a nonnegative source-completeness probability or an accepted selected-source count. No adoption is justified by its convergence alone.

The native 50% F444W probability boundary is correctly formed at source aperture flux = 5*catalog_denominator - background (`ali_domain_validation.py:88`), conditional on the Gaussian estimator model. Its boundary is not a physical support boundary. Scaled-shell totals are 28.706317 inside 0.001-1 times reach; 15.071112 in 1-1.5; 7.071083 in 1.5-2; 16.262648 in 2-5; and 13.470587 in 5-10. Material selection weight beyond reach rules out claiming that the original hard bound exhausts the soft selection. The positivity clamp on the boundary denominator is numerically defensive; the covered-map audit found no nonpositive denominator, so it does not create the diagnostic boundary.

`ali_domain_validation.py:81-86` uses the correct log-radius Jacobian r^3, the existing density components, and the grouped area-weighted selection probability. Source-free floor and source-inclusive counts are integrated with the same density and geometry, making their subtraction a coherent diagnostic. The positive exponent 0.23 is an asymptotic analytical result; numerical integration to 1 Mpc need not itself reach the exact asymptotic normalization.

## Can existing detection documentation supply a finite boundary?

The catalog paper describes a peak statistic on a combined detection image followed by hot/cold extraction. Its completeness assessment uses F277W-selected comparison objects. These do not establish a hard limit on intrinsic F444W flux or a source-associated detection kernel for the F444W-only PSF population. [COSMOS2025, sections 3.3-3.4](https://arxiv.org/html/2506.03243v1#S3.SS3).

A measured-image detection threshold is compatible with upward noise fluctuations from arbitrarily faint sources; it is not a true-flux cutoff. A finite catalog-generation/association model could change the counting problem, but that model is not present in the current population integral. The existing documented F277W galaxy completeness must not be silently reused as F444W-only brown-dwarf completeness.

The original `count_to_distance` convention supplies a finite user-specified domain, not a derivation of the domain from these native selection probabilities. Keeping that convention is acceptable for a clearly labeled depth-limited reconstruction; it cannot resolve the all-distance exact-native-selection request or validate the missing historical driver.

## Detection-footprint check and historical attribution

The reviewed detection-footprint output distinguishes:

- Original positive-WHT four-band centers: 0.533333202965 deg2.
- Native aperture-valid four bands: 0.532647460575 deg2.
- Native aperture-valid five bands: 0.532430830667 deg2.
- Requiring detection-image positivity on original four-band centers: 0.499244068620 deg2.
- Requiring positivity on native five-band coverage: 0.498388803518 deg2.

Detection-image positivity is an intensity/noise selection, not a coverage criterion; covered pixels can have nonpositive detection-image values. This audit identifies a plausible class of erroneous area reductions, but does not identify the historical cause. In particular, 0.499244068620 is not the archived 0.496306 deg2. No driver, exact historical mask, or equivalent source-selection configuration has thereby been recovered. Do not call the historical discrepancy explained, reproduced or fixed.

The recovered archive pass reports verification of all FITS bitmaps and grids after a trailing non-FITS member; unpersisted per-tile header diagnostics correctly remain null. That supports the stated artifact recovery, not fabricated complete per-tile header evidence.

## Review conclusion

The root's proposed conclusion is justified with precise wording: the all-distance forced-photometry population integral is divergent under the fixed untruncated halo and nonzero blank-selection floor; the requested finite source-associated/catalog prediction remains unidentified without a justified domain or suitable source-associated selection model. No valid finite all-distance domain was found in the existing model or documented catalog method. Preserve earlier counts as conditional diagnostics and retain final acceptance as unresolved. Historical mask/driver attribution remains unresolved independently of this divergence result.
