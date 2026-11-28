# Adopted depth-limited photometric foreground approximation

See the [scientific report and Section 4.6 paragraph](../../docs/ALI_DEPTH_LIMITED_20261128.md)
and [limiting-distance literature review](../../docs/ALI_DEPTH_LIMITED_LITERATURE_20261128.md).

`result.json` adopts the population with intrinsic F444W aperture signal at
least five empirical background sigma. The result is not an unrestricted
all-distance or exact catalog-selection prediction. No fixed kpc cap is applied.
The earlier divergence remains valid for the unrestricted estimand.

- `reused_central_execution.json` and `central_temperature_component_counts.csv`: completed native 512-cell run, reused without repeating it. Its `rmax_pc` metadata is ignored by the `depth` integration branch.
- `five_band_bins_and_effective_volumes.csv`: adopted five-band contributions and selected effective volume, N divided by the adopted total local density for that bin. This volume includes population-component fractions and selection probabilities.
- `sensitivity_summary.csv` and `sensitivity_temperature_bins.csv`: changed intrinsic depth, ACS noise/threshold and covariance diagnostics. They are not confidence bounds. `correlation` is the legacy r015 stress value 0.224419; `native_correlation` uses 0.1671.
- `*_distance_boundaries.csv`: area-weighted spatial 16/50/84 percentiles of the local limiting distances, not inferred individual-source distances. `r50_pc` here denotes a spatial median and should not be confused with the native-SNR 50%-pass radius scenario `native_r50`.
- `cdf_sampling_uncertainty.json`: sampling error in the finite empirical ACS CDF, conditional on the adopted calibrations and independent-scene approximation; shared CDF uncertainty across temperature bins is retained.
- `provenance.json`: environment and hashes of preserved small inputs. No new injections or raw-image depth reconstruction.
- SVG figures: depth-definition sensitivity and component contributions by temperature.

The formal measured-density error omits unquantified lower-limit-bin errors.
Noise/PSF/calibration, model and domain-definition effects are not absorbed into
a single reported confidence interval. Exact AUTO positivity, catalog detection,
morphology, quality flags and visual-vetting completeness are not modeled.

The calibration and 851,643-injection results remain in the earlier committed
F814W package and were not repeated. Large images and spatial caches remain on
the server; execution commands and input paths are documented in the report.
