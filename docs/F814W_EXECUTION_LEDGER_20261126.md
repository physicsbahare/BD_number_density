# F814W reconstruction execution ledger

**2026-11-28 update:** all twenty raw-image reconstructions and targeted tests
are complete. The preserved-convention conditional count is 23.84371388 over
0.5324308307 deg². All specified numerical runs and 23 selected regression
tests finished. Historical-baseline/exact-native-selection certification is
unresolved; see `data/ali_f814w_2026-11-26/REPORT.md`. Earlier references to an
accepted reconstruction below refer to processing stages, not scientific
acceptance of a final paper prediction. The primary continuation preserves
r015 NIRCam empirical-noise selection and uses native r010 F814W calibration;
all-native alternatives are diagnostics only.

The user authorized regeneration of the missing NIRCam caches on 2026-11-26.
The earlier missing-cache blocker is superseded. Survey images remain local
on E_202; no survey image download is attempted.

## Recovered and verified

- Original FC-ENZO solar-metallicity JWST/HST tables were joined on model
  name (7200 models). Their median sqrt(flux_zp25)*Radius_pc agrees with all
  26 archived F115W/F444W bin-reach values to their recorded precision.
- Original F814W response support is 6737.5–9778.5 Angstrom, entirely inside
  the Elf Owl 0.6–15 micron grid. Its pivot is 8045.528 Angstrom versus
  8047.256 in the current ACS header. This preserves the original synthetic
  photometry; no alternative atmosphere grid is introduced.
- NIRCam science is MJy/sr and uses PIXAR_SR to obtain nJy per pixel.
  HST science/ERR are nJy and WHT has inverse-nJy-squared units. WHT is used
  as an empirical depth coordinate; reciprocal total ERR variance is not
  assumed. Full images are streamed; only sampled scenes, grids and packed
  native validity masks are cached.
- Native 0.10 arcsec-radius fractional apertures are needed for the stored
  SNR columns. The original repository also specifies a 0.15 arcsec depth
  radius; it is reconstructed separately for historical comparison.
- Paired-source tests show that noise-equalized aperture sums reproduce
  catalog SNR. AUTO flux/error does not. Subpixel forced apertures reduce
  first-tile NIRCam held-out SNR median absolute errors to 0.03–0.06.
- ACS FK5 metadata creates a ~0.029 arcsec offset relative to NIRCam ICRS.
  A static paired-source registration independently finds ~(+0.89,+0.31)
  native pixels on A1, compensating this frame-label difference. This is a
  static alignment check, not a proper-motion or parallax measurement.

## Corrections made during reconstruction

- The streamed reader initially overflowed when applying the unsigned
  segmentation map's BZERO. A failing unsigned-FITS regression was added;
  conversion before scaling fixes it.
- Initial PSF controls incorrectly used flag_star as a stellar class.
  That flag is a bright-star contamination mask. Those controls were
  invalidated. Accepted controls use LePHARE type=1, unmasked sources,
  optical SNR and compactness cuts. Core encircled-energy comparisons use
  identical background treatment on models and observed scenes.
- Initial exploratory counts used a generic extended distance bound and
  a simple aperture-SNR scale ratio. They are not final results and are
  retained only as local diagnostic previews. The accepted reconstruction
  uses the existing depth-dependent maximum-distance model and a directly
  verified native noise-equalized SNR estimator.

## Numerical-domain ruling

The repository's count_to_distance integrates from zero to the detection-
depth-dependent dmax. The reconstruction retains that detectable-population
domain per spatial cell, using the empirical five-sigma F444W reach and
the same atmospheric normalization. Gaussian photometric probabilities are
applied inside that domain. A regression checks agreement with the existing
distance integral in a limiting case. Extended-distance tests are diagnostic:
an untruncated halo and a nonzero Gaussian upward-fluctuation tail otherwise
introduce dependence on an arbitrary outer distance. The depth restriction
must be explicit in the final method and interpretation; no reproduction
of an undocumented soft-bound implementation is claimed.

## Paths and stages

- Raw: /mnt/E_202/JWST/COSMOS-Web/NIRCam_mosaics and
  /mnt/E_202/JWST/COSMOS_Web/HST_F814W_v0.5.2026jul/tiles.
- Streamed inputs: /home/bahareh/results/ali_fiveband_reconstruction.
- Native valid-aperture footprint: /home/bahareh/results/ali_native_coverage.
- Original positive-weight-center footprint: /home/bahareh/results/ali_original_coverage.
- Final calibration: /home/bahareh/results/ali_calibration_final.
- Historical-radius calibration: /home/bahareh/results/ali_calibration_legacy_r015.
- Accepted-convention targeted tests: /home/bahareh/results/ali_f814w_final_recovery.

Directories named preview, catalogq, r010, or r015 without the final/legacy
qualifiers above contain intermediate diagnostics, not accepted final counts.
No all-field count or completion claim is made until all twenty fields,
matched native footprints, convergence and uncertainty checks finish.
