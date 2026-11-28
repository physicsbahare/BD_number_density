# Reproduce the executed F814W reconstruction

The execution environment and exact package versions are in
`data/ali_f814w_2026-11-26/execution_environment.json`. The server interpreter
is `/home/bahareh/miniforge3/envs/phd/bin/python`. No survey downloads are used.
The local directory was recovered as a source snapshot, not a complete git
checkout; the remote commit preserves all existing repository files.

The result is conditional. Read `data/ali_f814w_2026-11-26/REPORT.md` before
using any number in a paper. A completed numerical run does not signify a
passed historical-baseline validation.

## Inputs and retained products

- NIRCam SCI/ERR/WHT: `/mnt/E_202/JWST/COSMOS-Web/NIRCam_mosaics/`, twenty
  A1–A10/B1–B10 tiles, four filters, 30 mas, v1.0.
- HST SCI/ERR/WHT: `/mnt/E_202/JWST/COSMOS_Web/HST_F814W_v0.5.2026jul/tiles/`,
  `mosaic_cosmos_fulldepth_hst_acs_wfc_f814w_30mas_tile_{tile}_v0.5.2026jul_{drz,err,wht}.fits.gz`.
  The similarly named older failed-download directory contains an empty file
  and is not used.
- PSFs: `/mnt/E_202/JWST/COSMOS-Web/auxiliary/psfs_v5.0/`, per-tile native
  PSFEx HST F814W and JWST filter products. Full paths are in calibration JSON.
- Segmentation: `/mnt/E_202/JWST/COSMOS-Web/auxiliary/segmentation_maps/`,
  `detection_chi2pos_SWLW_{tile}_segmap_v1.3.fits.gz`.
- Catalog: `/home/bahareh/Data/Cosmos_Web/COSMOSWeb_mastercatalog_v1.1.fits`.
- JWST star masks: `/mnt/E_202/JWST/COSMOS-Web/auxiliary/cosmos_web_starmask_jwst.tar.gz`;
  optical star regions: adjacent `cosmos2020_starmask_hsc.reg`. Star-mask
  diagnostics are not silently added to the five photometric cuts.
- Original atmosphere photometry provenance: `data/elfowl_fiveband_provenance.json`;
  pinned upstream FC-ENZO solar HST/JWST tables, joined on model name. Server
  joined cache: `data/elfowl_solar_fiveband_original.tsv` (7200 rows).
  Adopted compact thirteen-bin reaches: `data/elfowl_fiveband_bin_reach.csv`.
  Original response: `data/original_F814W_response.dat`.
- Representative wavelength-support check: `/mnt/E_202/Brown_Dwarf_Models/Sonora_ElfOwl_v2_Y/teff_275_325/spectra_logzz_2.0_teff_275.0_grav_17.0_mh_1.0_co_1.0.nc`;
  h5py reads wavelength 0.6000015–14.9997608 micron. Response support is
  0.67375–0.97785 micron. This file is a wavelength-coordinate check, not a
  replacement population grid.

The 300-image inventory, sizes and timestamps are saved in
`survey_input_manifest.json`; no complete raw-FITS checksum is claimed.
Full WCS translations, native-mask definitions and calibration diagnostics
are also committed. Full image/stamp/packed-mask caches remain on the server
under `/home/bahareh/results/`, not in GitHub.

## Execution

From the repository root:

```sh
export OPENBLAS_NUM_THREADS=1
PY=/home/bahareh/miniforge3/envs/phd/bin/python
$PY -m src.depth_reconstruction
$PY -m src.ali_depth_analysis --out /home/bahareh/results/ali_calibration_final --radius .1
$PY -m src.ali_depth_analysis --out /home/bahareh/results/ali_calibration_legacy_r015 --radius .15
$PY -m src.ali_native_coverage --registration-dir /home/bahareh/results/ali_calibration_final
$PY -m src.ali_original_footprint
$PY -m src.f814w_targeted_recovery --calibration /home/bahareh/results/ali_calibration_final --out /home/bahareh/results/ali_f814w_final_recovery
$PY -m src.f814w_population_validation --calibration /home/bahareh/results/ali_calibration_final --out /home/bahareh/results/ali_f814w_population_validation
$PY -m src.ali_f814w_validation_report --recovery /home/bahareh/results/ali_f814w_final_recovery --population /home/bahareh/results/ali_f814w_population_validation --out data/ali_f814w_2026-11-26
$PY -m src.ali_empirical_validation --calibration /home/bahareh/results/ali_calibration_final --recovery /home/bahareh/results/ali_f814w_final_recovery --out data/ali_f814w_2026-11-26
$PY -m src.ali_count_run --calibration /home/bahareh/results/ali_calibration_final --nircam-calibration /home/bahareh/results/ali_calibration_legacy_r015 --snr-mode preserved_nircam --coverage /home/bahareh/results/ali_native_coverage/area_maps.npz --original-coverage /home/bahareh/results/ali_original_coverage/original_four_area_map.npz --out /home/bahareh/results/ali_preserved_nircam_counts --cells-per-tile 512
```

Image sampling uses 10,000 blank candidates per tile, a 100 native-pixel
(3 arcsec) spatial grid and 128-row gzip chunks. Existing COMPLETE markers
and per-tile caches are reused. Do not point a changed aperture radius at an
existing calibration directory; the two radii have separate directories.
The unsigned BZERO segmentation scaling is regression-tested. Blank scenes
are segmentation-clear and more than 0.75 arcsec from catalog objects.
Original even blank indices calibrate noise; the held-out scenes are split
again into even CDF-estimation and odd independent injection-validation
indices. RNG seeds are explicit in the sampling and injection modules.

The primary command retains reconstructed r015 NIRCam empirical-noise SNR
and zero-background probability convention. F814W uses the native r010
catalog-estimator scale, empirical residual CDF and source-noise convolution.
It is not an assertion that empirical-noise NIRCam SNR equals catalog SNR.
The intrinsic radial upper bound is the retained five-sigma empirical-noise
F444W reach. `--radial-domain catalog_depth` and `extended` are sensitivity
diagnostics, not adopted changes to the population model.

For the historical-radius comparison, use the r015 directory as
`--calibration`, omit `--nircam-calibration`, and set `--snr-mode empirical
--zero-background --selections four_original four_matched`.
For the all-native diagnostic, use the final r010 calibration with
`--snr-mode catalog` and no NIRCam override. Run names and exact flags are
recorded in `executed_runs.json`. The extended native diagnostics use
`--radial-domain extended --rmax 20000` and `50000`, with
`--selections four_matched five_gaussian`.

## Convergence, sensitivity and final collection

Repeat the primary count command with separate output directories and these
changes. Completed results are retained; no raw-image regeneration is needed.

| Server output directory | Change from primary command |
|---|---|
| `ali_preserved_convergence1024` | `--cells-per-tile 1024` |
| `ali_preserved_shot31` | `--shot-nodes 31 --selections five_empirical` |
| `ali_preserved_noise_low/high` | `--f814w-noise-scale .95/1.05 --selections five_empirical` |
| `ali_preserved_threshold_low/high` | `--f814w-threshold-shift -.5/.5 --selections five_empirical` |
| `ali_preserved_flux_low/high` | `--f814w-flux-scale .95/1.05 --selections five_empirical` |
| `ali_preserved_correlation` | `--nircam-correlation .224419 --selections four_matched five_empirical` |

Each `.95/1.05` entry represents two distinct runs, not literal shell syntax.
All primary sensitivity runs use 512 spatial leaves per tile. Native r010
correlation diagnostics use 0.1671; the retained r015 correlation diagnostic
uses its separately measured 0.224419.

```sh
$PY -m src.ali_selection_uncertainty --maps /home/bahareh/results/ali_preserved_nircam_counts/spatial_maps.npz --calibration /home/bahareh/results/ali_calibration_final --out /home/bahareh/results/ali_preserved_cdf_uncertainty256.json --nodes 256 --preserved-nircam
$PY -m src.ali_selection_uncertainty --maps /home/bahareh/results/ali_preserved_nircam_counts/spatial_maps.npz --calibration /home/bahareh/results/ali_calibration_final --out /home/bahareh/results/ali_preserved_cdf_uncertainty512.json --nodes 512 --preserved-nircam
$PY -m src.ali_footprint_phase_bound --coverage /home/bahareh/results/ali_native_coverage --out data/ali_f814w_2026-11-26/footprint_phase_bound.json
$PY -m src.ali_mask_audit --archive /mnt/E_202/JWST/COSMOS-Web/auxiliary/cosmos_web_starmask_jwst.tar.gz --coverage /home/bahareh/results/ali_native_coverage --out data/ali_f814w_2026-11-26/star_mask_audit.json
$PY -m src.ali_verify_atmosphere_cache --models data/elfowl_solar_fiveband_original.tsv --archived data/preliminary_cosmosweb_teff_counts.csv --out /home/bahareh/results/ali_regenerated_bin_reach.csv
$PY -m src.ali_f814w_count_figures --counts /home/bahareh/results/ali_preserved_nircam_counts --coverage /home/bahareh/results/ali_native_coverage/area_maps.npz --out data/ali_f814w_2026-11-26
$PY -m src.ali_f814w_extra_figures --calibration /home/bahareh/results/ali_calibration_final --out data/ali_f814w_2026-11-26
$PY -m src.ali_f814w_finalize --out data/ali_f814w_2026-11-26
$PY -m pytest -q tests/test_ali_probabilities.py tests/test_depth_reconstruction.py tests/test_native_coverage.py tests/test_geometry.py tests/test_injection_recovery_core.py tests/test_uncover_validation.py
```

The CDF uncertainty is conditional on calibration and uses a step ECDF; its
central value differs by 0.02439 objects from the main interpolated-ECDF count.
This difference is explicitly reported rather than absorbed into a small
quadrature error. The selected existing and new regression tests pass (23).
The unrelated dropout-noise test could not be collected in the partial local
snapshot because `pm_footprint` is absent; no full-repository test-suite claim
is made and no unrelated project work is undertaken.

The exact historical driver/mask/configuration remains unrecovered. Its
filename is not recorded. This prevents certification of the requested
historical-baseline and exact-five-band acceptance checks, despite successful
execution of the reconstructed calculations.
