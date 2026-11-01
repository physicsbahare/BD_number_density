# COSMOS2025 JWST starmask: independent WCS revalidation (2026-11-01)

## What is validated, and what is not

* The original 20-tile, 30mas JWST mask archive was previously processed from the user's attachment. SHA256: `a06107b1c654387bf3ab5284792201cb299b7942f52dc95e41089b0e2967cd52`, 19,370,234 bytes. **The binary tarball is not in GitHub. Never request a replacement or trigger a multi-GB WHT download.**
* Original native **code-1 counts** and **0/125 masked candidate centers** came from `jwst_starmask_audit.py`, which uses the correct row-scaled FITS PC/CDELT convention. They are *not* affected by the rasterizer bug.
* Independent COSMOS-Web WHT/epoch geometric result remains **0.314624457 / 0.528823685 = 59.495%** with 62/125 candidate positions having >=365d F444W+red (F356W) coverage. The local v13 analysis independently finds ~0.317401/0.533619 = 59.48% from 3.6-arcsec two-phase sampling. These values mean geometric imaging opportunities, NOT astrometric detection efficiency or selection completeness.
* The **old** 0.3/0.6-arcsec mask union, old weighted mask points, and derived tight mask-corrected fraction were computed before the PC/CDELT correction and are **not revalidated**. No exact spatial WHT × corrected-star-mask intersection has yet been measured. Do not use the earlier 0.0128-0.0131 deg² mask union as an exact published area.

## New independent safety check

`mask_astropy_parity.validate_fits_header_wcs` checks a full real FITS header at **nine actual native-pixel locations** (corners, edges, center), independently comparing the analytic TAN mapping to `astropy.wcs.WCS.all_pix2world(..., origin=0)`. It also compares the Astropy pixel transformation matrix, calculates the true TAN local pixel solid angle, and **refuses non-TAN or SIP-distorted images**. Tolerance: **0.001 arcsec** per location.

The mask rasterizer now executes this check on **each of the 20 real mask headers before depositing any flag pixels**, writing `astropy_wcs_parity` diagnostics into each tile entry. The provenance gate rejects an output missing any of the 20 nine-point checks, or showing WCS/solid-angle inconsistency. This guards subsequent epoch and sparse-WHT intersection scripts against old or wrong-WCS inputs.

Synthetic tests: `python -m unittest -v pm_footprint.test_mask_wcs pm_footprint.test_mask_astropy_parity pm_footprint.test_mask_product_integrity`. The lightweight GitHub Actions test workflow includes Astropy; neither test nor workflow downloads the star mask or WHT. Passing synthetic tests **does not substitute for reprocessing the actual archive**.

## Reproduce on the EXISTING archived file, only when already accessible

```bash
cd ~/BD_number_density
git pull
source .venv/bin/activate
# Use the previously supplied exact archive wherever it is already staged:
MASK_TAR="/existing/path/cosmos_web_starmask_jwst.tar.gz"
sha256sum "$MASK_TAR"
python -m pm_footprint.jwst_mask_union_repro \
  --archive "$MASK_TAR" \
  --output pm_footprint/results/jwst_starmask/corrected_wcs_union.json \
  --sample-out pm_footprint/results/jwst_starmask/corrected_wcs_points.csv \
  --samples 20000 --seed 20261101
python -m pm_footprint.jwst_mask_epoch_from_sample \
  --points pm_footprint/results/jwst_starmask/corrected_wcs_points.csv \
  --mask-summary pm_footprint/results/jwst_starmask/corrected_wcs_union.json \
  --output pm_footprint/results/jwst_starmask/corrected_epoch_intersection.json
```

**Acceptance checks:** archive SHA equals frozen SHA; exactly 20 tiles; each tile has a measured `astropy_wcs_parity` with nine positions and max WCS deviation <=0.001 arcsec; `wcs_pc_off_diagonal`, `legacy_column_scaled_pc_shift_arcsec_max_corners`, and `wcs_center_area_ratio_to_nominal_30mas` are recorded and reviewed; 0.3/0.6 grid area conservation and sky-union resolution/shift sensitivity reviewed; repeated 125-source center check remains consistent; code-1 exclusion semantics confirmed from originating reduction documentation. Preserve the output summary and weighted-point hash as a paired provenance unit.

For a **direct WHT-masked geometric overlap** (not selection completeness), additionally run `pm_footprint.jwst_mask_wht_sparse_intersection` against the already-existing local Zhen-Kai v13 WHT, which samples only four bytes per point. This is not the official 60mas WHT result; label the reduction used. Official DR1 full-WHT occupancy cells were not saved, and cannot be reconstructed exactly from only the scalar area summary.

## Manuscript relevance

The appropriate, currently supportable statement is approximately **59.5% geometric multi-epoch coverage** with **62/125** candidates having >=1-year red-band imaging. Do **not** claim 59.5% proper-motion completeness; do not apply unvalidated star-mask numerical corrections; do not revise the published model-predicted approximately two *photometrically selected* ultracool dwarfs using a geometric fraction.

PM source 387132 has robust repeated-exposure motion evidence. The 1651050 result is a significant **apparent** cross-epoch displacement, but intrinsic PM versus parallax/centroid systematics cannot be resolved from the two-epoch positional measurement alone. Significant, independently registered shifts can strongly favor Galactic rather than stationary distant-galaxy origins; the strength of a claim of specifically intrinsic PM and a **confirmed brown dwarf** must follow the actual astrometric and SED evidence, not wording preference.
