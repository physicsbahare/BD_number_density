# Multi-epoch footprint audit for Dave125 (separate from main number-density model)

**Question:** What fraction of the COSMOS-Web F444W footprint has independent multi-epoch JWST/NIRCam imaging, and how many of Dave's 125 dropout candidates lie in it?

This is a **geometrical opportunity** calculation, not a proper-motion detection probability. Neither two PM detections nor `2/123` / `2/125` measure sky coverage. The full 125 coordinates are frozen in `dave125_positions.csv`, copied from [proper_motion_86](https://github.com/physicsbahare/proper_motion_86/blob/main/data/candidates_step3_magnitude_cut_strict_125.csv) (including 387132 and 1651050).

## Approach
- Query MAST CAOM metadata restricted to `dataRights=PUBLIC` for JWST/NIRCam images within 0.9 degrees of COSMOS (RA 150.12°, Dec +2.20°); **no FITS downloads and no local workstation**.
- Build a primary F444W geometric footprint from COSMOS-Web JWST program 1727.
- Construct three conservative-to-permissive masks for **non-overlapping observation time intervals** (`later.t_min - earlier.t_max >= baseline`):
  1. **same_f444w:** repeated F444W at different epochs
  2. **f444w_plus_red:** F444W at one epoch plus red NIRCam F277W/F356W/F410M/F444W at another epoch
  3. **any_nircam:** any NIRCam filter at each of two epochs
- Intersect each with the primary survey footprint; compute areas (deg²), fractions relative to **archival F444W geometric area**, and numbers of 125 candidate positions within each mask. Save full MAST input metadata, MOC FITS masks, CSV candidate-membership table, JSON numeric summary.
- **Do not equate** any of these masks with an astrometric completeness map. Filters, detection S/N, matching epochs and registration determine whether PM can actually be measured. Cross-filter masks are especially permissive for F444W-only red dropout sources.
- The local estimate of usable F444W catalogue area (about 0.49672 deg²) is **external context only**, not silently substituted for the geometric denominator. Survey usable-area masks/variable depths must be incorporated to compare model and candidate surface density at paper quality.

## Automated implementation
`python -m pm_footprint.compute` with packages `astroquery, astropy, numpy, pandas, mocpy`. GitHub Actions workflow `cosmosweb_pm_footprint.yml` runs on commit, manually, and hourly until a successful result is saved. It runs HEALPix orders 14 and 15, plus a separate exact Shapely/Lambert equal-area geometry calculation for 30, 180, 365 and 730-day baselines. Source, MAST snapshot and candidate hashes are recorded to prevent stale-cache comparisons. `candidate_neighbour_qc.py` flags nearby catalog entries without automatically merging them.

Results directory: `pm_footprint/results/`. **Do not quote an area from this directory as publication-ready until QC is explicitly marked validated.**

References: [MAST astroquery observations](https://astroquery.readthedocs.io/en/stable/mast/mast_obsquery.html), [MOCPy MOC](https://cds-astro.github.io/mocpy/stubs/mocpy.MOC.html).

## Outstanding publication checks
- Verify CAOM observation-level polygons/time are detector/visit level, not deep mosaics that combine multiple epochs; reject malformed/aggregate observations.
- Quantify boundary/HEALPix resolution sensitivity at max_depth 14 vs 15.
- Compare primary coverage to the actual COSMOS-Web F444W catalogue mask and Dave selection/usable depth, including contamination and source deblending.
- Compute selection-matched expected counts in the same effective footprint; keep PM confirmation separately.

## Provisional public-data snapshot (2026-10-29)

The reproducible **PUBLIC-only** MAST query returned 3,168 JWST CAOM observations (2,976 usable imaging rows after excluding grism and invalid times). The exact spherical equal-area Shapely method yields a COSMOS-Web program-1727 archival F444W geometric footprint of **0.533899 deg²**, close to the published nominal 0.54 deg².

| Minimum disjoint epoch gap | F444W+F444W area | F444W+red (F277W/F356W/F410M/F444W) | Any NIRCam | Dave125 candidates inside (same/red/any) |
|---|---:|---:|---:|---:|
| 30 days | 0.074031 deg² (13.87%) | 0.351964 deg² (65.92%) | 0.361733 deg² (67.75%) | 1 / 72 / 75 |
| 180 days | 0.064184 deg² (12.02%) | 0.345759 deg² (64.76%) | 0.356028 deg² (66.68%) | 1 / 72 / 75 |
| 365 days | 0.038698 deg² (7.25%) | 0.317487 deg² (59.47%) | 0.331940 deg² (62.17%) | 0 / 62 / 64 |
| 730 days | 0.012079 deg² (2.26%) | 0.214880 deg² (40.25%) | 0.234175 deg² (43.86%) | 0 / 42 / 44 |

Fractions above use **the 0.533899 deg² archival geometric denominator**, NOT the selection-matched catalogue usable F444W area. For faint brown dwarfs, 365 days is a more physically meaningful first-cut astrometric baseline than 30 days, but even this does **not** guarantee measurable motion. Both confirmed moving-source IDs (387132 and 1651050) lie inside the 365-day red and any-NIRCam masks, not the repeated-F444W mask.

**Do not quote these as the fraction of the usable F444W selection footprint or as a PM confirmation probability.** The independent MOC/HEALPix approximation is boundary-sensitive, especially for narrow same-filter overlap strips (order-14 and order-15 candidate membership differ). The exact Shapely polygon union should be preferred for geometrical areas and source membership. The F444W selection-matched usable-depth mask, actual multi-epoch detection S/N and astrometric registration remain outstanding.

An external estimate of the usable F444W area is about 0.49672 deg², but its **spatial mask and provenance must be verified** before reporting a single usable-area fraction. If that area is a strict subset of the 0.533899 deg² archive footprint, set theory gives conservative *conditional* bounds for the 365-day F444W+red usable-area fraction: **56.43% to 63.92%**, without knowing which pixels are masked. These are bounds, not a measurement.

The 125 candidate entries also include close positional pairs: 115 groups at a 0.5 arcsec linking threshold (62 red-overlap entries at 365 days form 60 positional groups). Proximity alone does not prove duplicate astrophysical sources. Inspect deblending before using candidate counts in density estimates.

Numerical products: `results/exact/summary.json`, `results/exact/dave125_exact_membership.csv`, `results/exact/neighbour_qc.json`, and `results/resolution_qc.json`. See [issue #1](https://github.com/physicsbahare/BD_number_density/issues/1) for ongoing scientific QC and caveats.

## Official DR1 F444W WHT spatial-mask reconstruction (2026-10-29 follow-up)

The [official COSMOS-Web DR1 NIRCam download page](https://cosmos2025.iap.fr/nircam.html) confirms **20 individual F444W 60mas WHT tiles** and a combined **7.0 GB compressed** WHT tarball. The individual links follow the public pattern `https://cosmos2025.iap.fr/data/nircam/extensions/mosaic_nircam_f444w_COSMOS-Web_60mas_{tile}_v1.0_wht.fits.gz`; direct A1 and B10 files were independently checked as present. This is a route to measuring the spatial denominator without the user's workstation.

A new isolated runner `pm_footprint/wht_depth_mask.py` downloads each 60mas tile sequentially, decompresses only that tile, uses its WCS and positive WHT coverage to build a fractionally occupied **0.6-arcsec Lambert equal-area grid**, and intersects the grid with the already archived exact multi-epoch polygons at 30, 180, 365, 730 days. The GitHub Actions workflow `.github/workflows/cosmosweb_wht_area.yml` runs once on creation/push or by manual dispatch; it does **not** requery MAST or repeatedly download 7 GB hourly. It saves only the small numerical summary, source checksums, per-tile provenance and candidate membership. The runner includes WHT>0 and relative WHT>0.1 / >0.5 of each tile's positive median as **threshold sensitivity tests**. A synthetic FITS/WCS test is provided.

**This is an experimental implementation until its 20-tile run succeeds and spatial resolution / WHT threshold checks pass.** WHT>0 is *not* equivalent to a selection-complete F444W 26.2 AB mask; relative per-tile WHT thresholds cannot be relabeled as a magnitude limit without noise/PSF/aperture calibration. The previously quoted 0.49672 deg² may derive from a different PSF-matched WHT or effective-depth definition. Differences between those products must be reported rather than forcing the WHT threshold to reproduce 0.49672 deg². A 0.6-arcsec fractional-grid union of overlapping tiles is approximate; convergence and tile-edge tests are required.

No WHT-derived area fraction or candidate count should be quoted until the output file `pm_footprint/results/wht60/summary.json` exists and its provenance and QC have been reviewed.

### WHT FITS astrometry and rasterization audit (2026-10-29)

All three first WHT Actions attempts failed at A1 because the **standalone extracted WHT FITS primary HDU has no celestial WCS keywords** (verified actual HDU metadata: `EXTNAME=WHT`, `NAXIS1=9600`, `NAXIS2=12455`, and empty WCS key inventory). This is a genuine input-format property, not a missing Python package. The runner now searches the image/primary/combined headers and, when no WCS is present, **reconstructs a TAN transform from the four sky corners for each tile published on the [official DR1 page](https://cosmos2025.iap.fr/nircam.html)**. All 20 corners are frozen in `pm_footprint/dr1_tile_corners.csv`. The fit uses the published common tangent point, verifies ~0.060-arcsec pixel scales, correct parity, and <1-arcsec corner residuals; A1's published corners give ~0.060003 arcsec scales and ~0.00004 arcsec tangent-plane fit residuals for its 9600 x 12455 image. This is a **derived WCS, not a claim to have recovered the original SCI WCS**. Independent comparison to the full i2d SCI header or source matches is still needed before treating sub-arcsecond boundary membership as final.

The original coarse-grid assignment used nearest-cell maxima, which can undercount a rotated tile by leaving artificial holes. It was replaced with per-tile **bilinear fractional-area deposition**, clipping to unit occupancy, and a required <3% native-vs-gridded tile area conservation check for each threshold. This still approximates the union of overlapping tiles; cross-tile overlap and grid offset/resolution convergence remain to be checked. The workflow cancels superseded attempts to avoid redundant multi-GB downloads. None of these changes alters the saved PUBLIC-only exact MAST geometry.
