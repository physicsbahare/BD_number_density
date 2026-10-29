# Multi-epoch footprint audit for Dave125 (separate from main number-density model)

**Question:** What fraction of the COSMOS-Web F444W footprint has independent multi-epoch JWST/NIRCam imaging, and how many of Dave's 125 dropout candidates lie in it?

This is a **geometrical opportunity** calculation, not a proper-motion detection probability. Neither two PM detections nor `2/123` / `2/125` measure sky coverage. The full 125 coordinates are frozen in `dave125_positions.csv`, copied from [proper_motion_86](https://github.com/physicsbahare/proper_motion_86/blob/main/data/candidates_step3_magnitude_cut_strict_125.csv) (including 387132 and 1651050).

## Approach
- Query the public MAST CAOM observational metadata for JWST/NIRCam images within 0.9 degrees of COSMOS (RA 150.12°, Dec +2.20°); **no FITS downloads and no local workstation**.
- Build a primary F444W geometric footprint from COSMOS-Web JWST program 1727.
- Construct three conservative-to-permissive masks for genuinely separate timestamps (>=30 days):
  1. **same_f444w:** repeated F444W at different epochs
  2. **f444w_plus_red:** F444W at one epoch plus red NIRCam F277W/F356W/F410M/F444W at another epoch
  3. **any_nircam:** any NIRCam filter at each of two epochs
- Intersect each with the primary survey footprint; compute areas (deg²), fractions relative to **archival F444W geometric area**, and numbers of 125 candidate positions within each mask. Save full MAST input metadata, MOC FITS masks, CSV candidate-membership table, JSON numeric summary.
- **Do not equate** any of these masks with an astrometric completeness map. Filters, detection S/N, matching epochs and registration determine whether PM can actually be measured. Cross-filter masks are especially permissive for F444W-only red dropout sources.
- The local estimate of usable F444W catalogue area (about 0.49672 deg²) is **external context only**, not silently substituted for the geometric denominator. Survey usable-area masks/variable depths must be incorporated to compare model and candidate surface density at paper quality.

## Automated implementation
`python -m pm_footprint.compute` with packages `astroquery, astropy, numpy, pandas, mocpy`. GitHub Actions workflow `cosmosweb_pm_footprint.yml` runs on commit, manually, and hourly until a successful result is saved. The first output is provisional until the MOC resolution-sensitivity test and footprint/visit QC are done.

Results directory: `pm_footprint/results/`. **Do not quote an area from this directory as publication-ready until QC is explicitly marked validated.**

References: [MAST astroquery observations](https://astroquery.readthedocs.io/en/stable/mast/mast_obsquery.html), [MOCPy MOC](https://cds-astro.github.io/mocpy/stubs/mocpy.MOC.html).

## Outstanding publication checks
- Verify CAOM observation-level polygons/time are detector/visit level, not deep mosaics that combine multiple epochs; reject malformed/aggregate observations.
- Quantify boundary/HEALPix resolution sensitivity at max_depth 14 vs 15.
- Compare primary coverage to the actual COSMOS-Web F444W catalogue mask and Dave selection/usable depth, including contamination and source deblending.
- Compute selection-matched expected counts in the same effective footprint; keep PM confirmation separately.
