# Cached DR1 F444W WHT threshold-sensitivity audit (2026-10-30)

**No new WHT/FITS downloads, no MAST query, no new area measurement.**
All numbers below are direct arithmetic on the already-committed
[`summary.json`](./summary.json) (the 20-tile 60mas run). This is a
*threshold-sensitivity diagnostic* of the **same coarse 0.6-arcsec
fractional grid**, not a 26.2 AB completeness calibration.

| WHT threshold (relative to each tile's positive median) | F444W denominator (deg²) | >=365d F444W+red intersection (deg²) | Red geometric fraction | Dave entries in red mask |
|---|---:|---:|---:|---:|
| >0 | 0.528823685 | 0.314624457 | 59.4952% | 62/125 |
| >0.1 × tile median | 0.528743433 | 0.314587928 | 59.4973% | 62/125 |
| >0.5 × tile median | 0.517528902 | 0.308692384 | 59.6474% | 62/125 |

Relative to WHT>0, the 0.5-median cut reduces the gridded denominator by
**2.1358%**, the red-overlap numerator by **1.8854%**, and changes the
reported red fraction by **+0.1522 percentage points**. All 125 candidate
positions remain inside the F444W mask at these three thresholds; the 62
red-overlap entries and 64 any-NIRCam entries are unchanged.

This stability is **not evidence that the field is complete at F444W=26.2
AB**. Relative tile medians are not calibrated to aperture-noise 5-sigma
depth; the blue-band non-detection cuts are source-color/noise dependent;
the JWST bright-star mask and blends have not been applied; and a common
**~0.955–0.959% native-pixel-to-grid deposition/clipping loss per tile**
remains uncorrected. The 0.105% post-grid-resolution agreement is not a
measure of this systematic. Do not interpret the 0.528824 vs 0.533899
deg² difference as a physical gap or report sub-percent absolute accuracy.

Even the 0.5-median WHT area, **0.517529 deg²**, exceeds the externally
quoted 0.49672 deg² by **0.020809 deg²**. Therefore none of the tested
relative thresholds reproduces the external effective area, and there
is no justification for tuning a threshold to force agreement. The
0.49672 deg² figure needs spatial provenance and selection definition.

The independent cached candidate-by-filter pointcheck reports **all 62**
>=365d red-covered candidate entries have F356W geometric imaging
opportunity (zero are F277W-only), but it provides **no per-filter
area**. Do not substitute 0.314624 deg² (all-red WHT overlap) for the
F356W-only area until `filter_epoch_qc/summary.json` is generated and
cross-validated against the exact union.

**Pending**: official COSMOS2025 JWST spatial starmask, F444W noise-to-depth
calibration at 26.2 AB, F115W/F150W/F277W dropout selection function,
source-level flags/blending and F356W SNR/astrometric registration. The
official DR1 release announces starmasks but its catalog-download page
currently presents a form and captcha, so direct reproducible mask access
has not yet been established.

References:
- [Official COSMOS2025 catalog and flag definitions](https://cosmos2025.iap.fr/catalog.html)
- [Official DR1 release announcements](https://cosmos2025.iap.fr/)
- [Official catalog access form](https://cosmos2025.iap.fr/catalog_download.php)
- [Frozen exact public geometry](../exact/summary.json)
- [Independent candidate filter QC](../filter_epoch_qc/independent_candidate_pointcheck.json)
