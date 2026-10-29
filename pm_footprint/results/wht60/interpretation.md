# DR1 F444W WHT audit: successful 20-tile run, provisional interpretation

**2026-10-29 UTC, from [summary.json](summary.json), [grid_resolution_qc.json](grid_resolution_qc.json), [tile_progress.json](tile_progress.json), [candidate membership](dave125_wht_membership.csv).**

This analysis is a *geometric opportunity* audit. It **does not** estimate proper-motion detection efficiency or selection completeness for the F444W < 26.2 AB dropout sample.

## Reproducible WHT>0 intersections

The public 60 mas F444W DR1 WHT processing completed **all 20 tiles A1–B10**. Each tile's original `i2d SCI` WCS was read from a small header prefix; the standalone WHT extension has no WCS. MAST PUBLIC-only exact polygons and all 125 candidate positions were reused from frozen inputs.

- Exact MAST F444W archive polygon: **0.533898955 deg²**.
- F444W WHT>0 **inside that polygon** (fractional 0.6 arcsec equal-area grid): **0.528823685 deg²** (99.049% of archival area). This is a **numerical raster estimate**, not a demonstrated 0.95% physical coverage gap.
- F444W WHT > 0.5 × each tile's median positive WHT: **0.517528902 deg²**; this threshold is arbitrary and is **not** calibrated to F444W=26.2 AB completeness.

| Minimum independent epoch interval | Same F444W area (deg²), count | F444W+red area (deg²), count | Any NIRCam area (deg²), count |
|---|---|---|---|
| 30 days | 0.073347, 1/125 | 0.348777, 72/125 | 0.358448, 75/125 |
| 180 days | 0.063590, 1/125 | 0.342634, 72/125 | 0.352800, 75/125 |
| **365 days** | **0.038337, 0/125** | **0.314624, 62/125** | **0.328935, 64/125** |
| 730 days | 0.011966, 0/125 | 0.212973, 42/125 | 0.232086, 44/125 |

At **365 days**, the F444W+red opportunity area is **59.495% of WHT>0**, versus **59.466% of the archival polygon** using exact geometry alone. Any NIRCam is **62.201% of WHT>0**, repeated F444W **7.250%**. Both PM-confirmed IDs 387132 and 1651050 are in the red cross-filter mask, neither in repeated F444W. These are positions where imaging epochs exist, **not** proof that faint red objects are detected in both epochs.

The WHT>0, >0.1 median and >0.5 median masks produce the same candidate counts at all baselines; for the 365-day red mask their *area fractions* are 59.495%, 59.497%, and 59.647%, respectively. The 62 red-covered entries correspond to **60 groups at a 0.5 arcsec positional-link threshold**. This is positional QC, **not** evidence of 60 independent brown dwarfs.

## Precision and unresolved selection-area problem

- The post-rasterization 1.2 arcsec averaging and 0.3 arcsec shift alter tested area estimates by at most **0.1054%** and candidate counts by **0**. This is only *integration-grid* sensitivity, **not** an independent resampling of the original 60mas FITS pixels.
- **Systematic caution:** all 20 tiles lose **0.955–0.959%** of their positive-WHT area when bilinear deposition is clipped to 100% occupancy. This systematic per-tile area-conservation deficit is much larger than the 0.1054% integration-grid change. The WHT>0 estimate is **0.9506% below** the exact archival footprint, comparable to the deposition loss; do **not** interpret that difference as physical missing coverage or claim sub-percent absolute precision without a revised conservative rasterization or native-pixel test.
- The independently estimated *usable* F444W area **~0.49672 deg²** is **0.032104 deg² (6.46%) smaller** than WHT>0. Do not adjust thresholds to reproduce it: masks, depth, photometric extraction, and product provenance may differ.
- **Conditional only:** if 0.49672 deg² is truly a subset of this measured WHT>0 footprint, the 365-day red usable-area fraction must lie between **56.88% and 63.34%** (set-theoretic bounds). This does not identify its exact spatial overlap, and the 0.95% deposition bias also affects the bounds.
- The official [COSMOS2025 catalog documentation](https://cosmos2025.iap.fr/catalog.html) exposes `wht_f444w`, `snr_f444w`, `flag_star`, and `flag_blend`; the released JWST starmask and an **explicit calibrated selection/depth threshold** should be investigated. F444W>0 does not encode the F115W/F150W/F277W non-detection selection or cross-filter PM detectability.
- **Next scientific step:** recover a reproducible *selection-matched* spatial mask (depth, bright-star masking and source flags), validate native-pixel area conservation (without conflating MOC overfill with physical coverage), then intersect it with the saved exact epoch polygons. No local workstation is required for the current public-only WHT result. Do not download all 7 GB again without a justified algorithmic fix.

**Recommended wording for Tomo (provisional):** “About 59.5% of the publicly released COSMOS-Web F444W WHT-positive footprint has independent F444W-plus-red NIRCam imaging separated by at least one year (62/125 catalog entries lie in the geometric overlap). This is an imaging-opportunity fraction, not an astrometric completeness or dropout-selection-matched area; the latter is still being calibrated.”
