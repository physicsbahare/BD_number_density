# 2026-11-01: F444W-only quality-flag conflict and archive-clipped WHT-mask estimator

## New 125-object measurement (frozen v13 fixture, **not** official COSMOS2025 flags)

The official [COSMOS2025 catalog documentation](https://cosmos2025.iap.fr/catalog.html)
defines `warn_flag=4` as a significant detection in **only one NIRCam band**
(`delta_mag<0.2`), with other NIRCam+HST+MIRI measurements `delta_mag>0.5`
or undetected. This is a *potential conflict* with intentionally selected
F444W-only brown dwarfs: rejecting all `warn_flag=4` sources would be a
selection-dependent veto, not a neutral quality cut.

For a positive-flux source with approximately Gaussian small-error photometry,
`delta_mag ≈ (2.5/ln 10)/(F444W / sigma_F444W)`, so `delta_mag<0.2`
corresponds to `S/N_F444W > 5.428681`. Dave's `S/N>5` cut overlaps this.

Reproducible audit of the *already committed* 125-row Zhen-Kai v13 fixture,
using `new_e_F444W` (the corrected aperture error):

| Diagnostic | Count |
| --- | ---: |
| F444W approximate `delta_mag<0.2` | **124 / 125** |
| 365d red-epoch opportunities satisfying that condition | **62 / 62** |
| Outside 365d red epoch satisfying that condition | **62 / 63** |
| Only exception, 5 < S/N <= 5.42868 | **ID 9095** |
| PM-confirmed IDs satisfying F444W half-condition | **387132 and 1651050** |
| Signed blue S/N below -2, despite passing `S/N<2` | **2 IDs** |

The two negative-blue outliers are ID **142367**, F115W `S/N=-5.6303`,
and ID **203707**, F150W `S/N=-2.0082`. Negative flux is not itself a
detection; these deviations warrant local background/error QA, **not**
automatic removal. Signed `S/N<2` is not the same as `abs(S/N)<2`.

**Crucial limitation:** These 124 are **NOT** 124 measured
`warn_flag=4` entries. The official flag also requires checks of all other
bands and is derived from a **different** COSMOS2025 catalog. The v13 fixture
does not contain those official flags or all HST/MIRI data. No candidate is
rejected or reclassified here. Official `flag_star` is based on the fraction
of source pixels overlapping the JWST starmask; a zero center-pixel mask
hit for 125 objects cannot establish zero aperture/segmentation overlap.
`flag_blend` flags large F444W Kron corrections and needs independent review.

Reproduce with
`python -m pm_footprint.dropout_quality_flag_overlap`.
The frozen numerical result is
`pm_footprint/results/dropout_quality_flag_overlap_20261101.json`.
Tests: `python -m unittest pm_footprint.test_dropout_quality_flag_overlap -v`.

## Corrected direct mask × WHT × epoch estimator

The existing sparse v13-WHT routine originally counted `WHT>0` mask samples
**even outside** the COSMOS-Web archival F444W survey polygon, despite using
reference areas that were already clipped to that polygon. This can
over-subtract excluded area from the denominator. The estimator now defines

`M_W = mask ∩ WHT_positive ∩ archive_F444W` and
`M_R = M_W ∩ independent_red_epoch_365d`.

The estimator explicitly reports how many WHT-positive samples were outside
the archive and ignored. Two synthetic regression tests cover the mismatch.

**Scale diagnostic, not a corrected sky measurement:** In the old 1000-draw
mask-grid sample, 995 were inside the archive. If all five outside samples
were WHT-positive, the un-clipped estimator could over-subtract at most
`5/1000 × 0.0130711831 = 0.0000653559 deg²`. For the frozen v13 reference
areas, this corresponds to a worst-case red-fraction difference of
approximately **0.00765 percentage points**, under the fixed old raster
and area assumptions. The old mask sampling WCS has since been found
incorrect, so this bound is **only an illustration of the estimator bug**,
not a validated systematic uncertainty for COSMOS-Web.

## What remains unmeasured

1. **Official DR1 WHT** and JWST star-mask unique-sky intersection, with
   corrected PC/TAN WCS and native-pixel area conservation. The 0.0128–0.0131
   deg² previous union estimate is provisional after the WCS bug.
2. A persisted compact spatial WHT grid is missing. Previous WHT summary
   JSON contains area totals, not cell coordinates, so exact spatial
   intersection **cannot** be reconstructed from those totals alone.
   The existing `jwst_mask_wht_sparse_intersection.py` is for the
   **independent local v13 WHT**, not the official 20-tile DR1 WHT.
3. Spatially varying F444W=26.2 AB dropout selection completeness, local
   corrected F115W/F150W/F277W errors, `flag_star`/`flag_blend`, neighbor
   crowding, and F356W detection/centroid precision. None is inferred from
   the geometric 59.495% result.

**Valid current official result:** 0.314624457 / 0.528823685 deg² =
**59.495% geometric** F444W+red coverage at >=365 days, with **62/125**
candidate entries (60 proximity groups). It is **not** a PM-completeness
fraction. No WHT re-download, MAST re-query, or replacement starmask was
used for this audit.
