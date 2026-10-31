# Local v13 F444W WHT spatial QC (2026-10-31)

The user-provided header inventory is summarized in [results/local_wht_header_qc_20261031.json](results/local_wht_header_qc_20261031.json). It confirms a 50,810,063,040-byte, 115593-by-109890 v13 F444W image with native celestial WCS, and **80 COSMOS-Web 30mas compressed WHT tiles**: four filters (F115W/F150W/F277W/F444W), all 20 A/B tiles in each. The COSMOS-Web tile *WHT primary headers* do not contain celestial WCS; do not guess their placement.

The v13 entire **rectangular image bounds** imply approximately 0.882119 deg2 at 30mas. This is **NOT** a positive-weight footprint or selection area. It exceeds the previously measured official 60mas COSMOS-Web WHT>0 geometric area (0.528824 deg2) because the v13 FITS rectangle includes unreduced/unobserved space. Do not interpret this difference as catalogue completeness or differing survey depth.

## Execute on the user's already-cloned Ubuntu checkout

```bash
cd ~/BD_number_density
git pull
source .venv/bin/activate
python -m unittest pm_footprint.test_local_v13_wht_sample

python -m pm_footprint.local_v13_wht_sample \
  --v13 "/run/media/$USER/My Passport/Zhenkai/F444W_psf_matched_F444W_WHT_v13.fits" \
  --stride 120 \
  --out-prefix pm_footprint/results/local_v13_sample
```

The `--stride 120` default samples 3.6-arcsec cells at two offsets (0.25, 0.75 of each cell). It reads **two sets of native data row spans** using read-only `os.pread`, not an entire 51GB array or a memory map. It does not process the 121GB COSMOS-Web gzipped tile collection.

Outputs to **return for scientific review**:
- `pm_footprint/results/local_v13_sample.json`: two provisional WHT>0 area estimates, phase-offset disagreement and counts.
- `pm_footprint/results/local_v13_sample_candidates.csv`: *actual nearest native-pixel* v13 WHT at the frozen 125 candidate centers.
- `pm_footprint/results/local_v13_sample_grid.npz`: compact reusable two-phase boolean spatial grids with WCS metadata. Use this for prospective intersection with exact epoch polygons and the star-mask; **do not** call a 3.6-arcsec cell a native 30mas pixel.

Only transfer these compact outputs. Keep 50.8GB v13 FITS and 121GB COSMOS-Web compressed WHT tiles on the external disk. Do not upload outputs to the public repository unless visually inspected for provenance and obvious data/astrometric anomalies.

## Scientific limitations

1. Two offset grids are **not a rigorous convergence test** and their difference is not a confidence interval. Narrow star masks, faint-source neighborhoods and tile-edge features require finer/resolution tests.
2. The v13 WHT is a *different processing product* from COSMOS-Web WHT. Presence (WHT>0) is only geometric coverage. Depth/amplitude comparisons require variance convention, calibration, noise correlation, PSF/aperture and local empirical noise.
3. Neither v13 WHT nor four-band COSMOS-Web WHT alone provides F356W repeat-epoch detectability.
4. To report an *effective selection area* at F444W=26.2, reproduce Dave's flux/SNR/color and upper-limit selection, source quality flags, and inject/recover or empirically calibrate noise.
5. The current source-mask union is provisional; spatial intersection must use matched coordinates and deduplicated tile-overlap sky area, not subtract an unmatched 0.013 deg2.

No 7GB WHT download or repeated MAST query is required for this task.
