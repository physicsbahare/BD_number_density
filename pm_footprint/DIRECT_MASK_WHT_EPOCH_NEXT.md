# Direct JWST star-mask × v13 WHT × F356W epoch spatial intersection

**2026-10-31 status:** implementation and offline synthetic tests added. The final **direct spatial** result is **not yet measured** because the original 19,370,234-byte JWST star-mask archive is not mounted in this particular automation runtime, and the user's 50.8GB v13 WHT is on their external drive. The archive was ALREADY provided and processed in earlier analysis. **Do not ask for or redownload it again.** The committed archive-only Monte Carlo result (mask∩archive-red ~0.007516 deg²) and conditional WHT set bounds are separate.

## Efficient local run (no WHT download or full FITS loading)

The existing code `jwst_mask_union_repro` reads the **existing user-supplied** ~19MB compressed 20-tile JWST starmask, creates a 0.6-arcsec equal-area clipped-sum mask-union raster, and draws weighted Lambert cell centers. It does not save the full raster. The new script then queries the **original local** 50.8GB WHT at those points with 4-byte read-only `os.pread` calls and intersects them with the already-frozen 717 PUBLIC red imaging CAOM polygons/times.

```bash
cd ~/BD_number_density
git pull
source .venv/bin/activate

# Set to the EXISTING archive's local path, not a new download.
MASK_TAR="/path/to/your/existing/cosmos_web_starmask_jwst.tar.gz"
test -s "$MASK_TAR"

python -m pm_footprint.jwst_mask_union_repro \
  --archive "$MASK_TAR" \
  --output pm_footprint/results/jwst_starmask/local_mask_union_20k.json \
  --sample-out pm_footprint/results/jwst_starmask/local_mask_weighted_20k.csv \
  --samples 20000 --seed 20261031

python -m pm_footprint.jwst_mask_wht_sparse_intersection \
  --points pm_footprint/results/jwst_starmask/local_mask_weighted_20k.csv \
  --mask-summary pm_footprint/results/jwst_starmask/local_mask_union_20k.json \
  --wht "/run/media/$USER/My Passport/Zhenkai/F444W_psf_matched_F444W_WHT_v13.fits" \
  --metadata pm_footprint/results/mast_observations.csv \
  --out pm_footprint/results/jwst_starmask/direct_v13_wht_epoch_intersection.json
```

Upload **only** `local_mask_union_20k.json` and `direct_v13_wht_epoch_intersection.json` (and optionally the small `local_mask_weighted_20k.csv` for independent point-level QC); not the 50GB FITS, not the archive again. The new sparse WHT reader fetches at most `20000 × 4 = 80000` bytes of pixel values plus FITS header, although HDD random-access seeks can take time. The mask-raster creation is the main CPU/memory work, but its input is only the existing ~19MB tar.gz. The earlier 7GB 60mas COSMOS-Web WHT GitHub workflow must remain **manual only**.

## What is new versus previous bounds?

The earlier `epoch_mc_20261031.json` directly intersected mask samples with **archival polygons**, but NOT with WHT. Its `mask_within_archive_red_365d ≈ 0.007516 deg²` must NOT be subtracted from `WHT∩red` without spatial overlap. The new script measures, with one shared weighted sample:

- `|mask ∩ WHT|` directly from local native WHT pixels;
- `|mask ∩ WHT ∩ F444W+red >=365d|` directly from exact CAOM epoch membership;
- `(WHT_red - mask_WHT_red)/(WHT - mask_WHT)` and joint sampling-only uncertainty.

For v13 the reference two-phase sampled `WHT>0 ≈0.53361911 deg²`, `WHT∩red365≈0.31740121 deg²` (geometric). Those **3.6-arcsec** WHT reference areas have edge aliasing not represented by the reported MC 95% halfwidth. Mask 0.6-arcsec union lower/upper approximation (~0.01286–0.01307 deg²) is likewise a systematic. Do not report a high-precision survey area until higher-resolution convergence and mask semantics are checked.

## Separate F444W=26.2 AB selection science

The supplied original v13 notebook's `new_e_band=e_band × correction(use_aper)` uses factors F444W 2.42–3.63, F277W 2.74–4.62. Using nominal `e_*` instead of corrected `new_e_*` overstates nominal 5sigma depth by **0.96–1.40 mag (F444W)** and **1.09–1.66 mag (F277W)** if `new_e_*` is the intended noise model. See `results/v13_psf_noise_qc_20261031.json`. The independent F444W PSF 0.03-arcsec curve of growth normalized at 1 arcsec agrees with the table within 2.14% at tested apertures. Dave's actual use of `e_*` vs `new_e_*`, negative F277W ratio logic, source-quality flags, source matching and F356W detectability are still unverified. **The direct geometric mask result, once obtained, is still not dropout completeness.**
