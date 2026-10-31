# Dave 2026-10-29 paper: exact published dropout selection audit

**Source:** user-provided 12-page manuscript `Dave_BD_Paper(1).pdf` (SHA256 `5760a53a617a2cbc9eb7fd90bec748ec6a4b7e5a95117e12198b04fb405a5084`), uploaded `local_v13_catalog_idcheck_candidate_rows.csv` (SHA256 `93400baee3daccfc5020b01f93d05e9b65ef238afa62b5b691f4aa90728e2a40`, 125 rows independently position/ID matched to the 3,280,023-row Zhen-Kai v13 FITS), uploaded alignment and ID-check JSON. Only compact audited evidence is committed, not the original manuscript or catalogue.

## Clarification: non-detection is not lack of coverage

The paper's Section 2.2 **explicitly selects S/N<2 in F115W, F150W and F277W**, plus S/N>5 in F444W. Figure 1 illustrates this photometric non-detection. The previously reported **78/125 sources outside image bounds** obtained by mapping sky positions through blue SCI/ERR **header WCS** is a **coordinate-registration artifact**. Their blue headers are discrepant by ~67,565 pixels from the F444W WHT WCS. The user's new local alignment diagnostic instead finds **125/125 positive blue ERR values at the common F444W native-pixel locations** in all three blue bands (and F444W), as assumed by Zhen-Kai's forced-aperture code. This confirms valid ERR pixels at common indices, not source detection and not proof of correct image co-registration; compare nearby bright-source morphology in SCI before using the shared-index depth maps scientifically. Do not label the 78 as absent blue imaging.

## Published cuts reproduced on all *125 preselected* source rows

Paper Section 2.2 (pp. 2-3):
1. Signed `f_F444W/new_e_F444W>5`; signed `f_F115W/new_e_F115W<2`, `f_F150W/new_e_F150W<2`, `f_F277W/new_e_F277W<2`.
2. `f_F444W/f_F277W,eff>30`, where **positive** `f_F277W,eff = f_F277W`, and **zero/negative** `f_F277W,eff=2*new_e_F277W`. **Never** replace positive F277W with a 2sigma upper limit.
3. `mag_F444W<26.2 AB`.

The exact recipe selects **125/125** from the frozen 125-row export, with all separate cuts also 125/125. **99 F277W fluxes are positive; 26 are nonpositive**, and only the latter 26 use a 2sigma floor. Minimum corrected F444W S/N is **5.3593**, maximum corrected F277W S/N **1.97335**, minimum hybrid flux ratio **30.19163**, maximum mF444W **26.1980 AB**. Using *uncorrected* `e_*` produces just **87/125** all-cuts passing (F277W S/N<2 only 95/125), while the source notebook explicitly derives empirical corrected `new_e_*`. This compellingly resolves the previously ambiguous noise-error convention, though Dave's original script is still absent.

**Nuance:** the paper writes `S/N<2`, **not `|S/N|<2`**. The frozen sample includes one F115W strongly *negative* S/N (ID 142367, ~-5.63) and one F150W slightly below -2 (ID 203707, ~-2.008). Both pass the signed criterion; inspect the negative extreme for residual background artifacts in the source SCI rather than silently changing the selection.

**Quality flags:** all 125 `use_phot=1`; all have `combined_artifact_flag=0`. Yet 85 have nonzero `FLAG_AUTO_LW` (81 at 32; 4 at 128), and 60 have `FLAG_BAD_KRON_LW=1`. Do not impose blanket `FLAG_AUTO_LW==0`/Kron-good filters: that would remove selected sources without reproducing the paper. Aperture diameters: 112 at 0.32 arcsec, 6 at 0.48 and 7 at 0.70.

**Reproducible source:** [validate_dave_paper_selection.py](validate_dave_paper_selection.py) and [test_validate_dave_paper_selection.py](test_validate_dave_paper_selection.py); provenance and results in [results/dave125_paper_selection_reproduction_20261031.json](results/dave125_paper_selection_reproduction_20261031.json). On the user's laptop:
```bash
cd ~/BD_number_density
git pull
source .venv/bin/activate
python -m pm_footprint.validate_dave_paper_selection \
  --rows pm_footprint/results/local_v13_catalog_idcheck_candidate_rows.csv
```
If the CSV was saved elsewhere, supply its actual path to `--rows`.

### IMPORTANT: paper already contains a 0.5304 deg2 catalog-derived area

Section 6 (p.8) specifies **joint four-band depth zones using corrected errors** and gives their *catalog-derived geometric area* as **0.5304 deg²**. It then integrates selection probability over distance, Galactic component and depth zone, yielding a **model-predicted ~2 photometrically selected ultracool dwarfs**, **not** ~2 necessarily PM-confirmable objects. The local v13 F444W WHT>0 sampled geometry **0.533619 deg²** is only **0.003219 deg² (0.603%)** larger, but those footprints have different definitions and are not yet spatially matched. Neither number is PM effective area. We should **not present 0.5304 as a new measurement or claim its precise construction has been independently validated**.

Section 3 (pp.3-4) distinguishes 387132's multiple F444W/F410M exposures and more robust displacement from 1651050's sole usable F356W exposure, PSF sensitivity and parallax ambiguity, requiring independent epoch. Do not call both equally secure intrinsic-PM determinations.

## Next *original-catalog* validation: full 3.28-million-row stage counts

Paper explicitly states first SNR stage **2,997**, ratio stage **335**, magnitude stage **125**. The frozen 125 subset cannot test those counts. [v13_full_catalog_paper_stages.py](v13_full_catalog_paper_stages.py) now performs a **single read-only sequential ~6.57GB pass** (20,000 rows / ~40MB chunk) over original local FITS BINTABLE and reports those stages under `new_e_*` versus `e_*`, with and without `use_phot==1`, plus actual ID overlap with the frozen roster. It does not read any 50GB images or download anything.

```bash
cd ~/BD_number_density
git pull
source .venv/bin/activate
python -m pm_footprint.v13_full_catalog_paper_stages \
  --catalog "/run/media/$USER/My Passport/Zhenkai/F444W_catalog_update/cosmos_F444W_v13_catalog_v1.fits"
```

Return just `pm_footprint/results/v13_full_catalog_paper_stages.json`. Investigate discrepancies between counts before modeling completeness; do not fit thresholds by eye to force 2997/335/125.

## Remaining real science limitations

- Paper's `0.5304 deg²` is a depth-zone area, not a demonstrated source-recovery completeness to 26.2 AB nor a spatial mask for the >=365day PM geometry.
- Joint 4-band source detection/depth, selection probability, crowding and the **already supplied and processed** JWST star-mask need true spatial matched masks; approximate ratios cannot be promoted to full effective area.
- Actual F356W exposures or forced-photometry sensitivity and centroid errors are required to predict **measurable** PM area among the 62 archive-geometrically eligible sources.
- For F277W-positive flux sources, Dave uses the measured flux, not a 2sigma upper limit; a field-wide 5sigma F277W>28.9 AB cutoff would impose a different, stronger selection. The 28.9 conditional threshold only applies to demonstrating >30 for an approximately zero-flux upper-limit case at F444W=26.2.
