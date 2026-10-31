# Zhen-Kai v13 photometric selection: source-driven 2026-10-31 audit

## What is verified from the user's uploaded files

The original `JWST_photometry_v13_F444W_detection_NTHU.ipynb` (68 cells, SHA256 `1047f875cd5476804ab9c71e639468eb4abe3d067ffc87ad9cd75cd3b0a97733`) and `local_v13_photometry_probe.json` show that the provided 3,280,023-row, 285-column F444W-detection catalogue contains both:
- `f_F115W,e_F115W` through `f_F444W,e_F444W` in nJy;
- an additional `new_e_F115W` through `new_e_F444W` produced in notebook cell 63 from `e_band * corr_factor[use_aper]`.

For aperture **DIAMETERS** 0.32,0.48,0.70,1.0,1.4 arcsec, notebook cell 62's empirical correction factors are:
- F115W 1.48,1.71,2.00,2.37,2.72
- F150W 1.49,1.75,2.06,2.45,2.91
- F277W 2.74,3.15,3.58,4.08,4.62
- F444W 2.42,2.68,3.00,3.37,3.63.

Thus the *same* flux has **corrected S/N reduced by 2.42–3.63 in F444W**, depending on aperture. A corrected-F444W `SNR>5` criterion would correspond to raw `SNR>12.1–18.15`, so the choice of `e_` versus `new_e_` materially changes selection. These factors apply as implemented by this particular notebook; whether Dave used them is **NOT YET VERIFIED**.

**Important:** `use_phot` is created before cell 63. In this notebook it combines a F444W 0.32 arcsec aperture `SNR<3` low_snr flag, a bad-pixel flag and an artifact-near-bad-pixel flag. It is not a generalized proof of 5sigma F444W dropout detectability; blindly applying `use_phot==1` may change the roster, so Dave's original cuts must be checked.

The catalogue header exposes aperture-flux, AUTO, adaptive aperture, per-band native and revised errors, as well as `FLAG_BAD_KRON_LW` and `FLAG_AUTO_LW`; names `flag_star` and `flag_blend` were **not found as named columns** in this particular F444W-detection catalogue. This does not mean their physical quality issues are absent.

### Critical WCS discrepancy

All four SCI and ERR images have identical 115593 × 109890 array shapes. In the uploaded preflight, **the blue F115W, F150W, F277W SCI/ERR celestial WCS headers put the same frozen Dave sky positions ~67,565 native pixels (~33.8 arcmin) from the reference F444W WHT pixels**. F444W ERR matches the F444W WHT WCS and F444W SCI differs by only ~0.072 pixel. The initial WCS-based sampling produced 78 out-of-bounds candidate centers and only 33 positive ERR values in each blue band. **Neither number may be interpreted as a real blue-image coverage fraction or non-detection.** In notebook cell 42, the blue SCI/ERR images are processed at the *same F444W segmentation pixel positions* without any blue-WCS coordinate transform. This points to either stale/mismatched image headers or a reduction registration issue. Do not accept astrophysical conclusions until image-grid alignment is tested.

The user-uploaded F444W PSF file is a 201×201 normalized float32 image; its total pixels sum ~0.99999982. The user-supplied EE table gives at **radius 0.16 arcsec** F115W=0.8168, F150W=0.8034, F277W=0.7646, F444W=0.6446. The notebook normalizes the empirical F444W PSF curve-of-growth to EE at radius 1 arcsec, so the EE table alone is not a substitute for reproducing exact per-source aperture corrections.

At F444W=26.2, the flux is 120.226 nJy; `F444W/F277W>30` would require `F277W<4.0075 nJy`. For a **2sigma nondetection upper limit assuming zero measured flux**, that requires local sigma277 ≤2.0038 nJy, corresponding to an aperture-consistent 5sigma F277W depth ≥28.898 AB. These are conditional thresholds, **not measured completeness**.

## Small reproducible local diagnostics

The user's original 50 GB maps and 6.6 GB catalogue stay on their disk. Update the repository and use the existing environment:

```bash
cd ~/BD_number_density
git pull
source .venv/bin/activate
python -m unittest -v pm_footprint.test_local_v13_selection_probes

python -m pm_footprint.local_v13_alignment_audit \
  --root "/run/media/$USER/My Passport/Zhenkai"

python -m pm_footprint.local_v13_catalog_idcheck \
  --catalog "/run/media/$USER/My Passport/Zhenkai/F444W_catalog_update/cosmos_F444W_v13_catalog_v1.fits"
```

Return **only these small outputs** in chat:

- `pm_footprint/results/local_v13_alignment_audit.json`
- `pm_footprint/results/local_v13_alignment_audit_centers.csv`
- `pm_footprint/results/local_v13_catalog_idcheck.json`
- `pm_footprint/results/local_v13_catalog_idcheck_candidate_rows.csv` **if created**.

The alignment audit reads only 3×3 ERR pixel patches at the 125 reference-F444W image positions versus sky-WCS-mapped positions and indicates whether samples differ. **Positive ERR at shared pixels is necessary, not sufficient, to demonstrate correct astrometric co-registration:** the subsequent science step should use bright source morphology in the SCI maps or trusted WCS from aligned originals.

The catalog ID-check reads the BINTABLE header and **at most 125 sparse 2003-byte rows**, each at a hypothesized row index ID−1, and validates both FITS row ID and sky RA/Dec (<0.5 arcsec). This is deliberately not assumed to be an ID-matched catalogue. If any IDs/coordinates disagree, results will be marked unusable for source photometry; then build a bounded-memory *spatial* crossmatch through the 6.6 GB catalogue, rather than coercing unrelated IDs. It does not scan the full 6.6 GB file or attempt selection completeness.

Even if the 125 match, returned scenario results only test a **preselected subset** under several choices of `e` vs `new_e` and negative F277W handling; they are NOT the full 3.28-million-source selection function.

## Remaining requirements for a publication-ready selection-matched area

Recover Dave's *original selection implementation* and check negative flux handling, catalogue ID equivalence, aperture used, error convention, quality flags and contamination. Verify aligned blue SCI/ERR WCS and empirical correlated noise over spatial cells; then inject synthetic F444W-only and realistic SED sources to estimate a four-band selection area at F444W=26.2 and matched star exclusion. The separate F356W epoch detection and centroid precision remain necessary to translate geometric area into measurable PM completeness.

Current geometric overlap remains v13 WHT>0 ∩ >=365d archival red observations ≈0.31740 deg², not dropout selection or PM completeness. The JWST mask user archive has already been analyzed; do **not** ask for it again or rerun high-cost WHT.
