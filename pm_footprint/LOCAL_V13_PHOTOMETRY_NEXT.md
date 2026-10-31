# Next local data step: Zhen-Kai v13 four-band photometry preflight

**Status: source file locations identified by the user, not yet analyzed at the pixel or catalog-column level.** The user has the original uncompressed `cosmos_F444W_v13_catalog_v1.fits` (~6.6 GB), PSF-matched native 30mas ERR/SCI mosaics for F115W, F150W, F277W, F444W, F444W WHT, `Encircled_Energy_LW_ETCv2.txt` (~5kB), and a ~355kB `JWST_photometry_v13_F444W_detection_NTHU.ipynb`.

## Input inventory from screenshots

- `/run/media/$USER/My Passport/Zhenkai/F444W_catalog_update/cosmos_F444W_v13_catalog_v1.fits`
- `/run/media/$USER/My Passport/Zhenkai/{F115W,F150W,F277W,F444W}_psf_matched_F444W_ERR_v13.fits`
- `/run/media/$USER/My Passport/Zhenkai/{F115W,F150W,F277W}_psf_matched_F444W_SCI_v13.fits`
- `/run/media/$USER/My Passport/Zhenkai/F444W_psf_matched_F444W_SCI_v13_bkgsub.fits`
- `/run/media/$USER/My Passport/Zhenkai/F444W_psf_matched_F444W_WHT_v13.fits`
- `/run/media/$USER/My Passport/Zhenkai/Encircled_Energy_LW_ETCv2.txt`
- `/run/media/$USER/My Passport/Zhenkai/F444W_catalog_update/JWST_photometry_v13_F444W_detection_NTHU.ipynb`

## How to run on Ubuntu

```bash
cd ~/BD_number_density
git pull
source .venv/bin/activate
python -m unittest -v pm_footprint.test_local_v13_photometry_probe

python -m pm_footprint.local_v13_photometry_probe \
  --root "/run/media/$USER/My Passport/Zhenkai" \
  --catalog "/run/media/$USER/My Passport/Zhenkai/F444W_catalog_update/cosmos_F444W_v13_catalog_v1.fits" \
  --notebook "/run/media/$USER/My Passport/Zhenkai/F444W_catalog_update/JWST_photometry_v13_F444W_detection_NTHU.ipynb" \
  --ee "/run/media/$USER/My Passport/Zhenkai/Encircled_Energy_LW_ETCv2.txt" \
  --out-prefix pm_footprint/results/local_v13_photometry_probe
```

Small outputs **to upload privately into ChatGPT** (do not automatically publish to the GitHub repository):
- `pm_footprint/results/local_v13_photometry_probe.json`: v13 catalog FITS schema, names of all flux/ERR/flag fields, calibration metadata and 4-band FITS WCS sanity checks.
- `pm_footprint/results/local_v13_photometry_probe_err_candidates.csv`: 125 frozen Dave source centers and 3x3 native-pixel ERR samples (no source catalog flux or depth inferred yet).
- Optionally upload the ~355kB notebook and ~5kB encircled-energy table so the extraction, aperture/noise conventions can be verified independently **without executing** the notebook or publishing user-specific paths.

The probe reads FITS **headers** and tiny sections of four ERR mosaics using read-only `os.pread`. It does **not** load the ~6.6GB catalog data rows, 9×~50GB images, or downloaded COSMOS-Web WHT tiles; it does not modify the originals.

## Scientific interpretation, not just engineering

1. **Primary dataset for Dave:** his original candidate selection must be reconstructed against the Zhen-Kai v13 catalog (not the independent COSMOS-Web photometry). Identify exact columns corresponding to F444W flux/error, F115W/F150W/F277W non-detections, flux-ratio >30, negative or zero F277W flux and 2σ upper limits. Never silently impose `warn_flag==0`, `flag_star==0` or any flag not verified in his original selection.
2. **Photometric uncertainty:** a PSF-matched *pixel* ERR map is not an `AUTO` or aperture-flux uncertainty. BUNIT, pixel-to-flux calibration, correlated noise and aperture/encircled-energy corrections must be checked against catalog and notebook; measure empty-sky apertures on SCI and ERR maps to cross-check. Large-scale masked-area summaries based solely on WHT>0 remain geometric.
3. **F277W color discrimination:** `F444W/F277W >30` means `m277-m444>3.6928`. At F444W 26.2, verifying this at 2σ from a nondetection requires 5σ F277W depth ~28.90 AB (aperture-consistent), not just any positive WHT. A noisy measured flux ratio is not the same test.
4. **4-band area:** once original Dave selection and flux/error units are verified, construct spatially resolved detection/non-detection/completeness probability, jointly across all four bands and image masks; report the selection-weighted area and its uncertainty, comparing to the COSMOS-Web-depth-based proxy (with truly aligned 30mas SCI WCS).
5. **Proper motion:** there are no F356W epoch images in the user-supplied Zhen-Kai screenshot inventory. The previously measured 62/125 F356W geometric epoch opportunities do not establish F356W detection S/N or astrometric PM sensitivity. Separate F356W epochs or forced photometry will be required later.

Already-computed survey geometry from v13: WHT>0 ~0.533619deg2; WHT∩red >=365d ~0.317401deg2 (~59.48%); conditional code-1 starmask-corrected 59.50–59.54% **geometry only**. The original JWST starmask archive was already supplied/processed and does not need to be provided again.

## QC safeguards

- Candidate roster byte-for-byte SHA256 verified before any read.
- Checks every image WCS and pixel offsets from the v13 WHT; reports mismatches rather than assuming all mosaics share a grid.
- Source-centered ERR sample points have selection bias, so never infer field-wide usable area from their distribution alone.
- The small notebook is scanned as JSON for keyword positions only; no code, private filenames, or credentials from cells are copied to the probe output.
- No heavy GitHub workflow is triggered; synthetic tests only. If any image/cat format is unexpected, fail clearly rather than reading giant arrays.
