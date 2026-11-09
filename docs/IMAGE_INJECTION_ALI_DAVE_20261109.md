# Image-level injection/recovery for Dave125 and Ali45

**Status (2026-11-09): partial implementation, NOT a measurement.**
Branch includes tested building blocks for injecting synthetic sources into
calibrated image arrays and for applying each project's photometric cuts.
**No real COSMOS-Web SCI mosaics, aligned multi-band PSFs, HST F814W SCI, or
the original complete source-extraction chain are committed to this repo.**
There is therefore **no measured image-injection completeness** for Dave or
Ali yet. No published number counts or paper candidate flags were changed.

## Why two distinct experiments?

- **Dave (JWST F444W-only):** use the paper's corrected \`new_e_*\`
  signed S/N thresholds, the *positive flux versus 2σ* F277W rule,
  \`F444W/F277W_eff>30\`, and \`m_F444W<26.2\`.
  Code: \`dave_photo_selection\`.
  The current repository independently checks these conditions on
  **125 preselected catalogue rows**, NOT on newly injected sources.
  Paper's ~2 predicted photometric sources remain a model conditional on
  depth-zone selection; they do not imply ~2 measurable proper motions.
- **Ali (JWST+HST):** the initial catalogue preselection needs
  F444W S/N>5 and positive AUTO flux, F814W/F115W/F150W/F277W S/N<2
  and \`warn_flag==0\`. Code: \`ali_photo_selection\`.
  This is **preselection only**, reconstructed from existing project
  context. The final morphology/visual screening and the exact original
  candidate code need to be confirmed with Ali's frozen manuscript.
  Do not declare the 45 (or the subsequently rejected ID 569254)
  reproduced by these functions alone.

Both samples must be processed from the actual **same calibrated,
registered images and noise maps used by their own catalogue pipeline**.
Reusing another project's PSF, aperture, errors or masks will not measure
the published selection completeness.

## Stage 1 implemented in this branch

\`injection_recovery/core.py\` provides:

1. \`inject_psf_njy(SCI, VAR, PSF, x, y, flux_njy,
   electrons_per_njy, rng)\`: preserves scene/background noise, adds
   PSF-convolved point-source flux and the **extra source shot noise**
   to new SCI+VAR arrays; rejects edge-crossing injections and
   invalid units/calibration. All supplied arrays must already be in
   nJy/pixel, (nJy/pixel)^2, and same-grid PSF.
2. \`dave_photo_selection\`: exact checked Dave *photometric* rule,
   using signed, corrected flux/error. It does NOT create \`new_e\`
   from raw images or reproduce his whole catalogue.
3. \`ali_photo_selection\`: explicitly provisional Ali *photometric*
   preselection requiring a real F814W measurement and catalogue flags.
4. Synthetic unit tests in \`tests/test_injection_recovery_core.py\`.

Not implemented here: mosaic-scale image output, detection/deblending,
source extraction, fitted apertures, reproducing the \`new_e_*\`
noise correction on injected sources, cross-band registration, masks,
baseline-candidate rejection, final morphology/quality decisions,
or physical population-weighted number counts.

## Local inputs needed before a physical recovery run

| Input | Dave | Ali | Why |
|---|---|---|---|
| Aligned native F115W, F150W, F277W, F444W SCI | yes | yes | Inject into the *observed* scenes, not a noise-free catalogue |
| Pixel uncertainty/variance plus WHT and valid-pixel masks | yes | yes | Local noise, coverage, injection-safe area |
| Band- and position-matched PSFs plus flux calibration/exposure | yes | yes | Same source shape and source shot noise as extraction |
| Original JWST catalogue, extraction config, \`new_e_*\` convention | yes | as used | Repeat detection, flux and errors consistently |
| Original HST F814W SCI/uncertainty/PSF with verified WCS | no | **yes** | Ali's optical nondetection is part of the published preselection |
| COSMOS2025 master / v13 rows and candidate ID roster | yes | yes | Reproduce original selections and guard against mismatch |
| Original masks, field geometry and per-band co-registration checks | yes | yes | Correct spatial denominators and avoid WCS artifacts |

**Important:** the previously identified SCI/ERR vs F444W WHT header-WCS
discrepancy must be resolved with original i2d WCS/source matches before
placing sources across bands. A shared pixel index is NOT proof of
co-registration.

Avoid uploading many-GB imagery to GitHub. The scripts should run on the
user's workstation with read-only references to local FITS paths.
Run a small tile/patch pilot first, then expand by coverage/depth strata.

## Required next stages (not yet run)

1. Freeze exact Dave and Ali selection functions and source extractor
   against their *full* original catalogues (not only preselected rows).
2. Preflight per-filter SCI/VAR units, PSF flux normalization,
   WCS alignment, masks, exposure/count conversion and aperture methods.
3. Draw injected L/T/Y source spectra and apparent magnitudes using
   the existing Galactic model, with the same intrinsic SED across bands.
   Sample valid sky positions, magnitude, crowding, depth and subpixel
   offsets, with reproducible random seeds and no dense injections.
4. Rerun the actual detection, deblending, forced photometry and quality
   flags **blindly** on injected images. For Dave use the corrected
   \`new_e_*\` error logic; for Ali re-extract F814W and match the
   final morphology/visual eligibility rule. Scoring from injected
   *input* flux rather than recovered flux is not a recovery test.
5. Require reciprocal unique positional matching and no credit for
   a pre-existing candidate or a flux-dominated blend. Preserve all
   failure stages: not covered, not detected, bad photometry,
   colour rejection, morphology/quality rejection.
6. Compare recovered/placed by magnitude, spectral type, sky depth,
   crowding and filter availability. Integrate response with the
   Galactic population and **matched effective area**. Propagate
   binomial/Monte Carlo uncertainty and scientific model systematics.
7. Keep observed 125/45, model-predicted counts, and simulated
   recoverable counts as **different quantities**. No substitution
   of one for another without a documented matched definition.

## Minimal offline smoke test

\`\`\`bash
python -m pytest -q tests/test_injection_recovery_core.py
\`\`\`

This validates injected flux/variance bookkeeping and two preliminary
photometric selectors **only**. It does NOT prove astrophysical completeness.

**Gate before claiming a result:** same-source reproducibility, image
injection+recovery with the actual catalogue pipeline, masks/depth checked,
negative-control positions, and an uncertainty budget.
