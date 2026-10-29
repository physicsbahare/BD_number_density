# Scientific protocol: from WHT geometry to a selection-matched area

**Status (2026-10-30): design only; NO selection-completeness area has been measured.**

## What is already measured

- Exact PUBLIC-only JWST NIRCam archive polygons and separate-epoch overlaps at 30/180/365/730 days.
- Official COSMOS-Web DR1 F444W 60mas WHT>0 approximate geometric coverage (20 tiles); **do not rerun the 7 GB workflow** merely to check status.
- A 365-day F444W+red geometric opportunity of **0.314624457 deg²** within **0.528823685 deg²** WHT>0, or **59.495%**, and 62/125 source entries. The WHT rasterizer has a documented systematic native-to-grid area loss of about **0.96%**, so these are not sub-percent-accurate physical areas.

## Why there is not one universally correct 26.2 AB usable-area mask

For a *specified* source F444W magnitude m, brown-dwarf SED/model s, sky position x and measurement procedure p, define

    C_drop(x; m,s,p) =
      P(F444W detection passes threshold AND F115W,F150W,F277W
        non-detections pass threshold AND adopted F444W/F277W
        color/flux-ratio criterion passes AND quality flags pass
        | x,m,s,p).

This probability varies with source flux, SED and local depth. A binary mask from WHT>0, or a per-tile WHT relative-median cut, is **not** a substitute. If the 26.2 AB selection is described as an F444W *limiting magnitude*, it must be paired with a specified detection S/N definition and noise/aperture correction; do not equate an individual catalog object's native-resolution aperture SNR with the local limiting depth. Dave's original 125-candidate selection may use *auto-flux/error* rather than the catalog's native-resolution aperture SNR column. Verify the exact selection code before calibrating the noise model.

For an assumed source-population distribution p(m,s), a well-defined **effective dropout selection area** is

    A_eff_drop = integral_over_sky [ integral C_drop(x;m,s,p) p(m,s) dm ds ] dOmega.

For a fixed source model and magnitude, use A_eff_drop(m,s) without averaging. Injection/recovery with the same source extraction and selection is the preferred calibration. As a less complete approximation, construct a calibrated 5-sigma limiting-magnitude map from noise/ERR and empirical random-aperture noise, apply an explicit source model and dropout thresholds, and report the approximation as such.

For proper-motion **measurement opportunity** with epoch requirement T and filter policy k, define G_T,k(x) as the previously measured binary multi-epoch geometric mask. The relevant denominator and numerator are

    f_covered(T,k;m,s) =
      integral G_T,k(x) C_drop(x;m,s,p) dOmega
      / integral C_drop(x;m,s,p) dOmega.

This is a *selection-weighted geometrical overlap*, **not** the probability of measurable proper motion. To estimate PM detectability, additionally model/measure both-epoch detection probability, positional uncertainty, astrometric registration, time baseline and a PM significance threshold. Cross-filter detections of F444W-only red dropouts can be especially incomplete.

## Public-data acquisition plan (avoid local workstation first)

1. **Candidate quality and catalog version QC.** Query the public CDS/VizieR J/A+A/704/A339/phot table for the **125 known IDs only** (separate cosmosweb_catalog_qc.yml workflow). Check 0.5-arcsec coordinate agreement with frozen Dave125 positions, compare F444W aperture SNR/weights and blue-band values, and record flag-star, flag-blend, and quality columns when present. The CDS publication copy may not match the updated COSMOS2025 v1.1 catalog or Dave's selection catalog. Do not merge IDs without astrometric validation.
2. **Obtain the official spatial JWST starmask.** The release homepage states the JWST and HSC starmasks are downloadable under catalog/download, but the current download gateway presents a registration/captcha form. Do not invent a static file URL. If obtainable, inspect its coordinate system and semantics and rasterize *mask polygons*, not catalog flag-star values. Also inspect other exclusion masks and the original Dave selection cuts.
3. **F444W 26.2 depth.** Determine the actual PSF/aperture-based 5-sigma threshold from the official ERR/WHT + empirical noise methodology, *not* from arbitrary relative WHT fractions. Verify whether 26.2 is the catalog selection magnitude limit or merely a model truncation. The released DR1 WHT has spatial exposure weights, not an absolute point-source completeness curve.
4. **Blue non-detections.** Include F115W/F150W/F277W depth maps or source injection/noise simulations for brown-dwarf SEDs; a deeper blue image can *reduce* the chance of meeting a non-detection criterion for a fixed SED, so the joint selection footprint cannot be derived from F444W depth alone.
5. **Critical COSMOS2025 quality-flag selection trap.** The official [catalog documentation](https://cosmos2025.iap.fr/catalog.html) defines `warn_flag=4` as a *single-NIRCam-band detection*, classified as potentially instrumental (17,520 entries in the released catalog). Dave's **F444W-only** dropout selection can intentionally produce that same photometric pattern. Therefore, **do not automatically require `warn_flag=0` or remove `warn_flag=4` without inspecting the exact original selection and independent image/PM evidence**: doing so could systematically exclude the very population under study. Conversely, do not assume every single-band candidate is genuine. The published COSMOS2025 source detection is performed on a multi-band detection image, with both hot and cold modes; its `warn_flag` is not a spatial area mask. Keep `flag_star`, `flag_blend`, `warn_flag` and independent imaging checks as separate, source-level QC outputs. The same documentation says `flag_star` encodes a 0–80% overlap of source pixels with the JWST star mask and `flag_blend` identifies large F444W Kron aperture corrections (>3). Neither can be converted into a masked sky area by counting flagged sources.
6. **Crowding and quality.** The 125 entries form 115 positionally linked groups within 0.5 arcsec (62 red-covered entries form 60 groups). This is not an automatic duplicate-removal rule. Inspect flag-blend, flag-star, segmentation, and visual cutouts before treating each as an independent source.
7. **Area/raster precision.** Correct or independently bound the ~0.96% native-to-grid area deposition/clipping loss before claiming a precise WHT>0 denominator. Post-rasterization grid coarsening alone is not independent native 60mas convergence. Reprocess original 7 GB **only** if the revised algorithm has passed synthetic native-pixel area conservation tests and a clearly justified need for higher precision.
8. **Provenance.** Keep the official DR1 WHT mask, any alternate PSF-matched v13 WHT, and any 0.49672 deg² usable-area estimate separate unless their selection, WCS, processing and sky masks are demonstrated equivalent.

## What can be said now to Tomo

At >=365 days, 0.3146 deg² (approximately 59.5% of WHT>0 geometric F444W coverage) has an independent F444W+red NIRCam observation, and 62/125 catalog entries lie in that geometric overlap. These are *opportunities*, not PM-confirmed objects or guaranteed PM measurability. 2/125 is not a coverage fraction. The **usable F444W=26.2 dropout-selection-matched area is not identifiable** from the available geometric mask alone. Under the explicit and unverified assumption that the external 0.49672 deg² usable area is a subset of the measured WHT>0 region, its 365-day red-overlap fraction is bounded between **56.88% and 63.34%**; these are conditional set-theory bounds, not a point estimate.

Sources: [COSMOS2025 catalog documentation](https://cosmos2025.iap.fr/catalog.html), [CDS VizieR photometry table](https://vizier.cfa.harvard.edu/viz-bin/VizieR-3?-source=J/A%2BA/704/A339/phot), [official DR1 NIRCam WHT](https://cosmos2025.iap.fr/nircam.html).
