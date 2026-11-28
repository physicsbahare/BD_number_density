# Ali depth-limited prediction implementation plan

User specification: follow-up to d431c5e explicitly authorizes a sensitivity-defined detectable-population estimand. Continue inline without additional design approval.

Design: adopt F444W intrinsic aperture signal >= 5 times the calibrated empirical background noise, with PSF response and existing synthetic flux normalization. Apply native signed SNR probabilities within that finite local distance. This is a chosen survey-depth population definition, not a claim that catalog selection has hard support. Keep all model parameters, data products and probabilities unchanged. Test 4/6-sigma sensitivity limits; compare native nominal-SNR boundary as a separate convention diagnostic. No fixed-distance cap.

- [x] Review Ryan & Reid, Aganze, Chen and FC-ENZO primary methods; record differences rather than imply identical thresholds.
- [x] Parameterize only the intrinsic depth threshold; verify default equivalence and inverse-square distance scaling.
- [x] Reuse the completed 512-cell native five-sigma result; execute boundary, ACS-calibration and spatial-convergence checks from cached maps.
- [x] Calculate conditional empirical-CDF sampling uncertainty and retain measured-density errors/incomplete-bin labels.
- [x] Save numerical tables, figures, sources, reproduction commands and Section 4.6 paragraph; run ten focused tests and independent scientific review. Package for the requested GitHub commit.

Outcome: N=12.360844 on0.532430831deg²; measured-density sigma0.841443 and approximate conditionalCDFsigma0.025519. Four/six empirical-sigma definitions give16.085977/9.896737; no fixed-distance cap. Exactcatalog selection and all-distance count remain distinct and unclaimed. See docs/ALI_DEPTH_LIMITED_20261128.md for full qualifiers and results.

Files: minimal optional argument in src/ali_spatial_counts.py; cached-map driver src/ali_depth_limited.py; focused tests; compact output directory data/ali_depth_limited_2026-11-28; scientific report docs/ALI_DEPTH_LIMITED_20261128.md. Existing depth reconstruction and injection scripts are not executed.
