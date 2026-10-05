# COSMOS-Web Brown Dwarf Number Density

Research code for predicting brown-dwarf / ultracool-dwarf counts in the COSMOS-Web field.

## Scientific goals

This repository separates two quantities that must not be mixed:

1. **Intrinsic field counts / surface density**: expected Galactic UCDs in the COSMOS-Web line of sight as a function of spectral type, effective temperature, and apparent magnitude.
2. **Selection-matched counts**: expected subset that would pass a specific brown-dwarf/dropout selection used by Bahareh, Dave, or Ali.

The first target is the intrinsic population. Selection functions are applied only afterward.

## Baseline models

- FC-ENZO / Innala et al. (2026): Juric Galactic thin+thick+halo model, Aganze et al. scale heights, Kirkpatrick et al. local densities, synthetic Sonora spectra.
- Chen et al. (2025): COSMOS-Web empirical benchmark and Ryan & Reid (2016) double-exponential prediction.

## COSMOS-Web baseline

- Center: RA ~150.12 deg, Dec ~+2.2 deg (exact footprint treatment will be added)
- NIRCam: F115W, F150W, F277W, F444W
- Full survey nominal area: ~0.54-0.6 deg^2 depending on definition/data product
- Depth is spatially non-uniform, so final calculations must integrate over depth/coverage zones rather than assume one limiting magnitude.

## Status

Initial audit started 2026-09-26. Do not use preliminary outputs for a paper until the validation tests in `docs/methodology_audit.md` pass.


## Current validated baseline (2026-09-26)

See `docs/results_2026-09-26.md` and `data/preliminary_cosmosweb_teff_counts.csv`.

Key reference points:
- Literature conventional-T benchmark: Ryan & Reid prediction quoted by Chen et al. (2025): 21.4 T0-T5 dwarfs over 0.54 deg2; 13.1 +/- 7.9 over Chen's 0.243 deg2 search area at F115W=27.45.
- Cold-population baseline using official COSMOS-Web depth zones + Elf Owl reach + Kirkpatrick (2024) Teff-bin densities: about 52.5 detectable 450-750 K objects, and >=62.2 in 300-750 K, before project-specific selection/completeness cuts.
- Public FC-ENZO implementation discrepancies/bugs are documented and are not silently inherited.


## Independent UNCOVER validation (2026-10-05)

The Galactic line-of-sight model has now been tested on the independent
UNCOVER/A2744 JWST field using Li et al. (2026).

- Observed: 5 T dwarfs over 53.4 arcmin^2 (0.094 arcmin^-2).
- Current depth-only model over 450-1500 K: 2.88 objects
  (0.0539 arcmin^-2).
- The prediction lies inside the exact 68% Poisson interval for 5 observed
  objects, so the total-count check is encouraging.
- The 450-1050 K distribution is more late-T-rich than the fixed-component
  baseline predicts, motivating a Teff-dependent thick-disk sensitivity test.

See `docs/uncover_validation_2026-10-05.md`,
`data/uncover_validation_2026-10-05.csv`, and run
`python -m src.validate_uncover`.

This remains an intrinsic depth-only validation. The Li et al. colour/SED
selection must be forward-modelled before claiming selection-matched
agreement.
