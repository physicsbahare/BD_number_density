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
