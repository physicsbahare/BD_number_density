# COSMOS-Web Brown Dwarf Number Density

Research code for predicting brown-dwarf / ultracool-dwarf counts in the COSMOS-Web field.

## Ali F814W real-image continuation (2026-11-28)

The follow-up to `5c8cd56` adopts native catalog SNR and a valid five-band
geometric measurement footprint of **0.532430831 deg²**. The unchanged
untruncated halo and nonzero source-free forced-photometry selection probability
make the all-distance integral divergent. **No finite final paper count or
surface density is accepted.** A source-associated catalog detection kernel or
an independently justified population domain is required. The historical
0.496306 deg² mask remains unidentified; no area was tuned to recover it.

See the [executed validation resolution](docs/ALI_VALIDATION_RESOLUTION_20261128.md),
[new numerical tables and figures](data/ali_validation_2026-11-28), and
[acceptance status](data/ali_validation_2026-11-28/result_status.json).
Completed F814W calibration/injections, depth maps and atmosphere calculations
were reused. The previous 23.8437 result remains an explicitly conditional
diagnostic, superseded as a candidate full native-selection prediction.

See the [executed results and limitations](data/ali_f814w_2026-11-26/REPORT.md),
[machine-readable status](data/ali_f814w_2026-11-26/result_and_validation_status.json),
and [reproduction commands](docs/ALI_F814W_REPRODUCE_20261128.md).

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


## Independent UNCOVER audit (corrected 2026-10-05)

UNCOVER/A2744 is an independent validation target using Li et al. (2026).

- Observed: 5 T dwarfs over 53.4 arcmin^2 (0.094 arcmin^-2).
- Li et al. define the density depth using the faintest object at
  **F277W=29.24 AB**. F115W=28.03 is an auxiliary magnitude for that object.
- The repository does not yet contain atmosphere-based `d25_F277W_pc`
  values, so the **primary matched-depth validation is pending**.
- The previous 2.88-object result is retained only as an F115W proxy and must
  not be labelled a successful UNCOVER validation.
- Using Li et al.'s own late-T d_max values gives a separate
  Galactic-density/geometry cross-check of about 1.10 predicted versus
  4 observed objects in 450-1050 K.

See `docs/uncover_validation_2026-10-05.md`,
`data/uncover_validation_2026-10-05.csv`, and run
`python -m src.validate_uncover`.

Next requirement: generate F277W synthetic reach from the atmosphere grid,
then forward-model the Li et al. colour/SED selection before making a
selection-matched validation claim.

