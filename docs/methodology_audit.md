# Methodology audit

## What FC-ENZO actually provides

FC-ENZO can report the expected number of UCDs in a field without applying a high-z galaxy selection. Its number-density step integrates

N(SpT) = Omega * integral[rho(r,SpT) r^2 dr]

between magnitude-dependent distance limits.

The published method uses thin disk, thick disk and halo components, local densities tied to Kirkpatrick et al. (2021), and spectral-type-dependent scale heights from Aganze et al. (2022).

## Important distinction

FC-ENZO's final "contaminant" count is NOT the quantity needed first for our project. We first need the intrinsic expected UCD/BD population in COSMOS-Web. Then we apply each project's actual selection function.

## Two population products to report

### A. Conventional / warmer brown dwarfs
Report by spectral type and Teff, at minimum L, T and Y, plus useful sub-bins such as T0-T5. This can be benchmarked directly against Chen et al. (2025), who compared COSMOS-Web counts with Ryan & Reid (2016).

### B. Very cold UCDs
Report temperature bins extending into late-T/Y temperatures, with special attention below 700 K. The atmosphere grid and local luminosity-function completeness become increasingly important here. This population should not be inferred by simply extrapolating the Chen T0-T5 surface density.

## Critical implementation checks before trusting FC-ENZO numbers

The public FC-ENZO code requires an audit before reuse.

1. The paper states thick-disk/thin-disk normalization of 2%, while the public code inspected on 2026-09-26 contains `rho0_thick_factor = 0.12`. This must be resolved against the archived release/paper before reproducing published results.
2. The public `number_density` implementation appears to contain potentially ambiguous local-density unit handling. A fallback density array is already divided by 1000, while a later line divides the interpolated value by 1000 again. We must inspect `1Mj_fit.dat` and reproduce a published FC-ENZO field before changing anything.
3. The paper and public code use slightly different scale-height summaries in places. We will keep provenance for every parameter.
4. COSMOS-Web depth is non-uniform. A single-area/single-depth approximation is only a validation run, not the final result.
5. Selection completeness, morphology, blending, masking, and catalog detection are separate from the intrinsic Galactic count model.

## Validation ladder

1. Reproduce one published FC-ENZO test field.
2. Reproduce Chen et al. (2025) T0-T5 expected count at F115W 5-sigma = 27.45 over 0.243 deg^2: 13.1 +/- 7.9 expected from the Ryan & Reid model.
3. Run the same calculation over the nominal full COSMOS-Web area.
4. Replace uniform depth by COSMOS-Web exposure/depth zones.
5. Produce counts by SpT, Teff, magnitude, Galactic component and survey zone.
6. Apply project-specific selections separately.
7. Propagate parameter uncertainty plus Poisson uncertainty and model-systematic envelopes.

## Outputs required for papers

- expected N and surface density per arcmin^2
- counts versus F444W and F115W
- counts by L/T/Y and Teff bins
- thin/thick/halo contributions
- nominal full footprint and actual usable footprint
- uncertainty budget
- comparison with Chen et al. observed candidates and Ryan & Reid prediction
- separate selection-matched predictions for each project


## Independent-field audit: UNCOVER/A2744 (corrected 2026-10-05)

Li et al. (2026) define the density depth using the faintest sample object at
F277W=29.24 AB. F115W=28.03 is also quoted for that object, but it is not the
primary depth definition for their density calculation.

The repository currently lacks atmosphere-based F277W d25 values. Therefore
the primary matched-depth UNCOVER validation is **pending**. The earlier
F115W=28.03 calculation (2.88 objects over 450-1500 K) is retained only as
an auxiliary proxy and must not be cited as a successful validation.

A separate cross-check using Li et al.'s own Table-4 d_max values remains
valid as a Galactic-density/geometry test: the baseline model gives about
1.10 objects across 450-1050 K versus 4 observed. This lies at the lower edge
of the exact 95% Poisson interval and is retained as a temperature-
distribution tension to investigate.

Before publication:
- generate atmosphere-based F277W reach and add `d25_F277W_pc`;
- rerun the primary validation at F277W=29.24;
- forward-model the UNCOVER colour/SED selection;
- propagate local-density and atmosphere/distance uncertainties;
- test scale-height uncertainties and Teff-dependent thick-disk populations;
- retain the paper-text thick/thin=0.02 baseline, with 0.12 only as a
  documented sensitivity test.

Details: `docs/uncover_validation_2026-10-05.md`.
