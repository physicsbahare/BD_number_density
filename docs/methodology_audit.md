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


## Independent-field validation: UNCOVER/A2744 (2026-10-05)

A first independent sightline check has been completed using Li et al. (2026).

Using the UNCOVER effective area (53.4 arcmin^2), F115W=28.03 depth and the
repository's existing atmosphere-based F115W reach, the baseline Galactic
model predicts 2.88 objects across 450-1500 K versus 5 observed T dwarfs.
The model expectation is inside the exact central 68% Poisson interval for
five counts (2.849-8.365).

The colder 450-1050 K subset is more constraining: 4 objects are observed,
while the current F115W-depth model predicts 0.92; using Li et al.'s quoted
bin-specific d_max values gives 1.10. The latter sits at the lower edge of
the exact 95% Poisson interval for four counts. This is retained as a
temperature-distribution tension to investigate, not silently tuned away.

Before publication:
- forward-model the UNCOVER colour/SED selection;
- propagate local-density and atmosphere/distance uncertainties;
- test scale-height uncertainties and Teff-dependent thick-disk populations;
- retain the paper-text thick/thin=0.02 baseline, with 0.12 only as a
  documented sensitivity test.

Details: `docs/uncover_validation_2026-10-05.md`.
