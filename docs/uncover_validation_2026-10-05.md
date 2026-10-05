# UNCOVER/A2744 validation audit - corrected 2026-10-05

## Status

**Primary matched-depth validation is pending.**

Li et al. (2026), *Two late-T dwarfs at kiloparsec distances revealed by JWST UNCOVER survey* (MNRAS 547, stag227), define the depth for their surface- and space-density calculation using the faintest T dwarf in the sample at **F277W = 29.24 AB**. They also quote **F115W = 28.03 AB** for that same object.

The first repository implementation incorrectly treated F115W=28.03 as the main validation depth. That value can be kept as an auxiliary reference, but it must not replace the paper's primary F277W depth definition.

## Li et al. quantities used here

- effective common imaging area: **53.4 arcmin^2**;
- observed T dwarfs in that area: **5**;
- measured surface density: **0.094 arcmin^-2**;
- primary depth definition: **F277W = 29.24 AB**;
- auxiliary/reference magnitude for the same faintest object: **F115W = 28.03 AB**;
- Table-4 late-T counts in 450-600, 600-750, 750-900 and 900-1050 K: **2, 1, 0, 1**;
- Table-4 d_max values for those bins: **0.7, 1.3, 2.0 and 2.6 kpc**.

## What the repository can currently test

The COSMOS-Web temperature-bin table currently contains atmosphere-based d25 values for F115W and F444W, but **not F277W**. Therefore the repository does not yet contain the quantity needed for a like-for-like primary UNCOVER depth calculation:

`d25_F277W_pc(Teff)`

The validation script now requires `d25_F277W_pc` to be present, finite, and strictly positive for **all seven required 150-K bins from 450 to 1500 K**. Missing bins, duplicate required bins, NaN values, zero/negative reaches, or a missing F277W column all keep the primary validation in the pending state.

## Auxiliary F115W proxy

For reproducibility, the previous F115W calculation is retained but relabelled as an **auxiliary proxy**. Using F115W=28.03 and the existing F115W d25 values gives:

- N_proxy(450-1500 K) = **2.88** over 53.4 arcmin^2;
- Sigma_proxy = **0.0539 arcmin^-2**;
- observed Sigma = 5/53.4 = **0.0936 arcmin^-2**.

Although 2.88 happens to lie inside the exact central 68% Poisson interval for five observed counts, this must **not** be described as a successful primary UNCOVER validation because the detection band/reach is not matched to the paper's F277W depth definition.

## Table-4 d_max cross-check

A separate test remains useful: take Li et al.'s own late-T d_max values and test only the Galactic thin+thick+halo density integration. With the current baseline assumptions the four 450-1050 K bins give **N_model = 1.10** versus **4 observed** late-T dwarfs.

For four observed counts, the exact central 95% Poisson interval is about 1.09-10.24 counts, so the model lies at the lower edge. This is a useful temperature-distribution tension, not a replacement for the missing primary F277W calculation.

## Required next step

Generate synthetic **F277W** photometry from the same atmosphere grid used for the COSMOS-Web reach calculation and add a `d25_F277W_pc` column to the temperature-bin input table. Only then should the code report a primary F277W=29.24 UNCOVER prediction.

After that, forward-model the Li et al. photometric colour selection, SED-fit acceptance, completeness and source-selection effects.

## Publication rule

Until `d25_F277W_pc` exists and the primary run has been completed, do not state that the COSMOS-Web model reproduces or validates against the five UNCOVER T dwarfs. The defensible statement is that UNCOVER is an independent validation target, with an F115W proxy and a Table-4 d_max geometry check already completed while the matched F277W calculation remains pending.

Run the current audit with:

```bash
python -m src.validate_uncover
```

The generated snapshot is stored in `data/uncover_validation_2026-10-05.csv`.