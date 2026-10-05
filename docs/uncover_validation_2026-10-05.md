# Independent UNCOVER validation - 2026-10-05

## Goal

Use the same Galactic line-of-sight machinery developed for COSMOS-Web on an
independent JWST deep field: UNCOVER around Abell 2744.

The comparison source is Li et al. (2026), *Two late-T dwarfs at kiloparsec
distances revealed by JWST UNCOVER survey*, MNRAS 547, stag227,
doi:10.1093/mnras/stag227.

Li et al. report:
- effective common imaging area = 53.4 arcmin^2;
- five T dwarfs in that footprint;
- measured surface density = 0.094 arcmin^-2;
- limiting depth F277W = 29.24 AB, with F115W = 28.03 AB quoted from the
  faintest T dwarf;
- late-T temperature-bin counts of 2, 1, 0, 1 in 450-600, 600-750,
  750-900 and 900-1050 K, respectively;
- detection-cone d_max values of 0.7, 1.3, 2.0 and 2.6 kpc for those bins.

The fifth object is the previously known sdT1 candidate with a best-fit
temperature of about 1400 K, so the first-pass all-T comparison uses
450-1500 K.

## Method

The validation deliberately preserves the repository's baseline assumptions:

- Kirkpatrick et al. (2024) local 150-K temperature-bin densities;
- thin + thick disk + halo density integration;
- thin-disk scale height H = 187 pc;
- thick/thin local normalization = 0.02;
- halo/thin local normalization = 0.0025;
- the existing FC-ENZO/Elf-Owl-derived F115W d25 values.

For the primary check, each repository d25(F115W) value is scaled from
m_AB=25 to the Li et al. F115W depth of 28.03 via the distance modulus, then
integrated along the Abell 2744 line of sight over 53.4 arcmin^2.

A second, narrower check uses Li et al.'s own d_max values for the four
450-1050 K bins. This isolates the Galactic-density part of our calculation
from atmosphere-based detection reach.

This is a depth-only intrinsic-population validation. It does **not** yet
forward-model the Li et al. colour cuts, SED-fit acceptance, source
completeness, blending or morphology/point-source selection.

## Reproducible results

| Teff (K) | Observed | dmax from our F115W reach (pc) | Model N | Li dmax (pc) | Model N at Li dmax |
|---|---:|---:|---:|---:|---:|
| 450-600 | 2 | 484 | 0.132 | 700 | 0.214 |
| 600-750 | 1 | 1010 | 0.227 | 1300 | 0.285 |
| 750-900 | 0 | 1761 | 0.247 | 2000 | 0.285 |
| 900-1050 | 1 | 2601 | 0.314 | 2600 | 0.314 |
| 1050-1200 | 0 | 3574 | 0.426 | - | - |
| 1200-1350 | 0 | 4491 | 0.909 | - | - |
| 1350-1500 | 1 | 5673 | 0.625 | - | - |

### All T-like objects, 450-1500 K

The model predicts

N_model = 2.880

over 53.4 arcmin^2, corresponding to

Sigma_model = 0.0539 arcmin^-2.

Li et al. observe 5 objects,

Sigma_obs = 5 / 53.4 = 0.0936 arcmin^-2,

consistent with their rounded published value of 0.094 arcmin^-2.

For n=5, the exact central 68% Poisson interval is 2.849-8.365 counts.
The model expectation of 2.880 lies just inside this interval. Therefore the
independent field does **not** show a significant disagreement with the
current depth-only Galactic model. The observed/model ratio is 1.74, which
is modest given Poisson noise, atmosphere-to-distance systematics, and the
fact that the Li selection has not yet been forward modeled.

### Late-T check, 450-1050 K

Observed late-T count: 4.

Using our F115W-based detection reaches:

N_model = 0.920.

Using Li et al.'s quoted d_max values:

N_model = 1.098.

This temperature distribution is more late-T-rich than the baseline model.
For n=4, the exact 95% Poisson interval is 1.090-10.242 counts. The
Li-dmax prediction of 1.098 sits at the lower edge of that interval, while
the atmosphere-based 0.920 value is slightly below it.

This should be treated as a useful tension, not as a failed validation.
Li et al. themselves note that their T8-T8.9 density in the thick disk is
comparable to the local value and discuss cooling of an old thick-disk
population toward later spectral types. Our current model partitions a
single local Teff-bin density into thin/thick/halo components using fixed
normalizations. The UNCOVER result therefore motivates testing whether the
Teff distribution of the thick disk should be modeled separately.

## Interpretation

The first-pass independent-field test is encouraging:

1. The total 450-1500 K population prediction is consistent with the five
   observed T dwarfs within exact 68% Poisson uncertainty.
2. The temperature-bin distribution is not reproduced equally well: the
   observed sample contains more very late/cold T dwarfs than the baseline
   fixed-component model predicts.
3. This makes UNCOVER a stronger validation target than a single scalar
   count. It can constrain the component/temperature treatment before the
   final COSMOS-Web prediction is published.

## Next validation steps

- Forward-model the Li et al. photometric colour selection and SED-fit
  acceptance rather than comparing only depth-limited intrinsic counts.
- Add uncertainty propagation from Kirkpatrick density errors and
  atmosphere-based d_max.
- Test the Aganze scale-height uncertainty range.
- Test a Teff-dependent thick-disk population instead of assigning a fixed
  thick/thin ratio to every Teff bin.
- Keep the public FC-ENZO 0.12 thick-disk normalization only as a sensitivity
  test; do not silently replace the paper-text 0.02 baseline.
- Once these checks are stable, rerun COSMOS-Web and report UNCOVER as an
  independent external validation in the methods/results section.

Reproduce with:

```bash
python -m src.validate_uncover
```

The numerical snapshot is stored in
`data/uncover_validation_2026-10-05.csv`.
