# F277W hybrid flux-ratio cut: a discontinuous selection response (2026-11-01)

**Status:** source-verified mathematical diagnostic, **NOT** a measured effective area or proper-motion completeness.

The 2026-10-29 Dave manuscript Section 2.2, independently reproduced against all 3,280,023 v13 catalog rows (2997 -> 335 -> exact frozen 125), uses corrected `new_e_*` errors, signed blue-band `f/sigma<2`, and the F444W/F277W >30 rule. For **positive measured** F277W, the ratio uses measured `f277`; for **nonpositive** F277W it substitutes `2*sigma277`.

Define `T=f444/30`, `s=sigma277`, and measured `Y~Normal(mu277,s²)` at fixed positive F444W flux. The F277W-only pass region is:

- When `2s<T`: all `Y<2s` pass (negative and positive branches); `P=Phi(2-mu277/s)`.
- When `2s>=T`: **only** `0<Y<T` passes; `P=Phi((T-mu277)/s)-Phi(-mu277/s)`.

At `mu277=0`, the probability drops discontinuously from `Phi(2)=0.977249868` immediately below `s=T/2` to `Phi(2)-0.5=0.477249868` at/above the boundary. The 50-percentage-point jump comes from rejecting **all negative measured F277W fluxes** when `2s>=T`. This is a property of the *stated algorithm*, not an astrophysical effect. A realistic source population, varying F444W noise, correlated errors, background biases and the other selection stages will modify any survey-wide result.

At the **brighter-side limit** `m444->26.2^-` (the paper uses strict `m444<26.2`), `f444->120.22644 nJy`, `T->4.00755 nJy` and `s_crit->2.003774 nJy`. This is the same threshold that would require a 5sigma F277W depth of about 28.898 AB to *prove* a >30 flux ratio using a zero-flux 2sigma upper limit. It is **not** a universal depth cut: the actual paper has two distinct branches.

## New source-level evidence from already committed, sky-verified 125 v13 rows

The source-specific critical sigma is `f444/60`, **not** a constant 2.004 nJy, because candidate F444W fluxes differ.

| subset | sources | negative F277W measurement can pass? | negative branch forbidden |
| --- | ---: | ---: | ---: |
| all frozen candidates | 125 | **57** | **68** |
| >=365d F356W geometric epoch | 62 | **20** | **42** |
| without that epoch | 63 | **37** | **26** |

Of the **26 actually nonpositive F277W measurements**, all 26 are in the branch-permitted regime, as required by the already-reproduced selection. **Four candidates** have `sigma277/(f444/60)` within 10% of unity: **311339, 592641, 1687583, 1651050**. Candidate **1651050** is just above the boundary (ratio **1.09695**); candidate **387132** is farther above (**2.70944**). Their actual F277W measurements are positive. This does **not** alter their observed astrometric displacement evidence or validate intrinsic PM for 1651050.

For illustration only, the **median** conditional zero-true-F277W Gaussian pass probability across the *already-selected source positions* is **0.38544** in the 62 red-epoch sources and **0.97725** in the other 63. These are **NOT completeness estimates**: the source positions and noise values are selection-conditioned and do not sample survey sky uniformly. The strong contrast is a warning that multiplying the 59.48% geometric epoch fraction by a uniform dropout-selection probability is unjustified.

## Consequences and next steps

1. **Clarify the scientific intention** of treating positive and negative F277W flux asymmetrically. Do not silently modify the published rule while reproducing the paper; any revised rule requires a new full-catalog selection and new 125-source accounting.
2. To estimate `A_eff(m444=26.2)`, evaluate the full joint four-band selection probability on a **spatially representative**, calibrated ERR/depth/PSF/mask grid, with correlated noise and aperture errors; do not use the 125 preselected source errors as a sky-area distribution.
3. Separate photometric selection area from >=365d F356W geometric opportunity and from actual second-epoch SNR/astrometric PM sensitivity. Existing v13 WHT>0 ~0.533619 deg², WHT∩red365 ~0.317401 deg² (59.48%) remain *geometric*.
4. Direct star-mask unique-sky union ∩ WHT ∩ epoch remains unmeasured; existing mask tarball was already processed and must **not** be requested or redownloaded. The sparse local script is ready, but the 51GB original WHT is on the user's disk.

Reproduce with `python -m pm_footprint.hybrid_f277w_selection_probability`; targeted offline tests: `python -m unittest -v pm_footprint.test_hybrid_f277w_selection_probability`. Source data: `results/v13_epoch_candidate_photometry_fixture_20261031.csv`; compact evidence: `results/hybrid_f277w_selection_probability_20261101.json`.
