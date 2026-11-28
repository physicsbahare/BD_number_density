# Limiting-distance review for Ali's depth-limited estimand

Reviewed 2026-11-28. This is a definition of the population being counted, not a
new Galactic model or catalog-extraction simulation.

| Reference | Relevant definition | Implication for this calculation |
| --- | --- | --- |
| [Ryan & Reid (2016), AJ 151, 92, equations 8–10 and §4](https://arxiv.org/pdf/1510.05019) | Apparent-magnitude counts from the luminosity function and line-of-sight density. The formal distance integral extends to infinity while absolute luminosity constrains the distance contributing to each magnitude; cumulative counts can be restricted to a magnitude limit. They distinguish ideal population counts from completeness. | Supports a magnitude-defined population, not an arbitrary spatial cap or a universal five-sigma rule. Their disk parameters are not adopted. |
| [Aganze et al. (2022), ApJ 934, 73, §§2.2, 3.2–3.3, equations 5 and 7](https://escholarship.org/content/qt76c609m6/qt76c609m6.pdf) | Density-weighted effective volumes between luminosity-dependent distances determined from bright/faint field limits. Their simulation extends initially from 0.5 dmin to 2 dmax to address scattering and multiplicity. | Supports field- and luminosity-dependent boundaries. Their scatter extension shows that nominal dmax is not exact selection support. No new luminosity scatter or binaries are added here. |
| [Chen et al. (2025), PASA 42, e042, §5](https://arxiv.org/html/2503.18510v1) | A five-sigma flux limit and their color requirement imply an F444W magnitude limit, converted with a 1300 K model into a 4670 pc search distance. | Direct precedent for a sensitivity-to-model-distance conversion in COSMOS-Web; their color criterion, distance and Sonora-Bobcat model are not imported. |
| [Innala et al. (2026), FC-ENZO, §§III.1 and III.4](https://arxiv.org/html/2607.02354v1) | User-specified apparent-magnitude bins in a reference band define finite line-of-sight ranges; bin width and intrinsic luminosity scatter expand their bounds. | Supports separating the intrinsic magnitude-domain population count from the later selection fraction. Their parameters and extra scatter are not imported. |

The [public FC-ENZO implementation](https://github.com/onnalininn/Evaluator-of-Nearby-dwarf-stars-in-high-Z-photometrically-selected-Objects-ENZO/blob/3399b7598ccb56feca05f1f6b82e4051e1e1c6b0/FULL_synthesis_def.py)
was inspected as well: `number_density`, lines 1849–1856 in the inspected version,
uses the bright bin-edge flux for dmin and the faint edge for dmax, with inverse
square-root scaling. This resolves the reversed edge wording in the paper's
rendered §III.4 text. No upstream code is substituted into the preserved model.

## Adopted definition

At position x, require true background-limited F444W aperture signal

\[
F_{444,\mathrm{ap}}(T,r,x)\geq5\sigma_{444,\mathrm{emp}}(x),\qquad
r_{\max}(T,x)=d_{25,444}(T)
\sqrt{\frac{F_{25}\,\mathrm{EE}_{444}(x)}{5\sigma_{444,\mathrm{emp}}(x)}}.
\]

Here F25=10^((31.4−25)/2.5) nJy, and the existing d25, PSF response and correlated
noise calibration are retained. The equivalent total-flux depth is
`31.4 - 2.5*log10(5*sigma_emp/EE)` AB. No background offset is added to intrinsic
source luminosity. The calibrated net background offset belongs in the measured
native-SNR selection probability, separately from the sensitivity definition.

Within this population, integrate the unchanged density components times the
existing signed native-SNR selection probabilities, including the empirical
F814W CDF and source photon noise. This is the expectation of photometric
passes *among sources intrinsically above the defined depth*. It excludes
upward-scattered sources intrinsically fainter than that definition. It does
not assign a hard support boundary to the catalog, claim complete extraction,
or remove the earlier unrestricted-integral divergence by modifying the halo.

Five empirical sigma is a conventional background-limited sensitivity choice,
supported by the magnitude-limited literature and the calibrated survey noise;
it is not uniquely compelled by those papers. Native catalog SNR has a different
denominator, so native SNR=5 is not the same intrinsic depth. Boundary variants
4 and 6 empirical sigma and the native mean-SNR=5 radius quantify definition
sensitivity. They are alternative population definitions, not confidence limits.

The pre-existing 512-cell native-SNR calculation on this exact five-empirical-
sigma domain is reused as the central result. It was diagnostic under the
previous all-distance request; the user's follow-up explicitly authorizes its
finite population definition. No tuning to 23.843714, 26.500648 or 25.8331 occurs.
