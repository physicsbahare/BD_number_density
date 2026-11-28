# Review of the authorized depth-limited Ali estimand, 2026-11-28

Read-only review of `src/ali_depth_limited.py`, the depth-branch changes in `src/ali_spatial_counts.py`, `docs/ALI_DEPTH_LIMITED_20261128.md`, and `docs/ALI_DEPTH_LIMITED_LITERATURE_20261128.md`. No new computations or injections. Threshold/r50/CDF uncertainty products were still pending when this draft review was requested.

## Scientific assessment

The newly authorized estimand is scientifically defensible: expected signed-photometry passes among Galactic model UCDs whose true F444W aperture signal is at least five times the local empirical background noise. This is a population restriction, not an assertion about the support of the full catalog selection. The existing local intrinsic-depth boundary gives a finite integral without changing the Galactic model. Earlier all-distance divergence remains correct for the unrestricted estimand and does not invalidate this newly authorized finite definition.

The draft clearly states that intrinsic empirical five-sigma depth and measured native catalog SNR>5 are different quantities. Net background offsets correctly belong in measured-selection probabilities rather than stellar luminosities. Native probabilities are retained inside the depth domain. The formula's flux and distance normalization matches the reviewed core implementation; no fixed 20/50 kpc limit is applied by the depth branch.

The adopted 12.360844 prediction can therefore be reported as a **depth-limited photometric foreground approximation**, conditional on the preserved population model and calibration. It cannot be called an exact catalog-entry prediction, completeness-corrected recovery count, or reconstruction of the archived historical result. The draft supplies these qualifications, including the separate AUTO-positivity limitation. Missing historical mask provenance no longer blocks reporting this independently defined, explicitly authorized quantity; it still blocks historical reproduction claims.

## Driver review

- `ali_depth_limited.py:45-60` keeps measured noise, native thresholds and offsets separate from the altered population boundary. Its four/six-sigma and native-r50 options affect the integration domain as intended.
- The native-r50 denominator includes the net mean offset, `5*ref_F444-mu_F444`; with the default five-sigma parameter its boundary scale yields the stated mean-native-SNR=5 radius.
- The source-photon-noise coefficient and Hermite convolution are the previously reviewed ACS treatment. Source photon noise is explicitly included for ACS; the draft should not imply newly implemented photon-noise propagation in all four NIRCam bands.
- Spatial grouping, donor-ACS provenance and area weighting follow the previously reviewed implementation. New 1024-leaf comparison addresses the additional compression approximation; adaptive radial error is not a bound on compression or calibration error.
- The earlier correlation mismatch is resolved in the default scenario: native r010 correlation uses 0.1671. The retained r015 value 0.224419 is explicitly a legacy stress test, not a measured native covariance. Neither single-pair Gaussian sensitivity constitutes a full joint-band empirical model.
- Four/six empirical sigma and native r50 define different populations. The draft correctly calls them definition sensitivities, not confidence limits.

## Literature and claims

- [Ryan & Reid, equations 9-10](https://arxiv.org/pdf/1510.05019) use magnitude-count integrals with formally unbounded distance integration. They support magnitude-defined counts but do not prescribe a common finite distance cap. The literature note accurately distinguishes this from the adopted implementation.
- [Aganze et al., section 3.2](https://arxiv.org/pdf/2204.07621) define density-weighted effective volumes with field/luminosity-dependent limiting distances. This supports the general construction; their treatment of scatter and selection must not be presented as identical to the restricted representative-bin calculation. The draft explicitly retains that distinction.
- [Chen et al., section 5](https://arxiv.org/html/2503.18510v1) derive their model search distance from sensitivity and a color requirement. This is a relevant precedent for a model-dependent survey-depth conversion, not justification for importing their distance or atmosphere model.
- [FC-ENZO, sections III.1-III.4](https://arxiv.org/html/2607.02354v1) specifies reference-band magnitude bins and separates population counts from subsequent selection fractions. Its expanded luminosity/bin bounds are not identical to the present single representative flux per temperature bin. The draft correctly says the papers support this *type* of estimand, not a uniquely mandated empirical-five-sigma boundary.
- The upstream GitHub line-number citation points to mutable `main`. Prefer a commit-pinned permalink or record the inspected commit/hash for the asserted bright/faint-edge code correction. I could not independently refetch that GitHub page in this review; this is a provenance improvement, not evidence that the described inverse-square-root scaling is wrong.

## Draft numerical and interpretation checks

The central count, components, formal measured-density sigma and four-band matched prediction agree with the previously reviewed cached native-model results. The draft distinguishes reconstructed pixel-center four-band area from aperture-valid area. Its ACS-only comparison uses the identical common footprint and intrinsic F444W domain, rather than comparing conventions or footprints.

The old 23.843714 result changes NIRCam aperture/noise/SNR conventions, so the difference from 12.360844 is not an isolated ACS suppression. The draft states this. The lower-limit cold bin and unavailable model/calibration uncertainties prevent treating the formal measured-density error as a complete symmetric population interval.

Before marking the new result package complete, collect the pending threshold/r50/CDF products and ensure all referenced files exist and use native maps. Retain the conditional finite-CDF interpretation and separately report any step-versus-interpolated-CDF central-value difference. No extra expensive observations or modeling is needed for this reporting check.

## Conclusion

No new critical scientific or implementation objection was found to adopting the explicitly depth-limited approximation. The necessary distinctions are present: finite population definition versus full selection, empirical sensitivity versus native SNR, model expectation versus observed count, and definition sensitivity versus statistical uncertainty. Final numerical packaging remains contingent on completion and consistency review of the pending products; historical reproduction and exact catalog completeness remain unclaimed.

## Final completed-package review

Reviewed the completed `data/ali_depth_limited_2026-11-28/result.json`, sensitivity summary, per-bin selected-volume table, CDF uncertainty artifact, provenance/README, and updated scientific report. This was an artifact consistency review; no numerical science runs or injections were executed.

- The adopted count 12.360844100804112, common area 0.5324308306666868 deg2, surface density 0.006448852335390787 arcmin^-2 and components match the reused native-model central result. Matched four-band 14.861435786045444 minus five-band count gives the reported 2.500591685241332 reduction (16.8260%).
- Formal measured-density sigma 0.8414434272461976, conditional CDF sigma 0.025518640987744215 and their partial quadrature combination 0.8418302930482416 are consistently labeled. The lower-limit bins total 5.036238646467148, including a 2100-2250 K contribution 2.8445e-11; their unquantified uncertainties are not silently included as measured errors.
- CDF step diagnostic 12.348251077989005 differs from the adopted interpolated-CDF result by 0.012593022815106991. The report correctly separates this approximately 0.102% implementation difference from calibration error and from the approximate conditional finite-CDF SE. The artifact records 256 radial nodes, 15 source-noise nodes and 512 spatial groups per tile.
- Sensitivity summary agrees with the report: four/six-sigma domains 16.0859767748/9.8967367492; ACS noise scales 0.95/1.05 give 12.5327292783/12.2051604863; signed ACS thresholds 1.5/2.5 give 11.3592208360/13.1254184415; alternative native-r50 domain gives 28.7226824904. These remain alternative definitions or calibration stresses, not confidence limits.
- Native covariance sensitivity is correctly the r010 value 0.1671, producing 12.3602662282. The legacy 0.224419 stress is separately labeled. Spatial 1024-group count 12.3692314893 differs by 0.0083873885, distinct from the central adaptive radial error estimate.
- Per-bin selected effective volumes are explicitly N divided by total adopted local-bin density, incorporating population fractions and selection. Package README clarifies that boundary-file `r50_pc` is a spatial median, distinct from the mean-native-SNR=5 boundary scenario.
- Acceptance is explicitly limited to the newly authorized depth-limited population approximation. The documentation retains historical-mask uncertainty, all-distance divergence, AUTO-positivity and catalog-completeness exclusions. It does not promote any cutoff or conditional calculation into a full native catalog prediction.

**Completion verdict:** pending products from the draft review are now present and internally consistent with the adopted result and report. No material correction was identified before committing this completed package. The earlier scientific/model limitations remain appropriately explicit; acceptance applies only to the defined depth-limited photometric foreground approximation.
