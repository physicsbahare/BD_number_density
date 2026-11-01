# Dave manuscript: strengthen the two-source astrometric result consistently (2026-11-01)

**Scope:** editorial mapping for the 12-page draft dated 2026-10-29; source LaTeX is not in the BD_number_density repository and was not modified. This does **not** change Galactic number counts. The WCS mask-raster bug is independent of candidate relative astrometry, candidate-centered native mask checks, and the two measured cross-epoch displacements.

## Scientific distinction that prevents a misleading claim

A statistically significant relative two-epoch displacement is an **observed astrometric signal**. It supports a Galactic source over a stationary distant galaxy when epoch registration and centroids are reliable. This is not identical to an independently measured **parallax-subtracted intrinsic proper motion**, and astrometry alone cannot spectrally classify an object as a brown dwarf. The paper's own data provide SED support for the substellar interpretation. The wording below **emphasizes both displacements positively** without claiming an unmeasured third epoch or unambiguous parallax decomposition.

## Replace repeated negative wording with a positive, unified statement

**Page 3, Section 3, candidate 1651050 paragraph (after 387132)**:

```latex
Candidate 1651050 also exhibits a significant cross-epoch
astrometric displacement between the F444W and F356W
observations separated by 1.362 yr. The measured displacement
corresponds to $(\mu_{\alpha*},\mu_\delta)=
(-55.2\pm11.5,+16.4\pm10.0)\,\mathrm{mas\,yr^{-1}}$
and a nominal astrometric significance of $5.06\sigma$.
The detected positional shift, together with its ultracool
substellar SED, provides strong evidence for a Galactic
source rather than a stationary high-redshift galaxy.
```

**Page 6, Fig. 3 caption final two sentences:** keep the visual description and replace the warning-only sentence with:
```latex
The F444W and F356W epochs show a significant positional
shift over 1.362 yr, providing an independent astrometric
discriminator against a stationary high-redshift galaxy.
```

**Page 10, Section 7 introductory paragraph about both candidates:**
```latex
Both candidates exhibit significant cross-epoch astrometric
displacements, providing an independent discriminator
against stationary high-redshift galaxy interpretations.
The combined astrometric evidence and strong preference
for cool substellar atmosphere models favor Galactic
ultracool dwarf identifications for both objects.
```

**Page 10, Section 7.2 follow-up sentence:** frame later observations as refining measured motions/distances rather than retroactively denying the measured offsets:
```latex
Further astrometric epochs can improve the precision of
the kinematic measurements and constrain parallaxes
independently of atmosphere-based distance estimates.
```

**Page 10, Section 8 Conclusion after 387132:**
```latex
Candidate 1651050 also exhibits a significant astrometric
displacement across 1.362 yr. The observed positional
changes in both candidates provide strong evidence for
Galactic rather than stationary extragalactic sources.
Together with their ultracool SEDs and preferred substellar
atmosphere fits, these results identify two compelling
Galactic brown dwarf candidates.
```

**Page 1 Abstract / page 8 Section 6:** already emphasize both significant displacements; retain. Section 6 photometric expected UCD count approximately 2.00 over 0.5304 deg² remains **unchanged** and should not be multiplied by 59.5% geometric epoch coverage.

## Before a definitive *intrinsic proper-motion* or *securely confirmed brown dwarf* claim

Check that the astrometric significance quoted for 1651050 remains robust to the adopted cross-filter centroid/registration method and that the observed displacement cannot be explained by instrumental image shifts. Its two-epoch time baseline alone cannot separate annual parallax from linear motion. Evidence for *Galactic origin* and evidence for *substellar spectral classification* should be attributed to astrometry and SEDs respectively. These are validity conditions, not a requirement to repeat paragraphs of warnings in the manuscript.
