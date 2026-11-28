# Executed validation follow-up to 5c8cd56

Read the [scientific resolution and reproduction commands](../../docs/ALI_VALIDATION_RESOLUTION_20261128.md).

`result_status.json` is authoritative: no finite final count, surface density or
combined uncertainty is accepted. The all-distance forced-aperture integral
diverges under the fixed untruncated halo. Geometric five-band area is
0.532430830667 deg². Exact historical mask attribution remains unidentified.

- `footprint_ledger.json` and `detection_footprint_validation.json`: measurement coverage versus an intentionally invalid detection-intensity mask; recovered per-tile header omissions are explicit.
- `native_zero_source_selection.json`: existing held-out real backgrounds, with no new injection; reconstructed estimators, not catalog source detections.
- `domain_status.json`: asymptotic convergence check and modeled zero-source floor.
- `native_radial_shell_counts.csv`: executed component/bin radial integrals. Absolute-radius and local-r50-scaled rows are alternative diagnostic domains and must not be added together.
- `cutoff_diagnostics.csv` and `cutoff_bin_diagnostics.csv`: cumulative finite-cutoff diagnostics, never adopted survey counts. Formal density errors omit lower-limit-bin uncertainties and model systematics.
- `nominal_boundary_scatter.csv`: scattering beyond the local 50% F444W boundary, with its zero-source floor separately integrated.
- `relative_floor_convergence.json`: 64/128-node floor quadrature comparison.
- Three SVG figures visualize these diagnostic tables. The floor-subtracted excess is a different estimand and must not be quoted as the selected UCD population.
- `scientific_review.md`: independent review of the executed diagnostics.
- `validation_tests.txt`: targeted probability, coverage and tail tests.

All completed F814W recovery/calibration products remain in the preceding
`ali_f814w_2026-11-26` release. No raw FITS, new atmosphere grids or reconstructed
depth caches are included here.
