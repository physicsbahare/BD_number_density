# Reconstruct and execute the existing Ali count calculation

The user authorizes rebuilding missing intermediate depth maps from local
images. Historical counts and area are validation references, not fit targets.
No changes to density bins, Galactic components, normalizations, or atmosphere
grid are authorized. Implementation is in the existing Python repository.

1. Recover original solar-metallicity FC-ENZO Elf Owl five-band photometry,
   including its radius normalization. Verify median d25 against the existing
   thirteen-bin table, independently of number counts.
2. Stream/decompress one local tile at a time. Preserve science WCS, explicitly
   verify extracted uncertainty-array dimensions, and sample joint coverage
   on a common native grid with overlap removal. Save resumable tile caches.
3. Measure native-image blank apertures at radii 0.10 and 0.15 arcsec, masking
   detected sources and invalid pixels. Calibrate spatial noise empirically,
   verify the adopted 0.10 arcsec SNR against the actual catalog, and retain
   both conventions to diagnose the historical baseline.
4. Reuse the repository Galactic density functions and numerical quadrature.
   Apply signed Gaussian photometric probabilities; integrate the same
   population on four-band original and five-band common footprints.
5. Evaluate and validate the existing spatial PSF against local stars. Inject
   F814W flux into real backgrounds and recover using the verified aperture
   and noise convention. Test fluxes below, near and above 2 sigma; increase
   samples and retain empirical selection only if justified.
6. Check geometry and integration convergence, noise sensitivity, density
   uncertainties, matched-footprint monotonicity, and injection statistics.
   Save results, diagnostic figures, execution commands, and paper paragraph.

Recoverability finding: the original FC-ENZO solar grid contains 7200 models.
Its per-bin median sqrt(flux_zp25)*Radius_pc reproduces both existing d25
columns to their recorded precision. This establishes the early table's
synthetic flux normalization and bin representative statistic.

Pending scientific choices must be recorded explicitly; a failed benchmark
must be reported rather than repaired by fitting the old answer.
