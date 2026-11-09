"""Small offline tests: no COSMOS-Web science mosaics are implied."""
import numpy as np
import pytest

from injection_recovery.core import (
    ali_photo_selection, dave_photo_selection, inject_psf_njy,
)


def dave_baseline():
    f = dict(F115W=-30.0, F150W=0.0, F277W=-1.0, F444W=700.0)
    e = {b: 10.0 for b in f}
    return f, e


def test_dave_signed_blue_noise_and_nonpositive_f277():
    flux, error = dave_baseline()
    assert dave_photo_selection(flux, error, 25.0)
    # S/N<-2 must still pass the *signed*, not absolute-value, rule.
    flux["F115W"] = -100.0
    assert dave_photo_selection(flux, error, 25.0)
    # Positive F277W must NOT be replaced with 2*sigma.
    flux["F277W"] = 25.0
    assert not dave_photo_selection(flux, error, 25.0)
    flux["F277W"] = 10.0
    assert dave_photo_selection(flux, error, 25.0)


def test_dave_photometry_missing_is_not_dropout():
    flux, error = dave_baseline()
    del error["F115W"]
    assert not dave_photo_selection(flux, error, 25.0)


def test_dave_limits_are_strict():
    flux, error = dave_baseline()
    assert not dave_photo_selection(flux, error, 26.2)
    flux["F444W"] = 50.0
    assert not dave_photo_selection(flux, error, 25.0)


def test_ali_needs_hst_and_catalogue_quality():
    flux, error = dave_baseline()
    flux["F814W"], error["F814W"] = -5.0, 10.0
    assert ali_photo_selection(flux, error, f444w_auto=100.0, warn_flag=0)
    assert not ali_photo_selection(flux, error, f444w_auto=0.0, warn_flag=0)
    assert not ali_photo_selection(flux, error, f444w_auto=100.0, warn_flag=1)
    flux["F814W"] = 25.0
    assert not ali_photo_selection(flux, error, f444w_auto=100.0, warn_flag=0)
    del flux["F814W"]
    assert not ali_photo_selection(flux, error, f444w_auto=100.0, warn_flag=0)


def test_injection_keeps_originals_and_adds_source_variance():
    image = np.zeros((31, 31))
    variance = np.full((31, 31), 0.1)
    psf = np.ones((7, 7))
    out, vout = inject_psf_njy(
        image, variance, psf, x=15.2, y=14.8,
        flux_njy=100.0, electrons_per_njy=10.0,
        rng=np.random.default_rng(7),
    )
    assert np.all(image == 0)
    assert np.all(variance == 0.1)
    assert np.isclose((vout - variance).sum(), 100.0 / 10.0)
    assert np.any(out != image)
    assert np.isfinite(out).all()


def test_refuse_edge_or_invalid_calibration():
    image, variance, psf = np.zeros((21, 21)), np.ones((21, 21)), np.ones((5, 5))
    with pytest.raises(ValueError):
        inject_psf_njy(image, variance, psf, 0, 0, 100, 10, np.random.default_rng())
    with pytest.raises(ValueError):
        inject_psf_njy(image, variance, psf, 10, 10, 100, 0, np.random.default_rng())


def test_cutout_pilot_writes_fits_and_requires_ali_hst(tmp_path):
    from astropy.io import fits
    from injection_recovery.run_cutout_pilot import run

    entries = {}
    for band in ("F115W", "F150W", "F277W", "F444W"):
        sci = tmp_path / f"{band}_sci.fits"
        var = tmp_path / f"{band}_var.fits"
        psf = tmp_path / f"{band}_psf.fits"
        fits.writeto(sci, np.zeros((31, 31)), overwrite=False)
        fits.writeto(var, np.ones((31, 31)), overwrite=False)
        fits.writeto(psf, np.ones((5, 5)), overwrite=False)
        entries[band] = {
            "science_path": str(sci),
            "variance_path": str(var),
            "psf_path": str(psf),
            "x": 15.0,
            "y": 15.0,
            "flux_njy": 100.0,
            "electrons_per_njy": 10.0,
        }
    manifest = dict(project="dave", injection_id="fixture01", bands=entries,
                    wcs_alignment_verified=True, image_units_verified=True)
    result = run(manifest, tmp_path / "injected", seed=11)
    assert result["status"] == "IMAGE_INJECTION_ONLY_NOT_RECOVERED_NOT_COMPLETENESS"
    assert len(result["bands"]) == 4
    with fits.open(tmp_path / "injected" / "fixture01_F444W.fits") as hdus:
        assert hdus["SCI"].shape == (31, 31)
        assert np.isclose((hdus["VAR"].data - 1.0).sum(), 10.0, atol=1e-4)
    with pytest.raises(FileExistsError):
        run(manifest, tmp_path / "injected")
    manifest["project"] = "ali"
    with pytest.raises(ValueError, match="F814W"):
        run(manifest, tmp_path / "ali_injected")
