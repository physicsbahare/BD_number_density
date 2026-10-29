"""Offline WHT test: a tiny fake FITS replaces multi-GB public downloads."""
import gzip
import math
import shutil

import numpy as np
import pytest

pytest.importorskip("astropy")
pytest.importorskip("shapely")

from astropy.io import fits
from astropy.wcs import WCS
from pm_footprint import wht_depth_mask as wht


def test_official_tile_url():
    assert wht.url_for("A1").endswith(
        "/mosaic_nircam_f444w_COSMOS-Web_60mas_A1_v1.0_wht.fits.gz")
    assert len(wht.TILES) == 20
    with pytest.raises(ValueError):
        wht.url_for("A11")


def test_wht_fractional_grid_no_network(tmp_path, monkeypatch):
    hdr = WCS(naxis=2)
    hdr.wcs.crpix = [50.0, 50.0]
    hdr.wcs.crval = [150.12, 2.20]
    hdr.wcs.cdelt = [-0.06 / 3600, 0.06 / 3600]
    hdr.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    values = np.full((100, 100), 4.0, dtype=np.float32)
    values[:50, :50] = 0.0
    raw = tmp_path / "fake.fits"
    fits.PrimaryHDU(values, header=hdr.to_header()).writeto(raw)
    gz = tmp_path / "fake.fits.gz"
    with open(raw, "rb") as fin, gzip.open(gz, "wb") as fout:
        shutil.copyfileobj(fin, fout)

    def copy_instead_of_download(url, target):
        shutil.copyfile(gz, target)

    monkeypatch.setattr(wht, "download", copy_instead_of_download)
    block = 10
    step = math.radians(block * 0.06 / 3600)
    lo = math.radians(-0.08)
    grid = np.zeros((3, 960, 960), dtype=np.uint8)
    rec = wht.wht_tile("A1", grid, step, lo, block, tmp_path)
    assert rec["coarse_cells_mapped"] == 75
    assert rec["fits_shape"] == [100, 100]
    assert 6900 <= int(grid[0].sum()) <= 7800
    assert max(abs(v) for v in rec["grid_area_conservation_error_fraction"]) < 0.03
    assert np.array_equal(grid[0], grid[1])
    assert np.array_equal(grid[1], grid[2])


def test_extracted_image_uses_primary_header_wcs(tmp_path, monkeypatch):
    """An extracted extension can have image data but only primary WCS."""
    hdr = WCS(naxis=2)
    hdr.wcs.crpix = [50.0, 50.0]
    hdr.wcs.crval = [150.12, 2.20]
    hdr.wcs.cdelt = [-0.06 / 3600, 0.06 / 3600]
    hdr.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    raw = tmp_path / "primary_wcs.fits"
    fits.HDUList([
        fits.PrimaryHDU(header=hdr.to_header()),
        fits.ImageHDU(np.ones((100, 100), dtype=np.float32), name="WHT"),
    ]).writeto(raw)
    gz = tmp_path / "primary_wcs.fits.gz"
    with open(raw, "rb") as fin, gzip.open(gz, "wb") as fout:
        shutil.copyfileobj(fin, fout)
    monkeypatch.setattr(wht, "download",
                        lambda url, dest: shutil.copyfile(gz, dest))
    grid = np.zeros((3, 960, 960), dtype=np.uint8)
    rec = wht.wht_tile("A1", grid, math.radians(0.6 / 3600),
                       math.radians(-0.08), 10, tmp_path)
    assert rec["wcs_source"] == "primary"
    assert rec["coarse_cells_mapped"] == 100


def test_missing_celestial_wcs_fails_with_diagnostics(tmp_path, monkeypatch):
    raw = tmp_path / "no_wcs.fits"
    fits.PrimaryHDU(np.ones((100, 100), dtype=np.float32)).writeto(raw)
    gz = tmp_path / "no_wcs.fits.gz"
    with open(raw, "rb") as fin, gzip.open(gz, "wb") as fout:
        shutil.copyfileobj(fin, fout)
    monkeypatch.setattr(wht, "download",
                        lambda url, dest: shutil.copyfile(gz, dest))
    grid = np.zeros((3, 960, 960), dtype=np.uint8)
    with fits.open(raw) as hdul:
        with pytest.raises(ValueError, match="HDU header diagnostics"):
            wht.select_celestial_wcs(hdul, hdul[0], "A1")
    # Do not call the production i2d header-prefix network fallback here;
    # its gzip/header parser is separately mocked in test_pm_wcs_probe.py.


def test_official_dr1_a1_corners_reconstruct_missing_wht_wcs():
    # Real A1 WHT FITS has 9600 x 12455 pixels and zero WCS cards.
    derived, meta = wht.official_corner_wcs("A1", nx=9600, ny=12455)
    assert derived.has_celestial
    assert meta["wcs_corner_residual_arcsec_max"] < 1
    assert all(0.055 < x < 0.065 for x in meta["wcs_pixel_scales_arcsec"])
    assert "NOT_original_header" in meta["wcs_source"]


def test_corner_wcs_rejects_inconsistent_image_shape():
    with pytest.raises(ValueError, match="pixel scales"):
        wht.official_corner_wcs("A1", nx=100, ny=100)
