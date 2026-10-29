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
    assert int(grid[0].sum()) == 75 * 100
    assert int(grid[1].sum()) == 75 * 100
    assert int(grid[2].sum()) == 75 * 100
