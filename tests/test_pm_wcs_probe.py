"""Offline tests for streaming original SCI WCS without downloading image pixels."""
import gzip
import io

import numpy as np
import pytest

pytest.importorskip("astropy")
from astropy.io import fits
from astropy.wcs import WCS

from pm_footprint import probe_i2d_wcs as probe


def test_find_official_i2d_href_without_network():
    page = ('<html><a href="/data/nircam/mosaics/'
            'mosaic_nircam_f444w_COSMOS-Web_60mas_A1_v1.0_i2d.fits.gz">'
            'i2d</a></html>')
    url = probe.find_i2d_url("A1", page)
    assert url.endswith("A1_v1.0_i2d.fits.gz")
    with pytest.raises(ValueError):
        probe.find_i2d_url("A99", page)


def test_stream_original_sci_header_only(monkeypatch):
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [5, 5]
    wcs.wcs.crval = [150.12, 2.20]
    wcs.wcs.cdelt = [-0.06 / 3600, 0.06 / 3600]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    hdus = fits.HDUList([
        fits.PrimaryHDU(),
        fits.ImageHDU(np.arange(10000, dtype=np.float32).reshape(100, 100),
                      header=wcs.to_header(), name="SCI"),
    ])
    output = io.BytesIO()
    hdus.writeto(output)
    prefix = gzip.compress(output.getvalue())

    class Response(io.BytesIO):
        status = 206

    monkeypatch.setattr(probe, "urlopen",
                        lambda *args, **kwargs: Response(prefix))
    header, recovered, metadata = probe.i2d_sci_header("https://example.test/test.gz")
    assert header["NAXIS1"] == 100
    assert header["NAXIS2"] == 100
    assert recovered.has_celestial
    assert metadata["compressed_prefix_bytes"] == len(prefix)
