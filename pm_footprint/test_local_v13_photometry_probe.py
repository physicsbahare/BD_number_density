"""Synthetic tests for the local v13 photometry preflight; no real FITS needed."""
import json
import os
from pathlib import Path
import tempfile
import unittest

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from pm_footprint.local_v13_photometry_probe import (
    primary_header, read_catalog_schema, read_native_patch,
    summarize_patch, notebook_summary
)


class TestLocalV13PhotometryProbe(unittest.TestCase):
    def test_fits_header_and_sparse_err_pixels(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "F277W_psf_matched_F444W_ERR_v13.fits"
            w = WCS(naxis=2)
            w.wcs.crpix=[6,6]
            w.wcs.cdelt=[-0.03/3600,0.03/3600]
            w.wcs.crval=[150.12,2.2]
            w.wcs.ctype=["RA---TAN","DEC--TAN"]
            h = w.to_header()
            h["BUNIT"]="MJy/sr"
            a = np.arange(11*11,dtype=">f4").reshape((11,11))
            fits.PrimaryHDU(a,header=h).writeto(path)
            head,offset,info,wcs=primary_header(path)
            self.assertEqual(info["shape_yx"],[11,11])
            self.assertEqual(info["keywords"]["BUNIT"],"MJy/sr")
            self.assertTrue(info["has_celestial_wcs"])
            self.assertIsNotNone(wcs)
            self.assertEqual(offset%2880,0)
            fd=os.open(path,os.O_RDONLY)
            try:
                patch,inside=read_native_patch(fd,offset,head,5,5)
                edge,edgeinside=read_native_patch(fd,offset,head,0,0)
            finally:
                os.close(fd)
            self.assertTrue(inside)
            np.testing.assert_array_equal(patch,[48,49,50,59,60,61,70,71,72])
            self.assertFalse(edgeinside)
            self.assertIsNone(edge)
            stats=summarize_patch(patch)
            self.assertEqual(stats["center_err"],60)
            self.assertEqual(stats["positive_count"],9)

    def test_catalog_schema_no_data_scan(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d)/"cat.fits"
            c=[fits.Column(name="id",format="K",array=np.array([1,2],dtype=np.int64)),
               fits.Column(name="f_f444w",format="E",unit="uJy",
                           array=np.array([0.1,0.2],dtype=np.float32))]
            fits.HDUList([fits.PrimaryHDU(),
                          fits.BinTableHDU.from_columns(c)]).writeto(path)
            d=read_catalog_schema(path)
            self.assertEqual(d["rows"],2)
            self.assertEqual(d["column_count"],2)
            self.assertEqual(d["columns"][0]["name"],"id")
            self.assertEqual(d["columns"][1]["unit"],"uJy")

    def test_notebook_summary_has_no_source_text(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/"note.ipynb"
            path.write_text(json.dumps({"cells":[
                {"cell_type":"code","source":["PRIVATE_PATH = '/private/abc'\n",
                                               "flux_f277w = 1.0\n"]},
                {"cell_type":"markdown","source":["PSF correction"]}]}))
            r=notebook_summary(path)
            self.assertEqual(r["total_cells"],2)
            self.assertEqual(r["code_cells"],1)
            self.assertEqual(r["keyword_hits"][0]["cell_index"],0)
            self.assertIn("f277w",r["keyword_hits"][0]["keywords"])
            self.assertNotIn("PRIVATE_PATH",str(r))


if __name__=="__main__":
    unittest.main()
