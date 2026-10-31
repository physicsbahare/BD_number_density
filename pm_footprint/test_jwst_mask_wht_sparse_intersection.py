"""Synthetic tests for direct JWST mask × WHT × red-epoch estimator."""
import tempfile
import unittest
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from pm_footprint.exact_geometry import lambert
from pm_footprint.jwst_mask_wht_sparse_intersection import (
    lambert_inverse, sparse_wht_positive, estimate
)


class TestMaskWHTSparseIntersection(unittest.TestCase):
    def test_lambert_inverse_roundtrip(self):
        ra = np.array([150.12, 150.4, 149.8, 150.0])
        dec = np.array([2.2, 2.1, 2.4, 2.0])
        x,y = lambert(ra,dec)
        a,b = lambert_inverse(x,y)
        np.testing.assert_allclose(a,ra,atol=1e-10)
        np.testing.assert_allclose(b,dec,atol=1e-10)

    def test_sparse_primary_wht_uses_correct_big_endian_offsets(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d)/"WHT.fits"
            w = WCS(naxis=2)
            w.wcs.crpix=[5.,5.]
            w.wcs.crval=[150.12,2.2]
            w.wcs.cdelt=[-0.03/3600,0.03/3600]
            w.wcs.ctype=["RA---TAN","DEC--TAN"]
            a = np.zeros((9,9),dtype=np.float32)
            a[4,4]=12.5
            a[5,6]=np.nan
            fits.PrimaryHDU(a,header=w.to_header()).writeto(path)
            ra,dec=w.all_pix2world([4,0,6,30],[4,0,5,30],0)
            positive,qc=sparse_wht_positive(path,ra,dec)
            np.testing.assert_array_equal(positive,[True,False,False,False])
            self.assertEqual(qc["n_inside_native_rectangle"],3)
            self.assertEqual(qc["bytes_of_native_WHT_pixels_read"],12)
            self.assertEqual(qc["n_nonfinite_inside"],1)

    def test_nested_MC_covariance_and_known_ratio(self):
        n=1000
        archive=np.ones(n,dtype=bool)
        red=np.zeros(n,dtype=bool)
        red[:600]=True
        wht=np.ones(n,dtype=bool)
        wht[900:]=False
        d=estimate(0.01,archive,red,wht,0.50,0.30)
        self.assertAlmostEqual(d["direct_mask_WHT_intersection_MC_deg2"],0.009)
        self.assertAlmostEqual(d["direct_mask_WHT_red365_intersection_MC_deg2"],0.006)
        self.assertAlmostEqual(d["corrected_red_fraction"],0.294/0.491)
        self.assertGreater(d["MC_95pct_halfwidth_red_fraction"],0)
        self.assertLess(d["MC_95pct_halfwidth_red_fraction"],0.01)

    def test_invalid_nesting_fails(self):
        n=100
        archive=np.zeros(n,bool)
        red=np.ones(n,bool)
        wht=np.ones(n,bool)
        with self.assertRaises(ValueError):
            estimate(.01,archive,red,wht,.5,.3)


if __name__=="__main__":
    unittest.main()
