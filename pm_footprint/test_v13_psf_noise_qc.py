"""Offline synthetic tests for small PSF/EE and corrected-noise audit."""
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
from astropy.io import fits

from pm_footprint.v13_psf_noise_qc import (
    read_small_psf, enclosed_fraction, notebook_factors, run
)


class TestSmallPSFNoiseQC(unittest.TestCase):
    def test_psf_read_and_subpixel_integral(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"psf.fits"
            a=np.zeros((11,11),dtype=np.float32)
            a[5,5]=1
            fits.PrimaryHDU(a).writeto(p)
            x=read_small_psf(p)
            self.assertEqual(x.shape,(11,11))
            self.assertAlmostEqual(enclosed_fraction(x,.03),1.0)
            self.assertAlmostEqual(enclosed_fraction(x,.3),1.0)

    def test_noise_factor_literal_and_required_equation(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"n.ipynb"
            factors={b:[1.5,2,2.5,3,3.5]
                     for b in ("F115W","F150W","F277W","F444W")}
            src=("pixel_scale=0.03\nnorm_radius=1\n"
                 +"corr_factors = "+repr(factors)+"\n"
                 +"sp[f'new_e_{band}'] = sp[f'e_{band}'] * corr\n")
            p.write_text(json.dumps({"cells":[{"cell_type":"code","source":[src]}]}))
            self.assertEqual(notebook_factors(p)["F444W"][0],1.5)
            p.write_text(json.dumps({"cells":[{"cell_type":"code",
                                                "source":["corr_factors = "+repr(factors)]}]}))
            with self.assertRaises(ValueError):
                notebook_factors(p)

    def test_end_to_end_small_inputs(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d)
            a=np.zeros((101,101),dtype=np.float32)
            y,x=np.indices(a.shape)
            a[:]=np.exp(-((x-50)**2+(y-50)**2)/(2*5**2))
            a/=a.sum()
            fits.PrimaryHDU(a).writeto(d/"psf.fits")
            radii=np.array([.01,.1,.2,.5,1.,1.5])
            with (d/"ee.txt").open("w") as f:
                f.write("aper_radius F444W\n")
                for r in radii:
                    f.write(f"{r} {1-np.exp(-r*10):.8f}\n")
            factors={b:[2,2,2,2,2]
                     for b in ("F115W","F150W","F277W","F444W")}
            src=("pixel_scale=0.03\nnorm_radius=1\n"
                 +"corr_factors = "+repr(factors)+"\n"
                 +"sp[f'new_e_{band}'] = sp[f'e_{band}'] * corr\n")
            (d/"n.ipynb").write_text(json.dumps({
                "cells":[{"cell_type":"code","source":[src]}]}))
            r=run(d/"psf.fits",d/"ee.txt",d/"n.ipynb",d/"out.json")
            self.assertEqual(len(r["curve_of_growth_F444W"]),6)
            self.assertAlmostEqual(
                r["empirical_error_factor_by_aperture"]["F277W"][0]["nominal_depth_overstatement_mag_if_new_e_used"],
                2.5*np.log10(2))
            self.assertTrue((d/"out.json").exists())


if __name__=="__main__":
    unittest.main()
