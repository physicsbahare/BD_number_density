"""Regression tests for PC/TAN WCS parity and pixel area.

Run: python -m unittest pm_footprint.test_mask_wcs -v
No JWST archive, WHT files, network, or astropy required.
"""
import unittest
import numpy as np
from .mask_wcs import (cd_from_header, pixel_to_sky, pixel_to_lambert,
                       tan_pixel_area_deg2, legacy_column_scaling_error_arcsec)


def header(angle=0.0, sx=-0.03 / 3600, sy=0.03 / 3600):
    a = np.deg2rad(angle)
    return dict(CRPIX1=9600.5, CRPIX2=12455.5,
                CRVAL1=150.12, CRVAL2=2.2,
                CDELT1=sx, CDELT2=sy,
                PC1_1=np.cos(a), PC1_2=-np.sin(a),
                PC2_1=np.sin(a), PC2_2=np.cos(a))


class TestPCConvention(unittest.TestCase):
    def test_row_scaling_with_opposite_axis_signs(self):
        h = header(30)
        cd = cd_from_header(h)
        self.assertAlmostEqual(cd[0, 1], h["CDELT1"] * h["PC1_2"])
        self.assertAlmostEqual(cd[1, 0], h["CDELT2"] * h["PC2_1"])
        self.assertNotAlmostEqual(cd[0, 1], h["CDELT2"] * h["PC1_2"],
                                  places=12)

    def test_inverse_tan_roundtrip_independent_forward_projection(self):
        for angle in (0, 12, 30, -45):
            h = header(angle)
            x = np.array([0., 9600., 18000.])
            y = np.array([0., 12000., 24000.])
            ra, dec = pixel_to_sky(x, y, h)
            dra = (ra - np.deg2rad(h["CRVAL1"]) + np.pi) % (2*np.pi) - np.pi
            d0 = np.deg2rad(h["CRVAL2"])
            den = (np.sin(d0) * np.sin(dec)
                   + np.cos(d0) * np.cos(dec) * np.cos(dra))
            tangent = np.vstack((
                np.cos(dec) * np.sin(dra) / den,
                (np.cos(d0) * np.sin(dec)
                 - np.sin(d0) * np.cos(dec) * np.cos(dra)) / den
            )) * (180 / np.pi)
            offset = np.linalg.solve(cd_from_header(h), tangent)
            np.testing.assert_allclose(offset[0] + h["CRPIX1"] - 1, x,
                                       atol=1e-6)
            np.testing.assert_allclose(offset[1] + h["CRPIX2"] - 1, y,
                                       atol=1e-6)

    def test_lambert_center(self):
        h = header(30)
        x, y = pixel_to_lambert(h["CRPIX1"]-1, h["CRPIX2"]-1, h)
        self.assertAlmostEqual(float(x), 0, places=9)
        self.assertAlmostEqual(float(y), 0, places=9)

    def test_tan_area_and_projection(self):
        h = header(30)
        x, y = h["CRPIX1"]-1, h["CRPIX2"]-1
        self.assertAlmostEqual(float(tan_pixel_area_deg2(x, y, h)),
                               abs(np.linalg.det(cd_from_header(h))),
                               places=18)
        self.assertLess(float(tan_pixel_area_deg2(x+10000, y+10000, h)),
                        float(tan_pixel_area_deg2(x, y, h)))

    def test_legacy_error_zero_for_diagonal_pc(self):
        h = header(0)
        self.assertEqual(float(legacy_column_scaling_error_arcsec(h, 0, 0)),
                         0.0)

    def test_legacy_error_large_for_rotated_opposite_axis_scales(self):
        h = header(30)
        # Synthetic WCS only: this is NOT a measured COSMOS2025 offset.
        x, y = h["CRPIX1"]-1+1000, h["CRPIX2"]-1+1000
        self.assertAlmostEqual(
            float(legacy_column_scaling_error_arcsec(h, x, y)),
            42.426406871, places=5)

    def test_invalid_matrix_fails_closed(self):
        h = header()
        h["PC2_2"] = 0.
        h["PC1_1"] = 0.
        with self.assertRaises(ValueError):
            cd_from_header(h)


if __name__ == "__main__":
    unittest.main()
