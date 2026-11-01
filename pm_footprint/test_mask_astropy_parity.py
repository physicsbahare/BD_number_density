"""Independent WCS regression: test our mask TAN transform against Astropy."""
import importlib.util
import unittest
import numpy as np
from .mask_astropy_parity import validate_fits_header_wcs

HAS_ASTROPY = importlib.util.find_spec("astropy") is not None


def header(angle=31.0):
    a = np.deg2rad(angle)
    return {
        "NAXIS": 2, "NAXIS1": 19200, "NAXIS2": 24910,
        "CTYPE1": "RA---TAN", "CTYPE2": "DEC--TAN",
        "CUNIT1": "deg", "CUNIT2": "deg",
        "CRVAL1": 150.12, "CRVAL2": 2.2,
        "CRPIX1": 9600.5, "CRPIX2": 12455.5,
        "CDELT1": -0.03/3600, "CDELT2": 0.03/3600,
        "PC1_1": float(np.cos(a)), "PC1_2": float(-np.sin(a)),
        "PC2_1": float(np.sin(a)), "PC2_2": float(np.cos(a)),
    }


@unittest.skipUnless(HAS_ASTROPY, "Independent Astropy WCS dependency absent")
class MaskAstropyWcsTests(unittest.TestCase):
    def test_rotated_pc_opposite_cdelt_signs(self):
        for rotation in (0, 5, 31, -44, 90):
            with self.subTest(rotation=rotation):
                result = validate_fits_header_wcs(header(rotation), 19200, 24910)
                self.assertLess(result["max_angular_sep_arcsec"], 0.001)
                self.assertEqual(result["positions_checked"], 9)

    def test_tan_area_at_reference(self):
        result = validate_fits_header_wcs(header(), 19200, 24910)
        nominal = (0.03/3600)**2
        self.assertAlmostEqual(
            result["center_solid_angle_deg2_per_native_pixel"]/nominal, 1,
            places=6
        )

    def test_non_tan_projection_refused(self):
        h = header()
        h["CTYPE1"] = "RA---SIN"
        with self.assertRaisesRegex(ValueError, "undistorted"):
            validate_fits_header_wcs(h, 19200, 24910)

    def test_sip_distortion_refused(self):
        h = header()
        h["A_ORDER"] = 2
        with self.assertRaisesRegex(ValueError, "SIP"):
            validate_fits_header_wcs(h, 19200, 24910)

    def test_broken_cd_matrix_refused(self):
        h = header(0)
        h["PC1_1"] = 0.0
        h["PC2_2"] = 0.0
        with self.assertRaises((ValueError, RuntimeError)):
            validate_fits_header_wcs(h, 19200, 24910)


if __name__ == "__main__":
    unittest.main()
