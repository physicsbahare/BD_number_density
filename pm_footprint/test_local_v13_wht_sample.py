"""Small synthetic FITS tests; no external hard drive or archive downloads."""
import csv
import os
from pathlib import Path
import tempfile
import unittest

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from pm_footprint.local_v13_wht_sample import (
    read_header_and_offset, sample_raw, tan_cell_areas, run
)


class LocalV13WhtTests(unittest.TestCase):
    def test_big_endian_direct_read(self):
        with tempfile.TemporaryFile() as f:
            f.write(np.arange(60, dtype=">f4").tobytes())
            f.flush()
            np.testing.assert_array_equal(
                sample_raw(f.fileno(), 0, 10, 3, [0, 2, 5, 9]),
                [30, 32, 35, 39]
            )

    def test_tan_pixel_area_at_reference(self):
        cd = np.diag([-0.03/3600, 0.03/3600])
        a = tan_cell_areas(np.array([0]), 0, np.array([1]), 1,
                           cd, np.array([1, 1]))[0]
        self.assertAlmostEqual(a, (0.03/3600)**2, places=17)

    def test_synthetic_fits_and_candidate(self):
        with tempfile.TemporaryDirectory() as d:
            base = Path(d)
            w = WCS(naxis=2)
            w.wcs.crpix = [8, 6]
            w.wcs.cdelt = [-0.03/3600, 0.03/3600]
            w.wcs.crval = [150.12, 2.20]
            w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
            file = base / "tiny.fits"
            fits.PrimaryHDU(np.ones((12, 16), dtype=np.float32),
                            header=w.to_header()).writeto(file)
            header, offset, wcs = read_header_and_offset(file)
            self.assertEqual(offset % 2880, 0)
            ra, dec = wcs.pixel_to_world_values(7, 5)
            cat = base / "candidates.csv"
            with open(cat, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["candidate_id", "ra_deg", "dec_deg"])
                writer.writerow(["387132", ra, dec])
            result = run(file, cat, 4, base / "out")
            self.assertEqual(result["candidate_centers_wht_positive"], 1)
            self.assertEqual(result["candidate_centers_inside_image"], 1)
            self.assertEqual(result["candidate_count"], 1)
            self.assertLess(result["phase_difference_area_deg2"], 1e-12)
            area = 12 * 16 * (0.03 / 3600)**2
            self.assertAlmostEqual(result["sampled_mean_area_deg2"], area, delta=area*0.0001)
            self.assertTrue((base/"out_grid.npz").exists())
            self.assertTrue((base/"out_candidates.csv").exists())


if __name__ == "__main__":
    unittest.main()
