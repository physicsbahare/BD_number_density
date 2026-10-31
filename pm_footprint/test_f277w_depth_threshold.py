"""Offline tests for the F277W local-depth versus hybrid ratio threshold."""
import unittest
from pm_footprint.f277w_depth_threshold import m5_from_sigma, required_m5, summarize


class F277WDepthThresholdTests(unittest.TestCase):
    def test_nominal_faint_limit(self):
        self.assertAlmostEqual(required_m5(26.2), 28.897953115119062, places=10)
        self.assertGreater(28.9, required_m5(26.2))

    def test_flux_magnitude_equivalence(self):
        for m444 in (22, 24.5, 26.2):
            f444 = 10 ** ((31.4 - m444) / 2.5)
            for sigma in (0.5, 2, 10, 30):
                self.assertEqual(f444 / (2 * sigma) > 30,
                                 m5_from_sigma(sigma) > required_m5(m444))

    def test_invalid_sigma(self):
        for sigma in (-1, 0, float("nan")):
            with self.assertRaises(ValueError):
                m5_from_sigma(sigma)

    def test_frozen_roster_guard(self):
        with self.assertRaises(ValueError):
            summarize([])


if __name__ == "__main__":
    unittest.main()
