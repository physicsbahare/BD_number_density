"""python -m unittest pm_footprint.test_official_mask_placement_bounds"""
import unittest
from .official_mask_placement_bounds import fraction_bounds


class TestMaskBounds(unittest.TestCase):
    def test_zero_mask(self):
        self.assertEqual(fraction_bounds(1, .6, 0), (.6, .6))

    def test_red_or_nonred_extremes(self):
        low, high = fraction_bounds(1, .6, .1)
        self.assertAlmostEqual(low, .5/.9)
        self.assertAlmostEqual(high, .6/.9)

    def test_large_cap_clamps_strata(self):
        self.assertEqual(fraction_bounds(1, .6, .7), (0, 1))

    def test_invalid(self):
        with self.assertRaises(ValueError):
            fraction_bounds(1, 1.2, .1)

    def test_frozen_official_values(self):
        low, high = fraction_bounds(.528823684998808,
                                   .31462445722188803, .0141992195)
        self.assertAlmostEqual(low, .5837756614052463, places=12)
        self.assertAlmostEqual(high, .6113670808808773, places=12)


if __name__ == "__main__":
    unittest.main()
