import unittest
from pm_footprint.spatial_count_significance import conditional_test

class TestSpatialSignificance(unittest.TestCase):
    def test_symmetric(self):
        d = conditional_test(5, 10, 1, 2)
        self.assertAlmostEqual(d["conditional_exact_two_sided_pvalue"], 1)
        self.assertAlmostEqual(d["inside_outside_density_ratio"], 1)
    def test_frozen_entries(self):
        d = conditional_test(62, 125, .31462445722188803, .528823684998808)
        self.assertAlmostEqual(d["inside_outside_density_ratio"], .6700027133801733)
        self.assertAlmostEqual(d["conditional_exact_two_sided_pvalue"], .02844524875452935)
    def test_invalid(self):
        for args in [(10,9,1,2), (1,3,2,1), (1,0,1,2)]:
            with self.assertRaises(ValueError): conditional_test(*args)
    def test_grouped_ci(self):
        d = conditional_test(60, 115, .31462445722188803, .528823684998808)
        lo, hi = d["transformed_95pct_density_ratio_ci"]
        self.assertLess(lo, d["inside_outside_density_ratio"])
        self.assertGreater(hi, d["inside_outside_density_ratio"])
        self.assertGreater(d["conditional_exact_two_sided_pvalue"], .05)
if __name__ == "__main__":
    unittest.main()
