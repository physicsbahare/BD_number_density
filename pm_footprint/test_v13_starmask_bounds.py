"""Offline tests of exact algebraic set bounds and MC sensitivity semantics."""
import unittest
from pm_footprint.v13_starmask_bounds import bounds, analyze


class TestConditionalStarMaskBounds(unittest.TestCase):
    def test_if_wht_fills_archive_bounds_collapse(self):
        z = bounds(10.0, 6.0, 10.0, 6.0, 2.0, 1.0)
        self.assertEqual(z["masked_red_wht_area_bounds_deg2"], [1.0,1.0])
        self.assertEqual(z["masked_nonred_wht_area_bounds_deg2"], [1.0,1.0])
        self.assertAlmostEqual(z["geometric_red_fraction_after_mask_bounds"][0],5/8)
        self.assertAlmostEqual(z["geometric_red_fraction_after_mask_bounds"][1],5/8)

    def test_two_independent_red_nonred_allocations(self):
        # Within survey S=10: red 6, nonred 4.
        # WHT covers red 5, nonred 3: up to 1 deg2 of masked
        # red and 1 of masked nonred could lie outside WHT.
        z = bounds(10,6,8,5,2,1)
        self.assertEqual(z["masked_red_wht_area_bounds_deg2"],[0,1])
        self.assertEqual(z["masked_nonred_wht_area_bounds_deg2"],[0,1])
        self.assertAlmostEqual(z["geometric_red_fraction_after_mask_bounds"][0],4/7)
        self.assertAlmostEqual(z["geometric_red_fraction_after_mask_bounds"][1],5/7)

    def test_reject_inconsistent_areas(self):
        with self.assertRaises(ValueError):
            bounds(1,0.6,0.8,0.7,0.1,0.05)  # red WHT exceeds archive red
        with self.assertRaises(ValueError):
            bounds(1,0.6,0.8,0.5,0.1,0.2)  # mask red exceeds mask survey

    def test_mc_sensitivity_broader_than_nominal(self):
        o={"exact_archive_deg2":10.0,"exact_red_epoch_deg2":6.0,
           "phase_results":[{"phase":0.25,"wht_positive_area_deg2":8.0,
             "area_intersections_deg2":{"archive":8.0,"red_365d":5.0}}]}
        m={"measured_mask_intersections_deg2":{
            "mask_within_archive_f444w":2.0,
            "mask_within_archive_red_365d":1.0,
            "survey_mc_95pct_halfwidth_deg2":0.05,
            "red_365d_mc_95pct_halfwidth_deg2":0.1}}
        result=analyze(o,m)
        lo,hi=result["nominal_mc_mask_fraction_bounds_across_phases"]
        slo,shi=result["illustrative_mc_marginal_95pct_corner_sensitivity_envelope_across_phases"]
        self.assertLessEqual(slo,lo)
        self.assertGreaterEqual(shi,hi)


if __name__=="__main__":
    unittest.main()
