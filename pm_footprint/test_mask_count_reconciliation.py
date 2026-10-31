"""Offline tests for worst-case JWST mask area vs candidate count reconciliation."""
import random
import unittest

from pm_footprint.mask_count_reconciliation import reconcile, analyze


class TestMaskCountReconciliation(unittest.TestCase):
    def test_365day_entry_regression(self):
        d = reconcile(0.31462445722188803, 0.528823684998808,
                      62, 63, 0.014199219444444447)
        self.assertAlmostEqual(d["unmasked_density_ratio_red_over_nonred"], 0.6700027133801733)
        self.assertAlmostEqual(d["masked_density_ratio_absolute_bounds"][0], 0.6255883816648122)
        self.assertAlmostEqual(d["masked_density_ratio_absolute_bounds"][1], 0.7016695454547341)
        self.assertAlmostEqual(d["red_only_mask_area_required_to_equalize_densities_deg2"],
                               0.10382521718745882)
        self.assertFalse(d["can_star_mask_alone_equalize_densities"])

    def test_365day_proximity_groups(self):
        d = reconcile(0.31462445722188803, 0.528823684998808,
                      60, 55, 0.014199219444444447)
        self.assertAlmostEqual(d["red_only_mask_area_required_to_equalize_densities_deg2"],
                               0.08095257237433895)
        self.assertFalse(d["can_star_mask_alone_equalize_densities"])

    def test_extrema_cover_random_admissible_allocations(self):
        R, T, M = 6.0, 10.0, 1.0
        z = reconcile(R, T, 3, 4, M)
        rng = random.Random(31415)
        for _ in range(5000):
            x = M*rng.random()
            y = (M-x)*rng.random()
            actual = (3/(R-x))/(4/(T-R-y))
            self.assertLessEqual(z["masked_density_ratio_absolute_bounds"][0]-1e-12, actual)
            self.assertLessEqual(actual, z["masked_density_ratio_absolute_bounds"][1]+1e-12)

    def test_uniform_counts_already_equal(self):
        d = reconcile(6, 10, 6, 4, 0.5)
        self.assertAlmostEqual(d["unmasked_density_ratio_red_over_nonred"], 1.0)
        self.assertIsNone(d["red_only_mask_area_required_to_equalize_densities_deg2"])
        self.assertTrue(d["can_star_mask_alone_equalize_densities"])

    def test_bad_inputs(self):
        for args in [(6,10,2,2,4), (6,10,2,2,-1), (11,10,2,2,0),
                     (6,10,0,2,1), (6,10,2.0,2,1)]:
            with self.assertRaises(ValueError):
                reconcile(*args)

    def test_12_scenarios_and_mask_center_guard(self):
        wht = {"area_deg2": {}, "candidate_counts": {}}
        for t in ("0", "0.1", "0.5"):
            wht["area_deg2"][t] = {"archive_f444w": 0.53,
                **{f"f444w_plus_red_{d}d": a for d,a in
                   [(30,.35),(180,.34),(365,.315),(730,.213)]}}
            wht["candidate_counts"][t] = {"archive_f444w": 125,
                **{f"f444w_plus_red_{d}d": n for d,n in
                   [(30,72),(180,72),(365,62),(730,42)]}}
        mask = {"n_mask_fits":20, "candidate_centers_in_code1":[],
                "sum_per_tile_mask_code1_area_deg2_NOT_UNION":.0142}
        counts = {"positional_groups_0p5arcsec":{"n_inside":60,"n_outside":55}}
        result = analyze(wht,mask,counts)
        self.assertEqual(sum(map(len,result["wht_threshold_baseline_days"].values())),12)
        self.assertFalse(result["wht_threshold_baseline_days"]["0"]["365"]["can_star_mask_alone_equalize_densities"])
        mask["candidate_centers_in_code1"]=[387132]
        with self.assertRaises(ValueError):
            analyze(wht,mask,counts)


if __name__ == "__main__":
    unittest.main()
