"""Offline tests for hybrid F277W signed-SNR/ratio branch mathematics."""
import unittest
from pm_footprint.hybrid_f277w_selection_probability import (
    hybrid_pass_probability, normal_cdf, summarize)


class HybridF277WProbabilityTests(unittest.TestCase):
    def test_zero_flux_below_critical(self):
        self.assertAlmostEqual(hybrid_pass_probability(120,1.99),
                               normal_cdf(2),places=12)

    def test_exact_boundary_has_half_probability_jump(self):
        self.assertAlmostEqual(hybrid_pass_probability(120,2),
                               normal_cdf(2)-0.5,places=12)
        self.assertAlmostEqual(hybrid_pass_probability(120,1.999999999),
                               normal_cdf(2),places=12)

    def test_zero_flux_above_critical(self):
        self.assertAlmostEqual(hybrid_pass_probability(120,4),
                               normal_cdf(1)-0.5,places=12)

    def test_nonzero_true_flux(self):
        mu=1.5
        self.assertAlmostEqual(hybrid_pass_probability(120,1.5,mu),
                               normal_cdf(2-mu/1.5),places=12)
        self.assertAlmostEqual(hybrid_pass_probability(120,3,mu),
                               normal_cdf((4-mu)/3)-normal_cdf(-mu/3),
                               places=12)

    def test_invalid_inputs(self):
        for a,b,c in [(0,1,0),(1,0,0),(1,-1,0),(1,1,float("nan"))]:
            with self.assertRaises(ValueError):
                hybrid_pass_probability(a,b,c)

    def test_frozen_roster_guard(self):
        with self.assertRaises(ValueError):
            summarize([{"candidate_id":1}])

    def test_actual_committed_fixture_if_present(self):
        import csv
        from pathlib import Path
        path=Path("pm_footprint/results/v13_epoch_candidate_photometry_fixture_20261031.csv")
        if not path.exists():
            self.skipTest("Run from repository root to check actual frozen fixture")
        with path.open(newline="") as f:
            r=summarize(list(csv.DictReader(f)))
        self.assertEqual(r["all_candidates"]["n"],125)
        self.assertEqual(r["red365"]["n"],62)
        self.assertEqual(r["all_candidates"]["negative_branch_allowed_n"],57)
        self.assertEqual(r["all_candidates"]["observed_nonpositive_F277W_n"],26)
        self.assertEqual(r["no_red365"]["negative_branch_allowed_n"],37)
        self.assertEqual(r["red365"]["negative_branch_allowed_n"],20)


if __name__=="__main__":
    unittest.main()
