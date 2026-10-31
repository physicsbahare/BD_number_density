import unittest
from pm_footprint.pm_baseline_feasibility import analyze, mu_threshold, sigma_requirement


class TestBaselineFeasibility(unittest.TestCase):
    def test_inverse(self):
        g = 752.1011143171309
        self.assertAlmostEqual(mu_threshold(g, sigma_requirement(g, 50)), 50)

    def test_analytic(self):
        self.assertAlmostEqual(mu_threshold(365.25, 20), 3 * (2**.5) * 20)

    def test_monotonic(self):
        self.assertLess(mu_threshold(800, 20), mu_threshold(400, 20))
        self.assertGreater(sigma_requirement(800, 50), sigma_requirement(400, 50))

    def test_input_rejected(self):
        with self.assertRaises(ValueError): mu_threshold(0, 20)
        with self.assertRaises(ValueError): sigma_requirement(100, 0)
        with self.assertRaises(ValueError): analyze([{'candidate_id': '1'}])

    def test_counts(self):
        rows = [{'candidate_id': str(i), 'f356w_ge365': 'false', 'best_f356w_gap_days': ''} for i in range(125)]
        rows[0].update({'f356w_ge365': 'true', 'best_f356w_gap_days': '365.25'})
        rows[1].update({'f356w_ge365': 'true', 'best_f356w_gap_days': '730.5'})
        result = analyze(rows)
        self.assertEqual(result['n_geometric_f356w_ge365'], 2)
        self.assertEqual(result['conditional_baseline_sufficient_counts_by_sigma_then_mu']['20']['50'], 1)
        self.assertEqual(result['conditional_baseline_sufficient_counts_by_sigma_then_mu']['20']['100'], 2)


if __name__ == '__main__':
    unittest.main()
