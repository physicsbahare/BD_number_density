import unittest
from pm_footprint.candidate_spatial_diagnostic import exact_binomial_two_sided, positional_groups, diagnostic

class SpatialDiagnosticTests(unittest.TestCase):
    def test_exact_binomial(self):
        self.assertAlmostEqual(exact_binomial_two_sided(0, 2, .5), .5)
        self.assertAlmostEqual(exact_binomial_two_sided(1, 2, .5), 1.)
    def test_proximity_and_mixed(self):
        rows = [dict(candidate_id='a', ra_deg='150', dec_deg='2', in_archive_f444w='True', in_f444w_plus_red_365d='True'),
                dict(candidate_id='b', ra_deg='150.00002', dec_deg='2', in_archive_f444w='True', in_f444w_plus_red_365d='False'),
                dict(candidate_id='c', ra_deg='151', dec_deg='2', in_archive_f444w='True', in_f444w_plus_red_365d='False')]
        self.assertEqual(len(positional_groups(rows)), 2)
        self.assertEqual(diagnostic(rows, 2, 1)['mixed_coverage_groups'], [['a','b']])
    def test_bad_area(self):
        with self.assertRaises(ValueError): diagnostic([], 1, 2)

if __name__ == '__main__':
    unittest.main()
