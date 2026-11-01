"""Synthetic and frozen-repository tests for geometric tradeoff."""
import json
import unittest
from pathlib import Path
from .coverage_baseline_tradeoff import evaluate


class TestCoverageBaselineTradeoff(unittest.TestCase):
    def test_frozen_json_if_present(self):
        a=Path('pm_footprint/results/wht60/summary.json')
        b=Path('pm_footprint/results/exact/summary.json')
        if not (a.exists() and b.exists()):
            self.skipTest('Frozen GitHub fixtures not available in standalone checkout')
        r=evaluate(json.loads(a.read_text()),json.loads(b.read_text()))
        self.assertAlmostEqual(r['threshold_baseline_table'][0]['red_epoch_fraction'],.5949515238195074)
        self.assertEqual(r['730d_vs_365d_positive_wht']['red_candidates_lost'],20)
        self.assertEqual(r['730d_vs_365d_positive_wht']['confirmed_pm_ids_at_730d'],[387132])
        self.assertLess(abs(r['wht_50pct_median_cut_vs_positive']['red365_fraction_change_percentage_points']),.2)

    def test_roster_hash_mismatch_fails(self):
        with self.assertRaisesRegex(ValueError, 'roster SHA mismatch'):
            evaluate({'candidate_catalog_sha256':'a'}, {'candidate_catalog_sha256':'b'})


if __name__=='__main__':
    unittest.main()
