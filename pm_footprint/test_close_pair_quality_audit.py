"""Offline source-proximity audit tests, no survey FITS required."""
import csv
from pathlib import Path
import unittest
from pm_footprint.close_pair_quality_audit import audit, connected_components, separation_arcsec


class ClosePairQualityTests(unittest.TestCase):
    def test_separation_arcsec(self):
        a={"ra":150.,"dec":2.}
        b={"ra":150.,"dec":2.+0.4/3600}
        self.assertAlmostEqual(separation_arcsec(a,b),0.4,places=6)

    def test_transitive_components(self):
        a=[{"ra":150.,"dec":2.},
           {"ra":150.,"dec":2.+0.3/3600},
           {"ra":150.,"dec":2.+0.6/3600}]
        groups, edges=connected_components(a,0.5)
        self.assertEqual(len(groups),1)
        self.assertEqual(len(edges),2)

    def test_invalid_radius(self):
        with self.assertRaises(ValueError):
            connected_components([],0)

    def test_roster_guard(self):
        with self.assertRaises(ValueError):
            audit([])

    def test_frozen_fixture(self):
        path=Path("pm_footprint/results/v13_epoch_candidate_photometry_fixture_20261031.csv")
        if not path.exists():
            self.skipTest("Frozen fixture only available in repository root")
        with path.open(newline="") as f:
            x=audit(csv.DictReader(f))
        c=x["close_0p5arcsec"]
        self.assertEqual((c["n_multi_components"],c["n_close_members"],c["n_isolated"]),(9,19,106))
        self.assertEqual((c["bad_kron_close"],c["bad_kron_isolated"]),(17,43))
        self.assertEqual((c["red_close"],c["nonred_close"]),(4,15))
        self.assertEqual(x["radius_sweep"][-2]["n_components"],115)
        self.assertEqual(x["density_sensitivity"]["bright_m444_lt25_isolated_only"]["red"],23)
        self.assertEqual(x["density_sensitivity"]["bright_m444_lt25_isolated_only"]["nonred"],26)
        self.assertAlmostEqual(c["bad_kron_fisher_two_sided_exploratory_p"],8.730551621528084e-5,places=12)

if __name__=="__main__":unittest.main()
