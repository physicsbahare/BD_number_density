"""Synthetic tests: no giant FITS, network, or MAST required."""
import unittest
from pm_footprint.candidate_epoch_photometry_contrast import (
    groups, passes, permutation, median, load, analyze
)


class TestEpochPhotometryContrast(unittest.TestCase):
    def test_median_even_odd(self):
        self.assertEqual(median([3, 1, 2]), 2)
        self.assertEqual(median([4, 2, 1, 3]), 2.5)

    def test_group_transitivity(self):
        rows = [{"ra_deg": 150 + i * 0.00010, "dec_deg": 2.,
                 "red365": True} for i in range(3)]
        self.assertEqual(len(groups(rows)), 1)
        self.assertEqual(len(groups(rows, radius_arcsec=0.1)), 3)

    def test_negative_f277_2sigma_and_strict_ratio(self):
        r = {"f_F444W": 150., "new_e_F444W": 5.,
             "f_F277W": -3., "new_e_F277W": 2.5,
             "f_F115W": 0., "new_e_F115W": 5.,
             "f_F150W": 0., "new_e_F150W": 5.,
             "mag_F444W": 25.}
        self.assertFalse(passes(r))  # exactly 150/(2*2.5)=30
        r["f_F444W"] = 151.
        self.assertTrue(passes(r))
        self.assertFalse(passes(r, f277_offset_njy=10.))

    def test_seed_reproducible(self):
        rows = [{"red365": i < 3, "mag_F444W": float(i)}
                for i in range(6)]
        a = permutation(rows, reps=100, seed=12)
        self.assertEqual(a, permutation(rows, reps=100, seed=12))
        self.assertEqual(a["red_n"], 3)

    def test_real_small_fixture_invariants(self):
        from pathlib import Path
        p = Path("pm_footprint/results/v13_epoch_candidate_photometry_fixture_20261031.csv")
        if not p.exists():
            self.skipTest("Run from repository root")
        rows = load(p)
        self.assertEqual(len(rows), 125)
        self.assertEqual(sum(r["red365"] for r in rows), 62)
        self.assertTrue(all(passes(r) for r in rows))
        result = analyze(rows, reps=10)
        self.assertEqual(result["positional_groups_0p5arcsec"], 115)
        self.assertEqual(result["red_epoch_groups"], 60)
        self.assertEqual(result["selection_sensitivity"]["F277W_flux_plus_1nJy"]["passing"], 109)


if __name__ == "__main__":
    unittest.main()
