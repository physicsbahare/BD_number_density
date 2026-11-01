"""Offline regression tests for the frozen v13 F444W-only quality-flag audit."""
import csv
import math
import unittest
from pathlib import Path
from .dropout_quality_flag_overlap import audit, K


class TestDropoutFlagOverlap(unittest.TestCase):
    def setUp(self):
        fixture = Path(__file__).with_name("results") / "v13_epoch_candidate_photometry_fixture_20261031.csv"
        with fixture.open(newline="") as f:
            self.rows = list(csv.DictReader(f))

    def test_snr_magnitude_error_boundary(self):
        self.assertAlmostEqual(K / .2, 5.428681023790647, places=11)
        self.assertGreater(K / 5, .2)
        self.assertLess(K / 5.5, .2)

    def test_frozen_candidate_counts_and_exception(self):
        r = audit(self.rows)
        self.assertEqual(r["n_candidates"], 125)
        self.assertEqual(r["n_passing_f444w_half_of_warn4"], 124)
        self.assertEqual(r["f444w_half_exception_ids"], [9095])
        self.assertEqual(r["epoch_counts"]["red365"]["f444_mag_error_lt_0p2"], 62)
        self.assertEqual(r["epoch_counts"]["nonred365"]["f444_mag_error_lt_0p2"], 62)

    def test_strong_negative_flux_is_review_not_rejection(self):
        r = audit(self.rows)
        x = {(a["candidate_id"], a["band"]) for a in r["negative_blue_snr_below_minus2"]}
        self.assertEqual(x, {(142367, "F115W"), (203707, "F150W")})
        self.assertEqual(r["n_unique_candidates_with_negative_blue_snr_below_minus2"], 2)

    def test_confirmed_movers_pass_bright_condition_only(self):
        r = audit(self.rows)
        self.assertEqual(set(r["confirmed_movers"]), {"387132", "1651050"})
        self.assertTrue(all(x["passes_f444w_half_of_warn4"]
                            for x in r["confirmed_movers"].values()))

    def test_bad_input_fails_closed(self):
        with self.assertRaises(ValueError):
            audit(self.rows[:-1])
        broken = [dict(row) for row in self.rows]
        broken[0]["new_e_F277W"] = "0"
        with self.assertRaises(ValueError):
            audit(broken)


if __name__ == "__main__":
    unittest.main()
