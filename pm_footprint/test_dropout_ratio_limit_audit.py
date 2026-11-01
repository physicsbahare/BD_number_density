"""Regression tests for F277W dropout limit audit; stdlib only.

Run: python -m unittest pm_footprint.test_dropout_ratio_limit_audit
"""
import csv
import math
import unittest
from pathlib import Path

from .dropout_ratio_limit_audit import (
    AB_NJY_ZEROPOINT, audit, flux_njy, required_5sigma_f277_depth,
)

FIXTURE = Path(__file__).parent / "results" / "v13_epoch_candidate_photometry_fixture_20261031.csv"


class DropoutRatioLimitTests(unittest.TestCase):
    def test_flux_scale(self):
        self.assertAlmostEqual(flux_njy(AB_NJY_ZEROPOINT), 1.0)
        self.assertAlmostEqual(flux_njy(26.2), 120.22644346174131, places=9)

    def test_required_depth(self):
        self.assertAlmostEqual(required_5sigma_f277_depth(), 28.89795311511906, places=9)
        self.assertAlmostEqual(flux_njy(26.2) / 60.0, 2.0037740576956884, places=9)

    def test_fainter_depth_increases_ratio_certification(self):
        m444 = 26.2
        f = flux_njy(m444)
        sigma_shallow = 8.7
        sigma_deep = 1.9
        self.assertLess(f / (2 * sigma_shallow), 30)
        self.assertGreater(f / (2 * sigma_deep), 30)

    def test_negative_flux_is_not_an_automatic_color_failure(self):
        f444, f277, err277 = 120.0, -3.0, 5.0
        self.assertLess(f277 / err277, 2)
        self.assertEqual(f444 / (max(0, f277) + 2 * err277), 12.0)

    def test_real_frozen_sample(self):
        with FIXTURE.open(newline="") as f:
            result = audit(csv.DictReader(f))
        groups = result["groups"]
        self.assertEqual(groups["all"]["n"], 125)
        self.assertEqual(groups["red365"]["n"], 62)
        self.assertEqual(groups["nonred365"]["n"], 63)
        self.assertEqual(groups["f444w_25_to_26p2"]["n"], 59)
        self.assertEqual(groups["f444w_bright_lt25"]["n"], 66)
        self.assertEqual(groups["all"]["n_ratio_gt30_from_zero_centered_2sigma"], 57)
        self.assertEqual(groups["all"]["n_ratio_gt30_from_positive_flux_plus_2sigma"], 50)
        self.assertEqual(groups["f444w_25_to_26p2"]["n_ratio_gt30_from_zero_centered_2sigma"], 0)
        self.assertEqual(groups["red365"]["n_ratio_gt30_from_zero_centered_2sigma"], 20)
        self.assertEqual(groups["nonred365"]["n_ratio_gt30_from_zero_centered_2sigma"], 37)
        self.assertEqual(groups["all"]["f277_abs_snr_lt2"], 125)
        self.assertAlmostEqual(groups["all"]["median_f277_sigma_njy"], 8.681924948037663, places=9)
        self.assertAlmostEqual(groups["all"]["median_f277_5sigma_depth_ab"], 27.30603492114213, places=9)
        self.assertLess(result["confirmed_pm_objects"]["387132"]["ratio_lower_bound_zero_centered_2sigma"], 30)
        self.assertLess(result["confirmed_pm_objects"]["1651050"]["ratio_lower_bound_zero_centered_2sigma"], 30)

    def test_reject_inconsistent_mag_flux(self):
        with FIXTURE.open(newline="") as f:
            rows = list(csv.DictReader(f))
        rows[0]["mag_F444W"] = "20.0"
        with self.assertRaisesRegex(ValueError, "AB flux/magnitude mismatch"):
            audit(rows)

    def test_reject_duplicate_ids(self):
        with FIXTURE.open(newline="") as f:
            rows = list(csv.DictReader(f))
        rows[1]["candidate_id"] = rows[0]["candidate_id"]
        with self.assertRaisesRegex(ValueError, "125 unique IDs"):
            audit(rows)


if __name__ == "__main__":
    unittest.main()
