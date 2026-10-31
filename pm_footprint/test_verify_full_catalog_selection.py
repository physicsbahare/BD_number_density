"""Offline QA of compact user-provided whole-catalog selection evidence."""
import copy
import json
import unittest
from pathlib import Path

from pm_footprint.verify_full_catalog_selection import summarize

FIXTURE=Path("pm_footprint/results/v13_full_catalog_paper_stages_20261031.json")


class FullCatalogProofTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=json.loads(FIXTURE.read_text())

    def test_all_manuscript_stages_and_exact_frozen_roster_match(self):
        r=summarize(self.data)
        self.assertEqual(r["status"],"FULL_3P28M_V13_SELECTION_REPRODUCED")
        self.assertEqual(r["reported_vs_reproduced_stage_counts"]["new_e_with_use_phot"],
                         [2997,335,125])
        self.assertEqual(r["additional_candidates_excluded_by_use_phot_gate_per_stage"],
                         [1167,1047,1039])
        self.assertEqual(r["numeric_rows_excluded_by_use_phot_gate"],69622)
        self.assertEqual(r["candidate_ID_recovery_by_scenario"]["new_e_use_phot"],
                         {"in_frozen":125,"unexpected":0,"missing_frozen":0})
        self.assertAlmostEqual(
            r["fraction_of_ungated_corrected_final_retained_after_gate"],
            125/1164)

    def test_unfiltered_new_e_does_not_reproduce_manuscript(self):
        r=summarize(self.data)
        self.assertEqual(r["stage_counts_by_scenario"]["new_e_all"],[4164,1382,1164])
        self.assertEqual(r["candidate_ID_recovery_by_scenario"]["new_e_all"]["unexpected"],
                         1039)
        self.assertEqual(r["stage_counts_by_scenario"]["e_use_phot"],[446055,30176,119])

    def test_reject_broken_id_identity_even_if_counts_match(self):
        d=copy.deepcopy(self.data)
        d["scenarios"]["new_e_use_phot"]["not_in_frozen_125_count"]=1
        with self.assertRaises(ValueError):
            summarize(d)

    def test_reject_bad_stage_monotonicity(self):
        d=copy.deepcopy(self.data)
        d["scenarios"]["new_e_all"]["stage_counts"][1]=5000
        with self.assertRaises(ValueError):
            summarize(d)


if __name__=="__main__":
    unittest.main()
