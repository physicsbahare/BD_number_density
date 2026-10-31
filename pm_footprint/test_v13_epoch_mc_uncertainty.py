"""Frozen-input and mathematical tests for mask MC uncertainty."""
import json
import unittest
from pathlib import Path
from pm_footprint.v13_epoch_mc_uncertainty import analyze, wilson

BASE = Path("pm_footprint/results")


def inputs():
    return [json.loads((BASE/p).read_text()) for p in (
        "local_v13_epoch_overlay_20261031.json",
        "jwst_starmask/epoch_mc_20261031.json",
        "jwst_starmask/summary.json",
        "jwst_starmask/union_area_grid_qc.json")]


class TestMaskMCUncertainty(unittest.TestCase):
    def test_wilson_boundaries(self):
        self.assertEqual(wilson(0,100,.95)[0],0.)
        self.assertEqual(wilson(100,100,.95)[1],1.)

    def test_invalid_counts(self):
        for k,n in [(-1,10),(11,10),(0,0),(1,2.5)]:
            with self.assertRaises(ValueError):
                wilson(k,n,.95)
        with self.assertRaises(ValueError):
            wilson(5,10,1.)

    def test_nested_counts(self):
        src=inputs()
        src[1]["sample_membership_counts"]["red365"]=1001
        with self.assertRaises(ValueError):
            analyze(*src)

    def test_saved_results(self):
        x=analyze(*inputs())
        lo,hi=x["mc_fraction_envelope_all_phases"]
        self.assertAlmostEqual(lo,0.5938923434506636,places=10)
        self.assertAlmostEqual(hi,0.5963721866678143,places=10)
        self.assertEqual(x["n_mask_draws"],1000)
        self.assertLess(lo,x["phase_results"][0]["nominal_conditional_fraction"][0])
        self.assertGreater(hi,x["phase_results"][0]["nominal_conditional_fraction"][1])
        self.assertLess(x["placement_cap_fraction_envelope_all_phases"][0],lo)
        self.assertGreater(x["placement_cap_fraction_envelope_all_phases"][1],hi)

    def test_maskless_control(self):
        src=inputs()
        src[1]["sample_membership_counts"]["red365"]=0
        src[1]["sample_membership_counts"]["survey"]=0
        x=analyze(*src)
        self.assertTrue(0<x["phase_results"][0]["unmasked_fraction"]<1)


if __name__=="__main__":
    unittest.main()
