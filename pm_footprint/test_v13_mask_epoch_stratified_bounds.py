"""Regression and synthetic feasibility tests for epoch-stratified bounds."""
import unittest
from pm_footprint.v13_mask_epoch_stratified_bounds import phase_bounds


class TestEpochStratifiedBounds(unittest.TestCase):
    def test_full_overlap(self):
        x = phase_bounds(.6, .4, .6, .6, .4, .03, .02, .03)
        self.assertAlmostEqual(x["conditional_corrected_red_fraction"][0], .38/.57)
        self.assertAlmostEqual(x["conditional_corrected_red_fraction"][1], .38/.57)

    def test_saved_phase_nominal_cap(self):
        x = phase_bounds(
            .5338989551843558, .31748660889853736,
            .533609115207066, .53350011666025, .31740570799788875,
            .013005827189953518, .007515930285651531, .0141992195)
        self.assertAlmostEqual(x["archive_red_not_in_wht_deg2"], .00008090090064860611)
        self.assertAlmostEqual(x["archive_nonred_not_in_wht_deg2"], .00031793762345722065)
        self.assertAlmostEqual(x["conditional_corrected_red_fraction"][0], .5948879839387807)
        self.assertAlmostEqual(x["conditional_corrected_red_fraction"][1], .5954388229527625)

    def test_enumerated_feasible_allocations(self):
        S, R, W, WA, WR = .6, .35, .58, .57, .34
        MA, MR, MU = .03, .02, .035
        result = phase_bounds(S,R,W,WA,WR,MA,MR,MU)
        lo,hi = result["conditional_corrected_red_fraction"]
        amax = result["max_mask_red_outside_wht_deg2"]
        bmax = result["max_mask_nonred_outside_wht_deg2"]
        cmax = result["max_mask_outside_archive_in_wht_deg2"]
        for i in range(11):
            a = amax*i/10
            for j in range(11):
                b = bmax*j/10
                for k in range(11):
                    c = cmax*k/10
                    f=(WR-MR+a)/(W-MA+a+b-c)
                    self.assertLessEqual(lo-1e-12,f)
                    self.assertLessEqual(f,hi+1e-12)

    def test_red_cap_is_not_total_missing(self):
        x = phase_bounds(.6,.35,.58,.57,.34,.03,.02,.035)
        self.assertAlmostEqual(x["max_mask_red_outside_wht_deg2"], .01)
        self.assertAlmostEqual(x["max_mask_nonred_outside_wht_deg2"], .02)

    def test_bad_nesting(self):
        with self.assertRaises(ValueError):
            phase_bounds(.6,.35,.58,.57,.36,.03,.02,.035)
        with self.assertRaises(ValueError):
            phase_bounds(.6,.35,.58,.57,.34,.03,.04,.035)
        with self.assertRaises(ValueError):
            phase_bounds(.6,.35,.58,.57,.34,.03,.02,.025)


if __name__ == "__main__":
    unittest.main()
