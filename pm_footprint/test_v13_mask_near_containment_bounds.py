"""Tests for conditional JWST-mask/WHT near-containment set bounds."""
import unittest
from pm_footprint.v13_mask_near_containment_bounds import bound_phase


class TestNearContainment(unittest.TestCase):
    def test_full_containment(self):
        v = bound_phase(.6, .6, .6, .4, .03, .02, .03)
        self.assertAlmostEqual(v["conditional_corrected_red_fraction"][0], .38/.57)
        self.assertAlmostEqual(v["conditional_corrected_red_fraction"][1], .38/.57)

    def test_grid_of_feasible_assignments(self):
        A, W, I, R = .6, .58, .55, .3
        Ma, Mr, Mu = .03, .02, .035
        v = bound_phase(A, W, I, R, Ma, Mr, Mu)
        lo, hi = v["conditional_corrected_red_fraction"]
        for ia in range(11):
            a = min(A-I, Mr)*ia/10
            for ib in range(11):
                b = min(A-I-a, Ma-Mr)*ib/10
                for ic in range(11):
                    c = min(W-I, Mu-Ma)*ic/10
                    f = (R-Mr+a)/(W-Ma+a+b-c)
                    self.assertGreaterEqual(f+1e-12, lo)
                    self.assertLessEqual(f-1e-12, hi)

    def test_real_saved_v13_phase(self):
        v = bound_phase(.5338989551843558, .533609115207066,
                        .53350011666025, .31740570799788875,
                        .013005827189953518, .007515930285651531,
                        .013071183105480923)
        self.assertAlmostEqual(v["conditional_corrected_red_fraction"][0], .5947956100860958)
        self.assertAlmostEqual(v["conditional_corrected_red_fraction"][1], .5956358501314067)

    def test_bad_nesting_rejected(self):
        with self.assertRaises(ValueError):
            bound_phase(.5, .5, .5, .3, .04, .05, .06)
        with self.assertRaises(ValueError):
            bound_phase(.5, .5, .5, .51, .04, .02, .06)


if __name__ == "__main__":
    unittest.main()
