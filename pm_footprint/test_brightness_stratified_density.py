"""Offline synthetic tests for magnitude-stratified density."""
import unittest
from pm_footprint.brightness_stratified_density import analyze,rate_ratio,fisher_two_sided
def row(i,mag,red,ra=None):
    return {"candidate_id":i,"mag_F444W":mag,"red365":red,
            "ra_deg":150+i*.01 if ra is None else ra,"dec_deg":2}
class Tests(unittest.TestCase):
    def test_ratio(self):
        self.assertAlmostEqual(rate_ratio(1,1,1,3)["ratio"],2)
    def test_bad_area(self):
        with self.assertRaises(ValueError):rate_ratio(1,1,2,1)
    def test_fisher(self):
        self.assertAlmostEqual(fisher_two_sided([[1,1],[1,1]]),1)
        self.assertAlmostEqual(fisher_two_sided([[0,2],[2,0]]),1/3)
    def test_split(self):
        r=analyze([row(1,22.9,1),row(2,24.5,0),row(3,25,1),row(4,26.1,0)])
        self.assertEqual(r["bright_m444_lt25"]["red"],1)
        self.assertEqual(r["faint_m444_ge25"]["nonred"],1)
    def test_duplicate(self):
        with self.assertRaises(ValueError):analyze([row(1,24,1),row(1,25,0)])
    def test_groups(self):
        r=analyze([row(1,24,1,150),row(2,24.4,1,150+.1/3600),row(3,25,0,151)])
        self.assertEqual(r["proximity_groups_0p5arcsec"]["n_groups"],2)
    def test_bad_red(self):
        with self.assertRaises(ValueError):analyze([row(1,24,2)])
if __name__=="__main__":unittest.main()
