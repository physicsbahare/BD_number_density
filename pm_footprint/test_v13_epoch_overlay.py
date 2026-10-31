"""Offline synthetic WCS and polygon-overlay tests; do not open original WHT."""
import unittest
import numpy as np
from shapely.geometry import box, Polygon
from pm_footprint.exact_geometry import lambert
from pm_footprint.v13_epoch_overlay import cell_centers, reconstruct_wcs, evaluate_phase


class TestV13EpochOverlay(unittest.TestCase):
    def test_partial_edge_cells_and_sample_phases(self):
        x,y,wx,wy=cell_centers(3,8,7,0.25)
        np.testing.assert_array_equal(x,[0,3,6])
        np.testing.assert_array_equal(y,[0,3,6])
        np.testing.assert_array_equal(wx,[3,3,2])
        np.testing.assert_array_equal(wy,[3,3,1])
        x2,y2,_,_=cell_centers(3,8,7,0.75)
        np.testing.assert_array_equal(x2,[2,5,7])
        np.testing.assert_array_equal(y2,[2,5,6])

    def test_full_and_empty_exact_geometry(self):
        cd=np.diag([-0.03/3600,0.03/3600])
        g={"stride":np.array(4), "ncols":np.array(8), "nrows":np.array(8),
           "wcs_cd":cd, "wcs_crpix":np.array([4.,4.]),
           "wcs_crval":np.array([150.12,2.20]),
           "wht_positive_25":np.ones((2,2),bool)}
        w=reconstruct_wcs(g)
        ra,dec=w.all_pix2world(np.array([0,7]),np.array([0,7]),0)
        xx,yy=lambert(ra,dec)
        arch=box(min(xx)-1e-4,min(yy)-1e-4,max(xx)+1e-4,max(yy)+1e-4)
        result=evaluate_phase(g,.25,{"archive":arch,"red_365d":arch,
                                     "same_f444w_365d":Polygon()},chunk_size=3)
        self.assertEqual(result["wht_positive_centers_outside_archive"],0)
        self.assertEqual(result["wht_positive_cells"],4)
        self.assertAlmostEqual(result["area_intersections_deg2"]["red_365d"],
                               result["wht_positive_area_deg2"],places=14)
        self.assertEqual(result["area_intersections_deg2"]["same_f444w_365d"],0.)
        expected=64*(0.03/3600)**2
        self.assertAlmostEqual(result["wht_positive_area_deg2"],expected,
                               delta=expected*1e-5)


if __name__=="__main__":
    unittest.main()
