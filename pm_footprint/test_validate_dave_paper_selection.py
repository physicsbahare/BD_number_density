"""Offline regression tests for manuscript's signed dropouts and hybrid F277W denominator."""
import unittest
from pm_footprint.validate_dave_paper_selection import evaluate, analyze


def row(source_id=1, f277=2.0, e277=2.0, f444=150.0, f115=0., f150=0.,
        mag=25.96):
    r={"candidate_id":str(source_id),"crossmatch_separation_arcsec":"0",
       "use_phot":"1","use_aper":"0.32","FLAG_AUTO_LW":"0",
       "FLAG_BAD_KRON_LW":"0","mag_F444W":str(mag)}
    for band, flux, sigma in (("F115W",f115,2.),
                              ("F150W",f150,2.),
                              ("F277W",f277,e277),
                              ("F444W",f444,5.)):
        r["f_"+band]=str(flux)
        r["new_e_"+band]=str(sigma)
        r["e_"+band]=str(sigma/2)
    return r


class TestDaveSelection(unittest.TestCase):
    def test_positive_f277_does_not_use_two_sigma(self):
        a=evaluate([row(f277=2,e277=5,f444=90)])
        self.assertEqual(a["all_cuts_passed"],1)  # ratio 45, not 90/10=9

    def test_negative_f277_uses_two_sigma(self):
        a=evaluate([row(f277=-10,e277=2,f444=150)])
        self.assertEqual(a["all_cuts_passed"],1)  # 150/(2*2)=37.5

    def test_positive_f277_high_ratio_can_fail(self):
        a=evaluate([row(f277=4,e277=5,f444=90)])
        self.assertEqual(a["all_cuts_passed"],0)  # 22.5

    def test_blue_snr_is_signed_not_absolute(self):
        a=evaluate([row(f277=2,e277=2,f444=150,f115=-20)])
        self.assertEqual(a["all_cuts_passed"],1)
        self.assertEqual(a["negative_blue_SNR_at_most_minus2"]["F115W"],[1])

    def test_corrected_vs_raw_snr(self):
        a=analyze([row(source_id=i, f277=2,e277=2,f444=150) for i in range(1,126)])
        self.assertEqual(a["corrected_new_e"]["all_cuts_passed"],125)
        self.assertEqual(a["uncorrected_e"]["all_cuts_passed"],0)  # raw F277 SNR=2, fails strict <2
        b=evaluate([row(f277=0.4,e277=2,f444=28)])
        self.assertEqual(b["all_cuts_passed"],1)  # corrected SNR=5.6

    def test_invalid_sky_match_rejected(self):
        rows=[row(source_id=i) for i in range(1,126)]
        rows[-1]["crossmatch_separation_arcsec"]="1.0"
        with self.assertRaises(ValueError):
            analyze(rows)


if __name__=="__main__":
    unittest.main()
