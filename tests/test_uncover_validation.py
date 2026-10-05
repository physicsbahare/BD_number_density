import numpy as np

from src.validate_uncover import field_count, garwood_interval


def test_uncover_count_is_positive():
    n = field_count(
        dmax_pc=700.0,
        area_arcmin2=53.4,
        rho_local_total=4.20e-3,
        ra_deg=3.58125,
        dec_deg=-30.388661,
    )
    assert n > 0
    assert np.isclose(n, 0.214, atol=0.002)


def test_five_count_poisson_interval_contains_baseline_prediction():
    lo, hi = garwood_interval(5, confidence=0.68)
    assert lo < 2.880 < hi
