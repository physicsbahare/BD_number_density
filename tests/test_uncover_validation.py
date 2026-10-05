import pandas as pd
import numpy as np

from src.validate_uncover import field_count, has_primary_f277w_reach


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


def test_primary_validation_requires_f277w_reach():
    current = pd.DataFrame(
        {
            "teff_bin_K": ["450-600"],
            "d25_F115W_pc": [120.0],
            "d25_F444W_pc": [306.57],
        }
    )
    assert not has_primary_f277w_reach(current)

    ready = current.copy()
    ready["d25_F277W_pc"] = [200.0]
    assert has_primary_f277w_reach(ready)
