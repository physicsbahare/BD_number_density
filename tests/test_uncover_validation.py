import pandas as pd
import numpy as np

from src.validate_uncover import (
    PRIMARY_TEFF_BINS,
    field_count,
    has_primary_f277w_reach,
)


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


def _complete_table():
    return pd.DataFrame(
        {
            "teff_bin_K": list(PRIMARY_TEFF_BINS),
            "d25_F115W_pc": np.linspace(100.0, 700.0, len(PRIMARY_TEFF_BINS)),
            "d25_F444W_pc": np.linspace(200.0, 800.0, len(PRIMARY_TEFF_BINS)),
            "d25_F277W_pc": np.linspace(150.0, 750.0, len(PRIMARY_TEFF_BINS)),
        }
    )


def test_primary_validation_requires_all_f277w_bins():
    table = _complete_table()
    assert has_primary_f277w_reach(table)

    missing = table.iloc[:-1].copy()
    assert not has_primary_f277w_reach(missing)

    nan_value = table.copy()
    nan_value.loc[0, "d25_F277W_pc"] = np.nan
    assert not has_primary_f277w_reach(nan_value)

    zero_value = table.copy()
    zero_value.loc[0, "d25_F277W_pc"] = 0.0
    assert not has_primary_f277w_reach(zero_value)

    negative_value = table.copy()
    negative_value.loc[0, "d25_F277W_pc"] = -1.0
    assert not has_primary_f277w_reach(negative_value)


def test_primary_validation_rejects_duplicate_required_bin():
    table = _complete_table()
    duplicate = pd.concat([table, table.iloc[[0]]], ignore_index=True)
    assert not has_primary_f277w_reach(duplicate)


def test_primary_validation_rejects_missing_f277w_column():
    table = _complete_table().drop(columns=["d25_F277W_pc"])
    assert not has_primary_f277w_reach(table)
