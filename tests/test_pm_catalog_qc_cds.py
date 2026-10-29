import pandas as pd
import pytest
from pm_footprint.catalog_qc_cds import choose_columns, validate_matches, summarize


def roster():
    return pd.DataFrame({
        "candidate_id": [387132, 1651050],
        "ra_deg": [150.123456, 150.345678],
        "dec_deg": [2.123456, 2.345678],
    })


def cds():
    return pd.DataFrame({
        "Id": [387132, 1651050],
        "RAJ2000": [150.123456, 150.345678],
        "DEJ2000": [2.123456, 2.345678],
        "flag-star": [0, 3],
        "flag-blend": [0, 1],
        "SNR-F444W": [10, 7],
        "SNR-F277W": [1, 0.5],
    })


def test_column_normalization_and_flags():
    table, cols = validate_matches(roster(), cds())
    assert table.cds_coordinate_consistent.all()
    assert table.cds_flagstar.tolist() == [0, 3]
    assert table.cds_flagblend.tolist() == [0, 1]
    report = summarize(table, "test", cols)
    assert report["matched_by_id_and_position_within_0p5arcsec"] == 2
    assert report["cds_flagstar_distribution_valid_position_matches"] == {"0": 1, "3": 1}
    assert report["cds_snrf444w_n_snr_gt_5"] == 2


def test_reject_coordinate_mismatch():
    data = cds()
    data.loc[0, "RAJ2000"] += 0.01
    matched, _ = validate_matches(roster(), data)
    assert matched.cds_id_found.all()
    assert not matched.cds_coordinate_consistent.iloc[0]
    assert matched.cds_coordinate_consistent.iloc[1]


def test_duplicate_ids_rejected():
    data = pd.concat([cds(), cds().iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError, match="duplicate"):
        validate_matches(roster(), data)


def test_required_columns_fail_closed():
    with pytest.raises(ValueError, match="lacks columns"):
        choose_columns(["Id", "RAJ2000", "missing_dec"])
