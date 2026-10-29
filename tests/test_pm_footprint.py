import pandas as pd
import pytest
from pm_footprint.compute import calculate, filters, overlap, polygon, deg2


def rectangle(ra1, dec1, ra2, dec2):
    return (f"POLYGON {ra1} {dec1} {ra2} {dec1} "
            f"{ra2} {dec2} {ra1} {dec2}")


def sample():
    # 125 rows are required, including the two verified movers.
    ids = [387132, 1651050] + list(range(2000000, 2000123))
    return pd.DataFrame({
        "candidate_id": ids,
        "ra_deg": [150.04] * 125,
        "dec_deg": [2.04] * 125,
    })


def row(mjd, band, region, program="1727"):
    return {
        "obs_collection": "JWST",
        "instrument_name": "NIRCAM/IMAGE",
        "intentType": "science",
        "dataproduct_type": "image",
        "proposal_id": program,
        "filters": band,
        "t_min": mjd,
        "t_max": mjd + 0.02,
        "s_region": region,
    }


def test_separate_filter_regimes_are_nested():
    base = rectangle(149.8, 1.9, 150.4, 2.4)
    sub = rectangle(150.0, 2.0, 150.1, 2.1)
    red = rectangle(149.95, 1.95, 150.2, 2.2)
    blue = rectangle(149.90, 1.95, 150.3, 2.3)
    obs = pd.DataFrame([
        row(60000, "F444W", base),
        row(60010, "F444W", sub, program="9999"),  # not independent
        row(60040, "F444W", sub, program="9999"),  # independent F444W
        row(60050, "F356W", red, program="9999"),
        row(60060, "F115W", blue, program="9999"),
    ])
    result, membership, _ = calculate(obs, sample(), depth=11, gap=30)
    a = result["area_deg2"]
    assert 0.1 < a["cweb_f444w"] < 0.5
    assert 0 < a["same_f444w"] <= a["f444w_plus_red"] <= a["any_nircam"]
    assert result["candidates_with_geometric_coverage"]["same_f444w"] == 125
    assert membership["in_cweb_f444w"].all()
    assert result["status"].endswith("REQUIRES_QC")


def test_single_epoch_has_no_overlap():
    region = polygon(rectangle(149.8, 1.9, 150.4, 2.4), depth=11)
    x = overlap([(60000, region), (60001, region)], depth=11, gap=30)
    assert x.empty()
    assert deg2(region) > 0


def test_filters_parser_and_incomplete_roster():
    assert filters("F444W;F277W") == {"F444W", "F277W"}
    with pytest.raises(ValueError, match="125"):
        calculate(pd.DataFrame([row(60000, "F444W", rectangle(149.8, 1.9, 150.4, 2.4))]), sample().iloc[:5])
