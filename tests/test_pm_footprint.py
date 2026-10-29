import pandas as pd
import pytest

pytest.importorskip("mocpy", reason="Footprint geometry is tested in its dedicated workflow")
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
    x = overlap([(60000, 60000.1, region), (60001, 60001.1, region)], depth=11, gap=30)
    assert x.empty()
    assert deg2(region) > 0


def test_filters_parser_and_incomplete_roster():
    assert filters("F444W;F277W") == {"F444W", "F277W"}
    with pytest.raises(ValueError, match="125"):
        calculate(pd.DataFrame([row(60000, "F444W", rectangle(149.8, 1.9, 150.4, 2.4))]), sample().iloc[:5])

def test_mast_multi_polygon_is_union_not_bridged():
    a = rectangle(150.00, 2.00, 150.05, 2.05)
    b = rectangle(150.15, 2.15, 150.20, 2.20)
    m = polygon(a + " " + b, depth=13)
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    pts = SkyCoord([150.025, 150.175, 150.100] * u.deg,
                    [2.025, 2.175, 2.100] * u.deg)
    inside = m.contains_skycoords(pts)
    assert inside.tolist() == [True, True, False]


def test_grism_not_counted_as_imaging_epoch():
    base = rectangle(149.8, 1.9, 150.4, 2.4)
    raw = pd.DataFrame([
        row(60000, "F444W", base),
        dict(row(60100, "F444W;GRISMR", base, program="5398"),
             instrument_name="NIRCAM/GRISM"),
    ])
    summary, _, _ = calculate(raw, sample(), depth=11)
    assert summary["area_deg2"]["same_f444w"] == 0
    assert summary["converted_footprints"] == 1

def test_epoch_intervals_prevent_false_30_day_baseline():
    region = polygon(rectangle(149.8, 1.9, 150.4, 2.4), depth=11)
    # 34 days between t_min values but only 28 days between intervals.
    no = overlap([(60000, 60006, region), (60034, 60040, region)],
                 depth=11, gap=30)
    yes = overlap([(60000, 60006, region), (60036, 60042, region)],
                  depth=11, gap=30)
    assert no.empty()
    assert not yes.empty()


def test_multiday_second_epoch_is_retained():
    base = rectangle(149.8, 1.9, 150.4, 2.4)
    second = row(60100, "F356W", base, program="5893")
    second["t_max"] = 60106
    raw = pd.DataFrame([row(60000, "F444W", base), second])
    summary, _, _ = calculate(raw, sample(), depth=11)
    assert summary["multi_day_observation_rows_retained"] == 1
    assert summary["area_deg2"]["f444w_plus_red"] > 0
    assert summary["area_deg2"]["same_f444w"] == 0
