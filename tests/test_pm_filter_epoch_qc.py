"""Offline exact-filter audit tests; no MAST or WHT network access."""
import pandas as pd
import pytest

from pm_footprint.filter_epoch_qc import filter_epoch_masks


def rect(x1, y1, x2, y2):
    return (f"POLYGON ICRS {x1} {y1} {x2} {y1} "
            f"{x2} {y2} {x1} {y2}")


def row(start, band, footprint, pid=1727):
    return dict(obs_collection="JWST", instrument_name="NIRCAM/IMAGE",
                dataproduct_type="image", intentType="science", dataRights="PUBLIC",
                proposal_id=pid, filters=band, t_min=start, t_max=start + 1,
                s_region=footprint)


def roster():
    ids = [387132, 1651050] + list(range(2000000, 2000123))
    return pd.DataFrame({
        "candidate_id": ids,
        "ra_deg": [150.04, 150.18] + [150.04] * 123,
        "dec_deg": [2.04] * 125,
    })


def test_filter_specific_red_opportunities_and_epoch_baseline():
    base = rect(149.9, 1.9, 150.3, 2.2)
    left = rect(150.0, 2.0, 150.10, 2.10)
    right = rect(150.15, 2.0, 150.23, 2.10)
    obs = pd.DataFrame([
        row(60000, "F444W", base),
        row(60400, "F277W", left, 2000),
        row(60400, "F356W", right, 2001),
        row(60010, "F410M", base, 2002),  # too close in time
        row(60400, "F444W", left, 2003),
        row(60400, "F356W", left, 2004, )  # separate redder coverage
    ])
    result, membership, area = filter_epoch_masks(obs, roster(), baselines=(365, 730))
    assert 0 < area < 1
    a = result["365"]["area_deg2"]
    n = result["365"]["candidate_counts"]
    assert a["F277W"] > 0 and a["F356W"] > 0 and a["F444W"] > 0
    assert a["F410M"] == 0
    assert a["red_without_f277w"] <= a["red_union"]
    assert n["red_union"] == 125
    assert n["red_without_f277w"] == 125
    assert n["F277W"] == 124
    assert n["F356W"] == 125
    assert result["730"]["area_deg2"]["red_union"] == 0
    assert membership["in_red_union_365d"].all()


def test_blue_only_overlap_not_redder_detectability():
    obs = pd.DataFrame([
        row(60000, "F444W", rect(149.9, 1.9, 150.3, 2.2)),
        row(60400, "F277W", rect(149.9, 1.9, 150.3, 2.2), 2000),
    ])
    result, _, _ = filter_epoch_masks(obs, roster(), baselines=(365,))
    n = result["365"]["candidate_counts"]
    assert n["F277W"] == 125
    assert n["red_without_f277w"] == 0
    assert n["red_union"] == 125


def test_exclusive_access_does_not_count():
    obs = pd.DataFrame([
        row(60000, "F444W", rect(149.9, 1.9, 150.3, 2.2)),
        row(60400, "F356W", rect(149.9, 1.9, 150.3, 2.2), 2000),
    ])
    obs.loc[1, "dataRights"] = "EXCLUSIVE_ACCESS"
    result, _, _ = filter_epoch_masks(obs, roster(), baselines=(365,))
    assert result["365"]["candidate_counts"]["red_union"] == 0


def test_roster_guard():
    with pytest.raises(ValueError, match="125"):
        filter_epoch_masks(pd.DataFrame(), roster().iloc[:124])
