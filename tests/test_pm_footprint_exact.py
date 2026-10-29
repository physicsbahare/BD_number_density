import math

import pandas as pd
import pytest

pytest.importorskip("shapely", reason="Exact geometry is tested in dedicated footprint CI")
from pm_footprint.exact_geometry import compute, polygon_region, DEG2_PER_SR


def rect(ra1, dec1, ra2, dec2):
    return f"POLYGON {ra1} {dec1} {ra2} {dec1} {ra2} {dec2} {ra1} {dec2}"


def roster():
    return pd.DataFrame({
        "candidate_id": [387132, 1651050] + list(range(2000000, 2000123)),
        "ra_deg": [150.04] * 125,
        "dec_deg": [2.04] * 125,
    })


def row(t1, t2, band, region, pid=1727):
    return dict(
        obs_collection="JWST", instrument_name="NIRCAM/IMAGE",
        dataproduct_type="image", intentType="science", dataRights="PUBLIC",
        proposal_id=pid, filters=band,
        t_min=t1, t_max=t2, s_region=region,
    )


def test_equal_area_projection_and_disjoint_polygons():
    small = rect(150, 2, 150.1, 2.1)
    b = rect(150.3, 2.3, 150.4, 2.4)
    g = polygon_region(small + " " + b)
    assert 0.0198 < g.area * DEG2_PER_SR < 0.0201


def test_exact_overlap_respects_interval_baselines():
    base = rect(149.8, 1.9, 150.4, 2.4)
    # A six-day aggregate red observation is a genuine second epoch at +100d.
    df = pd.DataFrame([
        row(60000, 60000.1, "F444W", base),
        row(60100, 60106, "F356W", base, pid=5893),
    ])
    summary, members = compute(df, roster(), baselines=(30, 180))
    assert 0.29 < summary["archive_f444w_area_deg2"] < 0.31
    assert summary["baseline_results"]["30"]["candidate_counts"]["f444w_plus_red"] == 125
    assert summary["baseline_results"]["180"]["candidate_counts"]["f444w_plus_red"] == 0
    assert summary["baseline_results"]["30"]["candidate_counts"]["same_f444w"] == 0
    assert members["in_f444w_plus_red_30d"].all()
