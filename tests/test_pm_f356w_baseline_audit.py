"""Unit tests for F444W-to-F356W temporal opportunity audit."""
import numpy as np
import pandas as pd
from pm_footprint.pm_f356w_baseline_audit import (
    best_interval_pair, candidate_report, clean_observations,
)


def row(obsid, band, start, end, proposal="1727"):
    return {"obsid": obsid, "filters": band, "t_min": start, "t_max": end,
            "proposal_id": proposal}


def test_disjoint_interval_gap_is_not_midpoint_gap():
    best = best_interval_pair(
        [row("a", "F444W", 100, 102), row("b", "F356W", 500, 510)], "F356W")
    assert best["gap_days"] == 398
    assert best["midpoint_days"] == 404


def test_overlapping_intervals_do_not_count_as_independent_epochs():
    assert best_interval_pair(
        [row("a", "F444W", 100, 110), row("b", "F356W", 105, 120)], "F356W") is None


def test_same_obsid_cannot_be_paired_with_itself():
    assert best_interval_pair(
        [row("a", "F444W", 100, 102), row("a", "F356W", 500, 501)], "F356W") is None


def test_select_longest_gap_and_retain_provenance():
    records = [row("a", "F444W", 100, 101),
               row("b", "F356W", 400, 401, "5893"),
               row("c", "F356W", 800, 801, "6434")]
    r = candidate_report(387132, records)
    assert r["n_primary_cosmosweb_f444w"] == 1
    assert r["best"]["F356W"]["gap_days"] == 699
    assert r["best"]["F356W"]["other_proposal"] == "6434"
    assert r["best"]["F277W"] is None


def test_filter_public_and_valid_times():
    def obs(visibility, t0, t1):
        return {"obs_collection": "JWST", "instrument_name": "NIRCAM/IMAGE",
                "dataproduct_type": "image", "intentType": "science",
                "dataRights": visibility, "t_min": t0, "t_max": t1}
    raw = pd.DataFrame([obs("PUBLIC", 1, 2), obs("PROPRIETARY", 1, 2),
                        obs("PUBLIC", 2, 1), obs("PUBLIC", np.nan, 1)])
    assert len(clean_observations(raw)) == 1
