"""Offline set-theoretic regression tests: no WHT, FITS or MAST access."""
import pytest

from pm_footprint.wht_redder_bounds import intersection_bounds, analyze


def test_full_union_exact():
    assert intersection_bounds(10, 10, 6) == (6, 6)


def test_partial_subset_bounds():
    assert intersection_bounds(10, 8, 6) == (4, 6)
    assert intersection_bounds(10, 2, 6) == (0, 2)


def test_invalid_subset():
    with pytest.raises(ValueError):
        intersection_bounds(5, 6, 3)
    with pytest.raises(ValueError):
        intersection_bounds(5, 3, 6)


def test_cached_summary_and_candidate_crosscheck_synthetic():
    ids = list(range(1, 126))
    wht_rows = {
        i: {"wht_0_fraction": "1", "wht_0.1_fraction": "1",
            "wht_0.5_fraction": "1",
            **{f"in_f444w_plus_red_{gap}d": str(i <= 62)
               for gap in (30, 180, 365, 730)}} for i in ids
    }
    filter_rows = {
        i: {
            **{f"in_red_union_{gap}d": str(i <= 62)
               for gap in (30, 180, 365, 730)},
            **{f"in_red_without_f277w_{gap}d": str(i <= 62)
               for gap in (30, 180, 365, 730)},
            **{f"in_{f}_{gap}d": str(i <= 62 if f == "f356w" else False)
               for gap in (30, 180, 365, 730)
               for f in ("f277w", "f356w", "f410m", "f444w")},
        } for i in ids
    }
    wht = {
        "area_deg2": {
            t: {"archive_f444w": 0.53,
                **{f"f444w_plus_red_{gap}d": 0.30
                   for gap in (30, 180, 365, 730)}}
            for t in ("0", "0.1", "0.5")
        },
        "candidate_counts": {
            t: {f"f444w_plus_red_{gap}d": 62
                for gap in (30, 180, 365, 730)}
            for t in ("0", "0.1", "0.5")
        },
    }
    filt = {"baselines": {
        str(gap): {
            "area_deg2": {
                "red_union": 0.32, "red_without_f277w": 0.31,
                "F277W": 0.02, "F356W": 0.30,
                "F410M": 0.04, "F444W": 0.01,
            },
            "candidate_ids": {"red_union": ids[:62]},
        } for gap in (30, 180, 365, 730)
    }}
    result = analyze(wht, filt, wht_rows, filter_rows)
    redder = result["0"]["365"]["filter_intersection_bounds"]["red_without_f277w"]
    assert redder["wht_intersection_area_bounds_deg2"] == pytest.approx([0.29, 0.30])
    assert redder["candidate_count_from_existing_memberships"] == 62
    assert result["0"]["365"]["conditional_external_usable_area"]["redder_fraction_bounds"][0] > 0


def test_mismatched_roster_fails():
    with pytest.raises(ValueError, match="rosters differ"):
        analyze({}, {}, {1: {}}, {2: {}})
