import numpy as np
import pytest
from pm_footprint.jwst_mask_epoch_from_sample import epoch_flags, wht_fraction_bounds


def test_epoch_gap_is_interval_separation_not_start_to_start():
    # 100--101 and 466--467 are separated by exactly 365 days.
    survey, flags = epoch_flags([100, 466], [101, 467],
                                [True, False], [True, False])
    assert survey
    assert flags[365]["red"]
    assert not flags[365]["same_f444w"]
    assert not flags[730]["red"]
    # One day less does not pass.
    _, flags2 = epoch_flags([100, 465], [101, 466],
                            [True, False], [True, False])
    assert not flags2[365]["red"]


def test_repeated_f444w_counts_as_red_and_same():
    survey, flags = epoch_flags([100, 466], [101, 467],
                                [True, True], [True, False])
    assert survey and flags[365]["red"] and flags[365]["same_f444w"]


def test_no_primary_survey_f444w_does_not_count():
    survey, flags = epoch_flags([100, 466], [101, 467],
                                [True, False], [False, False])
    assert not survey and not flags[365]["red"]


def test_conditional_bounds_include_every_feasible_discrete_allocation():
    # A synthetic 100-cell survey: 60 red, 40 not red. WHT keeps 50 red
    # and 30 nonred. Star mask contains 12 red + 8 nonred cells.
    S,R,W,WR,b,c = 100.,60.,80.,50.,20.,12.
    out = wht_fraction_bounds(S,R,W,WR,b,c)
    assert out["fraction_lower"] <= out["fraction_upper"]
    for x in range(3,13): # mask in WHT red: 12 minus at most 10 missing red
        for y in range(0,9): # mask in WHT nonred: 8 minus at most 10 missing
            if x>WR or y>W-WR: continue
            frac=(WR-x)/(W-x-y)
            assert out["fraction_lower"]-1e-12<=frac<=out["fraction_upper"]+1e-12


def test_actual_cosmos_bounds_stay_in_plausible_range():
    out=wht_fraction_bounds(0.5338989551843558,0.31748660889853736,
        0.528823684998808,0.31462445722188803,
        0.013005827189953518,0.007515930285651531)
    assert 0.592 < out["fraction_lower"] < 0.594
    assert 0.597 < out["fraction_upper"] < 0.599


def test_inconsistent_areas_rejected():
    with pytest.raises(ValueError):
        wht_fraction_bounds(1.,0.4,0.8,0.6,0.1,0.05)
