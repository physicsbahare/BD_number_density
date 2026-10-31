"""Tests for worst-case JWST-mask geometric coverage fraction."""
import numpy as np
import pytest
from pm_footprint.mask_robust_bounds import bounds

def test_analytic_extremes_and_random_feasible_allocations():
    W, WR, M = .528823684998808, .31462445722188803, .014199219444444447
    b = bounds(W, WR, M)
    assert abs(b["lower"] - .5837756614501792) < 1e-12
    assert abs(b["upper"] - .611367080814878) < 1e-12
    assert b["lower"] < b["uncorrected"] < b["upper"]
    rng = np.random.default_rng(20261031)
    x = rng.random(10000) * M
    y = rng.random(10000) * (M-x)
    f = (WR-x)/(W-x-y)
    assert (f >= b["lower"]-1e-12).all()
    assert (f <= b["upper"]+1e-12).all()

def test_input_guardrails():
    with pytest.raises(ValueError): bounds(1, .5, .6)
    with pytest.raises(ValueError): bounds(1, 1, .1)
    with pytest.raises(ValueError): bounds(1, .5, -.1)
