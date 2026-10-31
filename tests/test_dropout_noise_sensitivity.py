"""Unit tests for conditional dropout noise diagnostic."""
import math
from pm_footprint.dropout_noise_sensitivity import evaluate, ab_ujy

def test_color_depth_conversion():
    r=evaluate()
    assert abs(r["f277_2sigma_limit_needed_ab"]-29.8928)<0.0001
    assert abs(r["f277_5sigma_depth_needed_ab"]-28.89795)<0.0001
    assert abs(r["f277_measured_upper_threshold_ujy"]-ab_ujy(26.2)/30)<1e-12

def test_noise_selection_policies():
    r=evaluate()
    for depth in r["depth_scenarios"]:
        z,r100,r30,r10=depth["true_ratio_scenarios"]
        assert abs(r30["p_measured_f277_below_f444_over30_allowing_negative"]-0.5)<1e-12
        assert 0<=z["p_positive_measured_f277_and_ratio_gt30_and_f277_snr_lt2"]<=0.5
        assert r10["p_positive_measured_f277_and_ratio_gt30_and_f277_snr_lt2"]<r100["p_positive_measured_f277_and_ratio_gt30_and_f277_snr_lt2"]
        for row in (z,r100,r30,r10):
            assert 0<=row["p_measured_f277_below_f444_over30_allowing_negative"]<=1
            assert 0<=row["p_positive_measured_f277_and_ratio_gt30_and_f277_snr_lt2"]<=row["p_measured_f277_below_f444_over30_allowing_negative"]

def test_analytic_gaussian_limits():
    d=evaluate(m5_depths=(27.71,))["depth_scenarios"][0]
    z=d["true_ratio_scenarios"][0]
    assert abs(z["p_measured_f277_below_f444_over30_allowing_negative"]-z["p_positive_measured_f277_and_ratio_gt30_and_f277_snr_lt2"]-0.5)<1e-12
