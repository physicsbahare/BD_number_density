import numpy as np
from src.galactic_density import DEG2_TO_SR, expected_count_component

def test_deg2_to_sr():
    assert np.isclose(DEG2_TO_SR, (np.pi/180.0)**2)

def test_uniform_density_volume():
    rho = 2.0e-3
    d1, d2 = 10.0, 20.0
    area = 1.0
    got, _ = expected_count_component(d1, d2, area, lambda d: rho)
    expected = area * DEG2_TO_SR * rho * (d2**3-d1**3)/3.0
    assert np.isclose(got, expected)
