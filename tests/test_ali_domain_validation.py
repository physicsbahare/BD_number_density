import numpy as np
from src.ali_domain_validation import zero_source_probability, halo_tail_power


def test_zero_flux_native_selection_is_not_zero():
    noise=np.ones((1,5));den=np.ones((1,5));den[0,4]=.43
    p=zero_source_probability(noise,den,np.zeros_like(noise))
    assert 0<p[0]<.05


def test_existing_halo_has_divergent_number_tail():
    assert np.isclose(halo_tail_power(2.77),.23)
    assert halo_tail_power(3.5)<0


def test_native_coverage_not_detection_positivity():
    from src.ali_detection_footprint_validation import packed_count_cells
    valid=np.ones((4,8),dtype=bool)
    detection=np.zeros_like(valid);detection[:,::2]=True
    full=packed_count_cells(np.packbits(valid,axis=1),8,4)
    positive=packed_count_cells(np.packbits(valid&detection,axis=1),8,4)
    assert np.array_equal(full,[16,16])
    assert np.array_equal(positive,[8,8])
