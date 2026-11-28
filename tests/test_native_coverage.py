import numpy as np
from scipy.ndimage import binary_erosion
from src.depth_reconstruction import aperture_kernel
from src.ali_native_coverage import aperture_validity


def test_fast_native_aperture_mask_matches_exact_fractional_support():
    rng=np.random.default_rng(42);mask=rng.random((123,157))>.002
    expected=binary_erosion(mask,structure=aperture_kernel(.1,.03)>0,border_value=0)
    np.testing.assert_array_equal(aperture_validity(mask),expected)
