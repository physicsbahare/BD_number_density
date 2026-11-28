import numpy as np
from src.ali_spatial_counts import integrate_cells


def test_default_boundary_and_flux_inverse_square_scaling():
    noise=np.ones((1,4));ref=noise.copy();ref[0,-1]=1e-8
    args=([0.,0.,0.,10.],.001,noise,ref,np.ones((1,4)),np.array([1e-5]))
    old,_=integrate_cells(*args)
    explicit,_=integrate_cells(*args,depth_sigma=5.)
    # Increasing the threshold by four shrinks reach by two: equivalently,
    # reduce only the intrinsic F444W d25 by two in the near-perfect-detection limit.
    shorter,_=integrate_cells(*args,depth_sigma=20.)
    equivalent,_=integrate_cells([0.,0.,0.,5.],*args[1:],depth_sigma=5.)
    np.testing.assert_allclose(old,explicit,rtol=1e-10)
    np.testing.assert_allclose(shorter,equivalent,rtol=1e-6)
    assert shorter.sum()<old.sum()
    scaled,_=integrate_cells(*args,boundary_flux_scale=np.array([4.]))
    np.testing.assert_allclose(shorter,scaled,rtol=1e-10)
