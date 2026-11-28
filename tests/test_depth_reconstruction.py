import gzip
import numpy as np
from astropy.io import fits
from src.depth_reconstruction import aperture_kernel, streamed_samples


def test_fractional_aperture_and_signed_flux():
    k = aperture_kernel(0.1, 0.03)
    assert np.isclose(k.sum(), np.pi*(0.1/0.03)**2)
    assert np.isclose((np.full(k.shape, -2.)*k).sum(), -2*k.sum())


def test_streamed_samples_preserve_rows_and_cutouts(tmp_path):
    a = np.arange(35*41, dtype='float32').reshape(35, 41)
    p = tmp_path/'a.fits.gz'
    fits.PrimaryHDU(a).writeto(p)
    xy = np.array([[9, 10], [27, 22]])
    h, centers, stamps = streamed_samples(p, xy, np.array([[15, 17]]), 3)
    np.testing.assert_array_equal(centers, a[xy[:,1], xy[:,0]])
    np.testing.assert_array_equal(stamps[0], a[14:21,12:19])
    assert h['NAXIS1'] == 41


def test_unsigned_segmentation_bzero(tmp_path):
    a = np.arange(35*41,dtype='uint32').reshape(35,41)
    p = tmp_path/'seg.fits.gz'
    fits.PrimaryHDU(a).writeto(p)
    _, centers, stamps = streamed_samples(p,np.array([[9,10]]),np.array([[15,17]]),3)
    assert centers[0] == a[10,9]
    np.testing.assert_array_equal(stamps[0],a[14:21,12:19])
