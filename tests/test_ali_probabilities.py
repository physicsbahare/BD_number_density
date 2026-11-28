import numpy as np
from scipy.special import ndtr
from src.ali_spatial_counts import signed_probability, density_components
from src.cosmosweb_counts import galactic_lb, disk_rho, halo_rho, local_component_densities
from src.ali_spatial_counts import integrate_cells
from src.cosmosweb_counts import count_to_distance


def test_signed_and_missing_measurements():
    assert np.isclose(signed_probability(0.,1.,1.,2.,False),ndtr(2))
    assert signed_probability(0.,np.nan,1.,2.,False)==0
    assert signed_probability(0.,1.,0.,2.,False)==0
    assert np.isclose(signed_probability(-3.,1.,1.,2.,False),ndtr(5))
    assert np.isclose(signed_probability(6.,1.,1.,5.,True),ndtr(1))


def test_density_components_reuse_existing_normalization():
    r=np.array([10.,100.,1000.]);l,b=galactic_lb();thin,thick,halo=local_component_densities(0.0042)
    expected=np.array([disk_rho(r,l,b,thin,2600,187),disk_rho(r,l,b,thick,3600,900),halo_rho(r,l,b,halo)])
    np.testing.assert_allclose(density_components(r,0.0042),expected)


def test_depth_bound_reproduces_existing_distance_integral():
    f25=10**((31.4-25)/2.5)
    noise=np.array([[1.,1.,1.,f25/(5*10**2)]])
    ref=np.array([[1.,1.,1.,1e-6]])
    area=np.array([(np.pi/180)**2])
    got,_=integrate_cells([0.,0.,0.,1.],.0042,noise,ref,np.ones((1,4)),area,radial_domain='depth')
    expected=count_to_distance(10.,3600.,.0042)*ndtr(2)**3
    assert np.isclose(got.sum(),expected,rtol=1e-6)


def test_f814w_interpolation_retains_donor_provenance_and_excludes_missing_coverage():
    from src.ali_count_run import fill_missing_f814w
    import pytest
    data=np.ones((3,25));data[1:,23]=0;data[1:,3]=np.nan
    xy=np.array([[0.,0.],[100.,0.],[200.,0.]])
    provenance,n=fill_missing_f814w(data,xy,np.array([4,7,8]),np.array([1.,1.,0.]),150.)
    assert n==1 and provenance[1]==4 and np.isfinite(data[1,3])
    assert np.isnan(data[2,3])
    data[1,23]=0;data[1,3]=np.nan
    with pytest.raises(ValueError,match='distance'):
        fill_missing_f814w(data,xy,np.array([4,7,8]),np.array([1.,1.,0.]),50.)


def test_source_noise_zero_limit_and_nonzero_empirical_convolution():
    from scipy.special import ndtr
    d25=np.array([1.,1.,1.,1.,10.])
    noise=np.ones((1,5));ref=noise.copy();ee=noise.copy();area=np.array([1e-5])
    base,_=integrate_cells(d25,.001,noise,ref,ee,area,empirical_cdf=ndtr)
    zero,_=integrate_cells(d25,.001,noise,ref,ee,area,empirical_cdf=lambda z,s:ndtr(z),f814w_source_coefficient=np.array([0.]))
    assert np.allclose(base,zero,rtol=1e-8)
    seen=[]
    def cdf(z,shot):
        seen.append(float(shot[0]));return ndtr(z/np.sqrt(1+shot**2))
    value,_=integrate_cells(d25,.001,noise,ref,ee,area,empirical_cdf=cdf,f814w_source_coefficient=np.array([.1]))
    assert max(seen)>0 and np.isfinite(value).all() and np.all(value>=0)


def test_joint_veto_detection_independent_limit_and_positive_covariance():
    from src.ali_spatial_counts import joint_veto_detection
    from scipy.special import ndtr
    v=joint_veto_detection(2.,5.,1.,1.,1.,1.,0.)
    assert np.isclose(v,.25)
    correlated=joint_veto_detection(2.,5.,1.,1.,1.,1.,.1671)
    assert np.isclose(correlated,.25-np.arcsin(.1671)/(2*np.pi),atol=1e-10)
