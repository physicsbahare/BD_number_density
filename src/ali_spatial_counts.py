"""Signed selection probabilities and existing Galactic component integrals."""
from __future__ import annotations

import numpy as np
from scipy.special import ndtr
from scipy.integrate import quad_vec
from numpy.polynomial.legendre import leggauss
from src.cosmosweb_counts import (galactic_lb, disk_rho, halo_rho,
                                local_component_densities)


def signed_probability(flux, sigma_noise, sigma_snr, threshold, detection):
    f,n,s=np.broadcast_arrays(flux,sigma_noise,sigma_snr)
    valid=np.isfinite(f)&np.isfinite(n)&np.isfinite(s)&(n>0)&(s>0)
    z=np.zeros(f.shape,float)
    np.divide(f-threshold*s,n,out=z,where=valid)
    p=ndtr(z if detection else -z)
    return np.where(valid,p,0.)


def density_components(distance,rho):
    l,b=galactic_lb()
    thin,thick,halo=local_component_densities(rho)
    return np.array([disk_rho(distance,l,b,thin,2600.,187.),
                     disk_rho(distance,l,b,thick,3600.,900.),
                     halo_rho(distance,l,b,halo)])


def joint_veto_detection(fv,fd,nv,nd,sv,sd,correlation):
    """Bivariate Gaussian veto/detection sensitivity via Plackett's identity."""
    a=(2*sv-fv)/nv;b=(5*sd-fd)/nd
    base=ndtr(a)*ndtr(-b)
    if correlation==0:return base
    nodes,weights=leggauss(12);rs=(nodes+1)*correlation/2
    correction=sum(w*np.exp(-(a*a-2*r*a*b+b*b)/(2*(1-r*r)))/(2*np.pi*np.sqrt(1-r*r)) for r,w in zip(rs,weights))*correlation/2
    return np.clip(base-correction,0,1)


def integrate_cells(d25,rho,sigma_noise,sigma_snr,encircled,area_sr,
                    rmax=50000.,epsrel=1e-6,background_mean=None,empirical_cdf=None,radial_domain='depth',f814w_source_coefficient=None,nircam_correlation=0.,depth_sigma=5.,boundary_flux_scale=1.):
    """Adaptive line-of-sight quadrature, preserving component normalization.

    Arrays are cell x band, ordered blue bands then F444W. Distances and
    area are never adjusted to reproduce a historical scalar answer.
    """
    sigma_noise=np.asarray(sigma_noise,float)
    sigma_snr=np.asarray(sigma_snr,float)
    ee=np.broadcast_to(encircled,sigma_noise.shape)
    mu=0 if background_mean is None else np.asarray(background_mean)
    # 25 AB in nJy, consistent with the original zp25 surface-flux tables.
    f25=10**((31.4-25.)/2.5)
    amp=f25*np.asarray(d25)**2*ee
    if radial_domain in ('depth','catalog_depth'):
        limit_noise=sigma_noise[:,-1] if radial_domain=='depth' else sigma_snr[:,-1]
        if not np.isfinite(depth_sigma) or depth_sigma<=0:raise ValueError('Depth threshold must be positive')
        boundary_flux_scale=np.asarray(boundary_flux_scale)
        if np.any(~np.isfinite(boundary_flux_scale)) or np.any(boundary_flux_scale<=0):raise ValueError('Boundary flux scale must be positive')
        reach=np.sqrt(amp[:,-1]/(depth_sigma*limit_noise*boundary_flux_scale))
        def depth_integrand(t):
            if t<=0:return np.zeros(3)
            r=t*reach
            flux=amp/r[:,None]**2+mu
            p=np.ones(len(sigma_noise))
            if nircam_correlation:p*=joint_veto_detection(flux[:,2],flux[:,-1],sigma_noise[:,2],sigma_noise[:,-1],sigma_snr[:,2],sigma_snr[:,-1],nircam_correlation)
            for k in range(sigma_noise.shape[1]):
                if nircam_correlation and k in [2,sigma_noise.shape[1]-1]:continue
                det=k==sigma_noise.shape[1]-1
                shot=None
                if f814w_source_coefficient is not None and k==sigma_noise.shape[1]-2:
                    shot=np.sqrt(np.maximum(amp[:,k]/r**2*f814w_source_coefficient,0))
                if empirical_cdf is not None and k==sigma_noise.shape[1]-2:
                    z=(2*sigma_snr[:,k]-flux[:,k])/sigma_noise[:,k]
                    p*=empirical_cdf(z) if shot is None else empirical_cdf(z,shot/sigma_noise[:,k])
                else:p*=signed_probability(flux[:,k],sigma_noise[:,k] if shot is None else np.hypot(sigma_noise[:,k],shot),sigma_snr[:,k],5. if det else 2.,det)
            return np.sum(density_components(r,rho)*(r*r*reach*p*area_sr)[None,:],axis=1)
        value,error=quad_vec(depth_integrand,0.,1.,epsabs=1e-10,epsrel=epsrel,limit=1000,norm='max')
        return value,float(error)
    if radial_domain!='extended':raise ValueError('Unknown radial domain')
    def integrand(logr):
        r=np.exp(logr)
        flux=amp/(r*r)+mu
        p=np.ones(len(sigma_noise))
        if nircam_correlation:p*=joint_veto_detection(flux[:,2],flux[:,-1],sigma_noise[:,2],sigma_noise[:,-1],sigma_snr[:,2],sigma_snr[:,-1],nircam_correlation)
        for k in range(sigma_noise.shape[1]):
            if nircam_correlation and k in [2,sigma_noise.shape[1]-1]:continue
            det=(k==sigma_noise.shape[1]-1)
            shot=None
            if f814w_source_coefficient is not None and k==sigma_noise.shape[1]-2:
                shot=np.sqrt(np.maximum(amp[:,k]/r**2*f814w_source_coefficient,0))
            if empirical_cdf is not None and k==sigma_noise.shape[1]-2:
                z=(2*sigma_snr[:,k]-flux[:,k])/sigma_noise[:,k]
                p*=empirical_cdf(z) if shot is None else empirical_cdf(z,shot/sigma_noise[:,k])
            else:p*=signed_probability(flux[:,k],sigma_noise[:,k] if shot is None else np.hypot(sigma_noise[:,k],shot),sigma_snr[:,k],5. if det else 2.,det)
        return density_components(r,rho)*(r**3*float(np.sum(area_sr*p)))
    value,error=quad_vec(integrand,np.log(.1),np.log(rmax),epsabs=1e-8,
                         epsrel=epsrel,limit=1000,norm='max')
    return value,float(error)
