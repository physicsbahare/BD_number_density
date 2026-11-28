"""Empirical aperture calibration and PSF evaluation for local Ali tile caches."""
from __future__ import annotations

import argparse
import time
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy.units as u
from photutils.aperture import CircularAperture
from scipy.ndimage import map_coordinates
from scipy.special import ndtri
from scipy.optimize import least_squares
from src.depth_reconstruction import aperture_kernel,read_header

BANDS=['F115W','F150W','F277W','F814W','F444W']
CAT_BANDS=['F814W','F115W','F150W','F277W','F444W']


def padded_kernel(radius,shape=15):
    k=aperture_kernel(radius,.03)
    if shape<k.shape[0]:raise ValueError('Aperture exceeds cached scene')
    return np.pad(k,(shape-k.shape[0])//2)


def psf_basis(path):
    with fits.open(path,memmap=True) as h:
        header=h[1].header.copy()
        basis=np.array(h[1].data['PSF_MASK'][0],float)
    if header['POLDEG1']!=1 or basis.shape[0]!=3:
        raise ValueError('Only the available first-degree PSFEx basis is supported')
    samp=header['PSF_SAMP']
    half=int(np.floor((basis.shape[1]-1)*samp/2))
    yy,xx=np.mgrid[-half:half+1,-half:half+1]
    native=np.array([map_coordinates(b,[(yy/samp+(basis.shape[1]-1)/2).ravel(),
                                      (xx/samp+(basis.shape[2]-1)/2).ravel()],
                                     order=3,mode='constant').reshape(xx.shape) for b in basis])
    return header,native


def psf_at(header,basis,xy,clip=False):
    x,y=xy
    # PSFEx uses one-based X_IMAGE/Y_IMAGE context coordinates.
    cx=(x+1-header['POLZERO1'])/header['POLSCAL1']
    cy=(y+1-header['POLZERO2'])/header['POLSCAL2']
    p=basis[0]+cx*basis[1]+cy*basis[2]
    if clip:p=np.maximum(p,0)
    return p/p.sum()


def encircled_map(header,basis,xy,radius,kernel_power=1):
    half=(basis.shape[1]-1)//2
    k=aperture_kernel(radius,.03)
    sh=k.shape[0]//2
    coeff=(basis[:,half-sh:half+sh+1,half-sh:half+sh+1]*k**kernel_power).sum(axis=(1,2))
    total=basis.sum(axis=(1,2))
    context=np.column_stack([np.ones(len(xy)),
                             (xy[:,0]+1-header['POLZERO1'])/header['POLSCAL1'],
                             (xy[:,1]+1-header['POLZERO2'])/header['POLSCAL2']])
    ee=(context@coeff)/(context@total)
    if np.any((ee<=0)|(ee>1.05)):raise ValueError('Invalid PSF encircled energy')
    return ee


def robust_sigma(flux):
    q=np.percentile(flux,[16,50,84])
    return float(q[1]),float((q[1]-q[0])/(-ndtri(.16)))


def noise_curve(flux,weight,valid,nbins=8):
    """Measure noise as a free function of WHT, without imposing inverse variance."""
    ids=np.flatnonzero(valid)
    ids=ids[np.argsort(weight[ids])]
    rows=[]
    for sub in np.array_split(ids,nbins):
        if len(sub)<50:continue
        mu,sig=robust_sigma(flux[sub])
        rows.append([float(np.median(weight[sub])),mu,sig,len(sub)])
    rows=np.array(rows,float)
    if len(rows)<2:raise ValueError('Insufficient valid blank apertures')
    if np.any(rows[:,2]<=0):raise ValueError('Nonpositive empirical noise')
    return rows


def eval_curve(rows,weight,column):
    # Collapse repeated median weights before interpolation.
    x,inv=np.unique(np.log(rows[:,0]),return_inverse=True)
    vals=np.array([np.mean(rows[inv==i,column]) for i in range(len(x))])
    w=np.asarray(weight,float)
    if column==2:
        y=np.log(vals)
        if len(x)>1:
            slope,intercept=np.polyfit(x,y,1)
            out=np.exp(np.interp(np.log(np.maximum(w,1e-300)),x,y))
            outside=(np.log(np.maximum(w,1e-300))<x[0])|(np.log(np.maximum(w,1e-300))>x[-1])
            out[outside]=np.exp(np.clip(intercept+slope*np.log(np.maximum(w[outside],1e-300)),-30,40))
        else:out=np.full(w.shape,vals[0])
    else:out=np.interp(np.log(np.maximum(w,1e-300)),x,vals)
    return np.where(np.isfinite(w)&(w>0),out,np.nan)


def calibrate_tile(tile,root,psfdir,radius):
    pos=np.load(root/tile/'positions.npz')
    n=len(pos['blank_xy'])
    clear=np.load(root/tile/'blank_seg_clear.npy')
    with fits.open('/home/bahareh/Data/Cosmos_Web/COSMOSWeb_mastercatalog_v1.1.fits',memmap=True) as h:
        ids=np.array(h[1].data['id']);lookup={int(v):i for i,v in enumerate(ids)}
        ci=np.array([lookup[int(v)] for v in pos['catalog_id']])
        catalog_sky=SkyCoord(np.asarray(h[1].data['ra'][ci])*u.deg,np.asarray(h[1].data['dec'][ci])*u.deg,frame='icrs')
    refhead=fits.Header.fromstring(json.loads((root/tile/'COMPLETE.json').read_text())['science_wcs'],sep='\n')
    random_phase=np.random.default_rng(20261200+int(tile[1:])).uniform(-.5,.5,(n,2))
    arrays={};report={}
    for band in BANDS:
        with np.load(root/tile/(band+'.npz')) as data:
            sci=np.array(data['sci_stamps'],float)
            err=np.array(data['err_stamps'],float)
            wht=np.array(data['wht_stamps'],float)
            sg=np.array(data['sci_grid']);eg=np.array(data['err_grid']);wg=np.array(data['wht_grid'])
            native_grid=np.array(data['native_grid_xy']) if 'native_grid_xy' in data else pos['grid_xy']
            native_stamps=np.array(data['native_stamp_xy']) if 'native_stamp_xy' in data else np.vstack([pos['blank_xy'],pos['catalog_xy']])
        if band=='F814W':
            nh=read_header(next(Path('/mnt/E_202/JWST/COSMOS_Web/HST_F814W_v0.5.2026jul/tiles').glob('*30mas_tile_'+tile+'_*_drz.fits.gz')))
        else:nh=refhead
        fx,fy=WCS(nh).celestial.world_to_pixel(catalog_sky)
        catalog_phase=np.column_stack([fx,fy])-native_stamps[n:]
        if np.max(np.abs(catalog_phase))>1.5:raise ValueError('Source positions exceed the cached subpixel region')
        registration=np.zeros(2)
        if band=='F814W' and radius<=.101:
            sn_fit=pos['catalog_snr'][:,CAT_BANDS.index(band)]
            pair_fit=(np.isfinite(sn_fit)&(sn_fit>-5)&(sn_fit<30)&
                np.isfinite(sci[n:]).all(axis=(1,2))&(wht[n:]>0).all(axis=(1,2))&
                (np.arange(len(sn_fit))%2==0))
            images=sci[n:][pair_fit]*np.sqrt(wht[n:][pair_fit]);cp_fit=catalog_phase[pair_fit];sf=sn_fit[pair_fit]
            ux,uy=WCS(nh).celestial.world_to_pixel_values(catalog_sky.ra.deg,catalog_sky.dec.deg)
            seed=np.clip(np.median(np.column_stack([ux-fx,uy-fy]),axis=0),-1.4,1.4)
            k0=np.array([CircularAperture((7+p[0]+seed[0],7+p[1]+seed[1]),radius/.03).to_mask(method='exact').to_image((15,15)) for p in cp_fit])
            q0=np.sum(images*k0,axis=(1,2));hi=sf>5;a0=np.median(q0[hi]/sf[hi])
            def astrometric_residual(par):
                kernels=np.array([CircularAperture((7+p[0]+par[0],7+p[1]+par[1]),radius/.03).to_mask(method='exact').to_image((15,15)) for p in cp_fit])
                return np.sum(images*kernels,axis=(1,2))/np.exp(par[2])-par[3]-sf
            registration_fit=least_squares(astrometric_residual,[seed[0],seed[1],np.log(a0),0.],
                bounds=([-1.5,-1.5,np.log(a0/4),-5.],[1.5,1.5,np.log(a0*4),5.]),
                loss='soft_l1',f_scale=.2,max_nfev=60,diff_step=.002)
            registration=registration_fit.x[:2]
            catalog_phase=catalog_phase+registration
        phases=np.vstack([random_phase,catalog_phase])
        phases[:n]+=registration
        k=np.array([CircularAperture((7+p[0],7+p[1]),radius/.03).to_mask(method='exact').to_image((15,15)) for p in phases])
        footprint=np.any(k>0,axis=0)
        wt=wht[:,7,7]
        qflux=np.sum(sci*np.sqrt(np.maximum(wht,0))*k,axis=(1,2))
        flux=np.divide(qflux,np.sqrt(wt),out=np.full(len(wt),np.nan),where=wt>0)
        finite=np.isfinite(sci[:,footprint]).all(axis=1)&np.isfinite(err[:,footprint]).all(axis=1)
        coverage=np.all(wht[:,footprint]>0,axis=1)&np.all(err[:,footprint]>0,axis=1)
        valid=clear&finite[:n]&coverage[:n]&np.isfinite(wt[:n])&(wt[:n]>0)
        train=valid.copy();train[np.arange(n)%2==1]=False
        test=valid&~train
        curve=noise_curve(flux[:n],wt[:n],train)
        sig=eval_curve(curve,wt,2);mu=eval_curve(curve,wt,1)
        z=(flux[:n][test]-mu[:n][test])/sig[:n][test]
        sn=pos['catalog_snr'][:,CAT_BANDS.index(band)]
        fc=flux[n:]
        cp=finite[n:]&coverage[n:]&np.isfinite(sn)&(sn>5)&(sn<30)&(fc>0)
        # This is a measured cross-check, not an assumption that AUTO defines SNR.
        ratio=fc[cp]/sn[cp]/sig[n:][cp]
        # Native catalog SNR is reproduced by the noise-equalized aperture
        # sum, not by an AUTO error or a rescaled empirically calibrated error.
        pair=finite[n:]&coverage[n:]&np.isfinite(sn)&(sn>-5)&(sn<30)
        train_pair=pair&(np.arange(len(sn))%2==0);test_pair=pair&~train_pair
        initial=np.median(qflux[n:][cp]/sn[cp])
        fit=least_squares(lambda par:qflux[n:][train_pair]/np.exp(par[0])-par[1]-sn[train_pair],
            [np.log(initial),0.],loss='soft_l1',f_scale=.2)
        scale=float(np.exp(fit.x[0]));offset=float(fit.x[1]*scale)
        validation=(qflux[n:][test_pair]-offset)/scale-sn[test_pair]
        ref_scale=float(np.median((scale/np.sqrt(wt[n:][cp]))/sig[n:][cp]))
        ref_scatter=float((np.percentile(ratio,84)-np.percentile(ratio,16))/2)
        psfname='hst_f814w' if band=='F814W' else band.lower()
        psfpath=next(psfdir.glob('*_'+psfname+'_'+tile+'_*.psf'))
        ph,pb=psf_basis(psfpath)
        ee=encircled_map(ph,pb,native_grid,radius)
        gridvalid=np.isfinite(sg)&np.isfinite(eg)&(eg>0)&np.isfinite(wg)&(wg>0)
        sigma_grid=eval_curve(curve,wg,2)
        mean_grid=eval_curve(curve,wg,1)-np.divide(offset,np.sqrt(wg),out=np.full(len(wg),np.nan),where=wg>0)
        arrays[band+'_valid']=gridvalid
        arrays[band+'_sigma_noise_njy']=sigma_grid
        arrays[band+'_sigma_snr_njy']=np.divide(scale,np.sqrt(wg),out=np.full(len(wg),np.nan),where=wg>0)
        arrays[band+'_background_njy']=mean_grid
        arrays[band+'_ee']=ee
        arrays[band+'_test_z']=z
        arrays[band+'_test_flux_njy']=flux[:n][test]
        arrays[band+'_test_sigma_njy']=sig[:n][test]
        arrays[band+'_test_background_njy']=mu[:n][test]
        arrays[band+'_test_xy']=native_stamps[:n][test]
        if band=='F814W':
            arrays['F814W_test_stamps_njy']=sci[:n][test]
            arrays['F814W_test_var_njy2']=err[:n][test]**2
            arrays['F814W_test_wht']=wt[:n][test]
            arrays['F814W_test_weight_ratio']=np.sqrt(np.maximum(wht[:n][test],0))/np.sqrt(wt[:n][test])[:,None,None]
            arrays['F814W_test_aperture_kernel']=k[:n][test]
            arrays['F814W_test_phase_xy']=phases[:n][test]
            arrays['F814W_test_snr_denom_njy']=scale/np.sqrt(wt[:n][test])
            arrays['F814W_test_snr_offset_njy']=offset/np.sqrt(wt[:n][test])
            extra=err**2-np.divide(1.,wht,out=np.full(wht.shape,np.nan),where=wht>0)
            shot_good=np.isfinite(sci)&np.isfinite(extra)&(sci>5*err)&(extra>0)&(wht>0)
            shot_coeff=(extra*wht/np.where(sci>0,sci,np.nan))[shot_good]
            arrays['F814W_shot_coefficient']=np.median(shot_coeff)
        report[band]=dict(n_train=int(train.sum()),n_test=int(test.sum()),
            aperture_radius_arcsec=radius,noise_curve_wht_mean_sigma_n=curve.tolist(),
            weight_log_slope=float(np.polyfit(np.log(curve[:,0]),np.log(curve[:,2]),1)[0]),
            blank_z_p016_p50_p84_p977=np.percentile(z,[16,50,84,97.7249868]).tolist(),
            empirical_blank_p_snr_lt_2=float(np.mean(z<2)),
            catalog_n=int(cp.sum()),catalog_sigma_over_empirical_sigma_median=ref_scale,
            catalog_sigma_ratio_scatter=ref_scatter,
            catalog_snr_estimator='(sum(SCI_nJy * sqrt(WHT_native) * fractional_aperture) - offset) / scale',
            catalog_noise_equalized_scale=scale,catalog_noise_equalized_offset=offset,
            catalog_snr_fit_n=int(train_pair.sum()),catalog_snr_validation_n=int(test_pair.sum()),
            catalog_snr_validation_residual_p16_p50_p84=np.percentile(validation,[16,50,84]).tolist(),
            catalog_snr_validation_median_absolute_error=float(np.median(np.abs(validation))),
            native_static_registration_offset_xy_pixels=registration.tolist(),
            catalog_sigma_ratio_p16_p50_p84=np.percentile(ratio,[16,50,84]).tolist(),
            ee_p0_p50_p100=np.percentile(ee,[0,50,100]).tolist(),psf_path=str(psfpath))
        if band=='F814W':
            report[band]['source_variance_fit']=dict(n_bright_pixels=int(shot_good.sum()),
                coefficient_p16_p50_p84=np.percentile(shot_coeff,[16,50,84]).tolist(),
                definition='ERR^2 - 1/WHT = coefficient * SCI/WHT; calibrated on bright source pixels; electrons_per_njy = WHT/coefficient')
    arrays['global_xy']=pos['global_xy'];arrays['grid_xy']=pos['grid_xy']
    return arrays,report


def run(args):
    args.out.mkdir(parents=True,exist_ok=True)
    summary={}
    pending=set(args.tiles)
    while pending:
        for tile in sorted(pending):
            if not (args.cache/tile/'COMPLETE.json').exists():continue
            if (args.out/(tile+'.json')).exists() and (args.out/(tile+'.npz')).exists():
                summary[tile]=json.loads((args.out/(tile+'.json')).read_text());pending.remove(tile);continue
            print('Calibrate:',tile,flush=True)
            arrays,report=calibrate_tile(tile,args.cache,args.psfdir,args.radius)
            np.savez_compressed(args.out/(tile+'.npz'),**arrays)
            summary[tile]=report;pending.remove(tile)
            (args.out/(tile+'.json')).write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
            (args.out/'calibration.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
        if not args.watch:break
        if pending:time.sleep(15)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cache',type=Path,default=Path('/home/bahareh/results/ali_fiveband_reconstruction'))
    p.add_argument('--psfdir',type=Path,default=Path('/mnt/E_202/JWST/COSMOS-Web/auxiliary/psfs_v5.0'))
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--radius',type=float,default=.1)
    p.add_argument('--tiles',nargs='+',default=[f'{a}{i}' for a in 'AB' for i in range(1,11)])
    p.add_argument('--watch',action='store_true')
    run(p.parse_args())
