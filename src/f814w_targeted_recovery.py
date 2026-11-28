"""Targeted native-aperture F814W injection/recovery on real held-out backgrounds.

This does not measure catalog detection, deblending, morphology, visual cuts,
or quality flags. Cropped PSF flux is retained explicitly: injections into the
15x15 observed scenes have the same in-aperture flux as a full PSF injection.
"""
from __future__ import annotations
import argparse,json,time
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.ndimage import shift
from scipy.special import ndtr
from photutils.aperture import CircularAperture
from src.ali_depth_analysis import psf_basis,psf_at,encircled_map,padded_kernel
from injection_recovery.core import inject_psf_njy

LEVELS=np.array([0.,.001,.01,.1,.5,1.,2.,3.,5.])


def wilson(k,n):
    z=1.95996398454;p=k/n;den=1+z*z/n
    mid=(p+z*z/(2*n))/den;half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return [float(mid-half),float(mid+half)]


def star_check(tile,cache,ph,pb):
    pos=np.load(cache/tile/'positions.npz');n=len(pos['blank_xy'])
    with fits.open('/home/bahareh/Data/Cosmos_Web/COSMOSWeb_mastercatalog_v1.1.fits',memmap=True) as h:
        d=h[1].data;lookup={int(v):i for i,v in enumerate(d['id'])}
        ci=np.array([lookup[int(v)] for v in pos['catalog_id']])
        select=(h[2].data['type'][ci]==1)&(~d['flag_star'][ci])&(d['snr_hst-f814w'][ci]>30)&(d['snr_hst-f814w'][ci]<2000)&(d['fwhm'][ci]<.6)
    with np.load(cache/tile/'F814W.npz') as d:
        stars=d['sci_stamps'][n:][select]
        xy=d['native_stamp_xy'][n:][select] if 'native_stamp_xy' in d else pos['catalog_xy'][select]
    control_ids=pos['catalog_id'][select]
    out=[];yy,xx=np.mgrid[:15,:15];r=np.hypot(xx-7,yy-7)
    for i,a in enumerate(stars):
        if not np.isfinite(a).all():continue
        back=np.median(a[r>6.5]);v=np.maximum(a-back,0)*(r<4)
        if v.sum()<=0:continue
        cx,cy=(xx*v).sum()/v.sum(),(yy*v).sum()/v.sum()
        if np.hypot(cx-7,cy-7)>2:continue
        ksmall=CircularAperture((cx,cy),.1/.03).to_mask(method='exact').to_image((15,15))
        klarge=CircularAperture((cx,cy),.18/.03).to_mask(method='exact').to_image((15,15))
        small=np.sum((a-back)*ksmall);large=np.sum((a-back)*klarge)
        if large<=0:continue
        pp=psf_at(ph,pb,xy[i]);pc=(pp.shape[0]-1)/2
        from scipy.ndimage import map_coordinates
        model_image=map_coordinates(pp,[yy-cy+pc,xx-cx+pc],order=3)
        mb=np.median(model_image[r>6.5]);model=np.sum((model_image-mb)*ksmall)/np.sum((model_image-mb)*klarge)
        out.append(dict(tile=tile,catalog_id=int(control_ids[i]),centroid_dx=float(cx-7),centroid_dy=float(cy-7),
             observed_ee010_over_ee018=float(small/large),model_ee010_over_ee018=float(model),ratio=float(small/large/model),
             control_definition='LePHARE type=1, unmasked, 30<SNR_F814W<2000, FWHM<0.6 arcsec; native core EEF comparison with identical annular background treatment'))
    return out


def recover_tile(tile,args):
    f=args.calibration/(tile+'.npz');r=json.loads((args.calibration/(tile+'.json')).read_text())['F814W']
    ph,pb=psf_basis(Path(r['psf_path']))
    star_rows=star_check(tile,args.cache,ph,pb)
    with np.load(f) as d:
        scenes=np.array(d['F814W_test_stamps_njy']);variance=np.array(d['F814W_test_var_njy2'])
        sigma=d['F814W_test_sigma_njy'];mean=d['F814W_test_background_njy'];xy=d['F814W_test_xy'];wht=d['F814W_test_wht']
        shot=float(d['F814W_shot_coefficient'])
        weight_ratio=d['F814W_test_weight_ratio'];kernels=d['F814W_test_aperture_kernel'];phases=d['F814W_test_phase_xy']
        snr_denom=d['F814W_test_snr_denom_njy'];snr_offset=d['F814W_test_snr_offset_njy']
    factor=r['catalog_sigma_over_empirical_sigma_median']
    if args.max_backgrounds:
        n=min(len(scenes),args.max_backgrounds);scenes,variance,sigma,mean,xy,wht,weight_ratio,kernels,phases,snr_denom,snr_offset=[a[:n] for a in [scenes,variance,sigma,mean,xy,wht,weight_ratio,kernels,phases,snr_denom,snr_offset]]
    n=len(scenes);rng=np.random.default_rng(20261126+int(tile[1:])+(100 if tile[0]=='B' else 0))
    recovered=np.zeros((len(LEVELS),n));injected=np.zeros_like(recovered);expect=np.zeros_like(recovered)
    analytic=np.zeros_like(recovered);source_only=np.zeros_like(recovered);ee_values=np.zeros(n)
    baseline=np.sum(scenes*kernels*weight_ratio,axis=(1,2))
    for i in range(n):
        k=kernels[i]*weight_ratio[i]
        p=psf_at(ph,pb,xy[i]);phase=phases[i][::-1]
        p=shift(p,shift=phase,order=1,mode='constant',cval=0.)
        p/=p.sum();half=p.shape[0]//2;crop=p[half-7:half+8,half-7:half+8]
        negative_mass=float(-crop[crop<0].sum())
        if negative_mass>1e-4:raise ValueError('PSF has appreciable negative flux inside injection scene')
        crop=np.maximum(crop,0);mass=crop.sum();ee=float(np.sum(crop*k));ee_values[i]=ee
        gain=float(wht[i]/shot)
        if gain<=0 or not np.isfinite(gain):raise ValueError('Invalid empirically inferred source-noise gain')
        for j,level in enumerate(LEVELS):
            total_flux=float(level*sigma[i]/ee)
            image,var=inject_psf_njy(scenes[i],variance[i],crop,7.,7.,total_flux*mass,gain,rng)
            measured=float(np.sum(image*k));recovered[j,i]=measured;injected[j,i]=total_flux
            expect[j,i]=total_flux*ee;source_only[j,i]=measured-baseline[i]-expect[j,i]
            added_variance=float(total_flux*np.sum(crop*k*k)/gain)
            analytic[j,i]=ndtr((2*snr_denom[i]+snr_offset[i]-mean[i]-expect[j,i])/np.sqrt(sigma[i]**2+added_variance))
    snr=(recovered-snr_offset[None,:])/snr_denom[None,:]
    passed=snr<2.
    records=[]
    for j,level in enumerate(LEVELS):
        successes=int(passed[j].sum());p=successes/n
        records.append(dict(tile=tile,source_aperture_snr_true=float(level),n=n,n_pass=successes,
            empirical_probability=p,binomial_standard_error=float(np.sqrt(p*(1-p)/n)),wilson95=wilson(successes,n),
            gaussian_probability=float(analytic[j].mean()),difference_empirical_minus_gaussian=float(p-analytic[j].mean()),
            injected_total_flux_njy_p16_p50_p84=np.percentile(injected[j],[16,50,84]).tolist(),
            recovered_aperture_flux_njy_mean=float(recovered[j].mean()),
            aperture_flux_bias_njy=float(np.mean(recovered[j]-mean-expect[j])),
            source_only_flux_bias_njy=float(source_only[j].mean()),
            source_only_bias_standard_error=float(source_only[j].std(ddof=1)/np.sqrt(n)),
            convergence_probability_by_n={str(m):float(passed[j,:m].mean()) for m in [250,1000,n] if m<=n}))
    np.savez_compressed(args.out/(tile+'_recovery.npz'),levels=LEVELS,injected_total_njy=injected,
        expected_aperture_njy=expect,recovered_aperture_njy=recovered,measured_snr=snr,
        gaussian_probability=analytic,ee=ee_values,blank_sigma_njy=sigma,blank_mean_njy=mean)
    return records,star_rows


def run(args):
    args.out.mkdir(parents=True,exist_ok=True);pending=set(args.tiles);allrows=[];stars=[]
    while pending:
        for tile in sorted(pending):
            if not (args.calibration/(tile+'.json')).exists():continue
            marker=args.out/(tile+'.json')
            if marker.exists():
                saved=json.loads(marker.read_text());records,st=saved['recovery'],saved['stars']
            else:
                print('F814W recovery',tile,flush=True);records,st=recover_tile(tile,args)
                marker.write_text(json.dumps(dict(recovery=records,stars=st),indent=2)+'\n')
            allrows.extend(records);stars.extend(st);pending.remove(tile)
            (args.out/'recovery_summary.json').write_text(json.dumps(dict(status='TARGETED_APERTURE_RECOVERY_ONLY',
                tiles=sorted(set(x['tile'] for x in allrows)),records=allrows,star_controls=stars,
                warning='No catalog detection, blending, morphology, visual-selection or quality-flag completeness is measured.'),indent=2)+'\n')
        if not args.watch:break
        if pending:time.sleep(15)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--calibration',type=Path,default=Path('/home/bahareh/results/ali_calibration_r010'))
    p.add_argument('--cache',type=Path,default=Path('/home/bahareh/results/ali_fiveband_reconstruction'))
    p.add_argument('--out',type=Path,default=Path('/home/bahareh/results/ali_f814w_recovery'))
    p.add_argument('--tiles',nargs='+',default=[f'{a}{i}' for a in 'AB' for i in range(1,11)])
    p.add_argument('--max-backgrounds',type=int,default=0)
    p.add_argument('--watch',action='store_true')
    run(p.parse_args())
