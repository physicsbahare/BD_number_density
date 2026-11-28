"""Validate F814W recovery at fluxes taken directly from the adopted population."""
import argparse,json,time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.ndimage import shift
from scipy.spatial import cKDTree
from scipy.special import ndtr
from src.ali_depth_analysis import psf_basis,psf_at
from src.f814w_targeted_recovery import wilson
from injection_recovery.core import inject_psf_njy


def validate(tile,args):
    report=json.loads((args.calibration/(tile+'.json')).read_text());r=report['F814W']
    with np.load(args.calibration/(tile+'.npz')) as d:
        scene=d['F814W_test_stamps_njy'];var=d['F814W_test_var_njy2'];w=d['F814W_test_weight_ratio']
        kernel=d['F814W_test_aperture_kernel'];phase=d['F814W_test_phase_xy'];xy=d['F814W_test_xy']
        sig=d['F814W_test_sigma_njy'];mu=d['F814W_test_background_njy'];den=d['F814W_test_snr_denom_njy'];off=d['F814W_test_snr_offset_njy']
        wht=d['F814W_test_wht'];shot=float(d['F814W_shot_coefficient'])
        good=d['F115W_valid']&d['F150W_valid']&d['F277W_valid']&d['F444W_valid']
        distance,ind=cKDTree(d['grid_xy'][good]).query(xy+np.array(r['native_static_registration_offset_xy_pixels']))
        indices=np.flatnonzero(good)[ind];nir_sigma=d['F444W_sigma_noise_njy'][indices];nir_ee=d['F444W_ee'][indices]
    table=pd.read_csv(Path(__file__).resolve().parents[1]/'data/elfowl_fiveband_bin_reach.csv')
    ph,pb=psf_basis(Path(r['psf_path']));rng=np.random.default_rng(20261300+int(tile[1:])+(100 if tile[0]=='B' else 0))
    rows=[];f25=10**((31.4-25)/2.5)
    for i in np.flatnonzero((distance<=150)&np.isfinite(nir_sigma)&(nir_sigma>0)):
        model=table.iloc[i%len(table)];s444=[5.,8.,15.][i%3]
        dpc=float(model.d25_F444W_pc*np.sqrt(f25*nir_ee[i]/(s444*nir_sigma[i])))
        flux=float(f25*(model.d25_F814W_pc/dpc)**2)
        p=psf_at(ph,pb,xy[i]);p=shift(p,phase[i][::-1],order=1,mode='constant',cval=0.);p/=p.sum()
        c=p.shape[0]//2;crop=p[c-7:c+8,c-7:c+8]
        if -crop[crop<0].sum()>1e-4:raise ValueError('Invalid PSF core')
        crop=np.maximum(crop,0);mass=crop.sum();k=kernel[i]*w[i];gain=wht[i]/shot
        injected,_=inject_psf_njy(scene[i],var[i],crop,7.,7.,flux*mass,gain,rng)
        measured=float(np.sum(injected*k));snr=(measured-off[i])/den[i]
        true_ap=float(flux*np.sum(crop*k));sv=float(flux*np.sum(crop*k*k)/gain)
        pred=float(ndtr((2*den[i]+off[i]-mu[i]-true_ap)/np.sqrt(sig[i]**2+sv)))
        rows.append([int(model.teff_lo_K),dpc,s444,flux,true_ap/sig[i],measured,snr,int(snr<2),pred])
    a=np.array(rows);summary=[]
    for lo in table.teff_lo_K:
        v=a[a[:,0]==lo];n=len(v)
        if not n:continue
        npass=int(v[:,7].sum());p=npass/n
        summary.append(dict(tile=tile,teff_lo_K=int(lo),teff_hi_K=int(lo+150),n=n,n_pass=npass,
            empirical_probability=p,gaussian_probability=float(v[:,8].mean()),wilson95=wilson(npass,n),
            injected_total_F814W_njy_p16_p50_p84=np.percentile(v[:,3],[16,50,84]).tolist(),
            true_F814W_aperture_snr_p16_p50_p84=np.percentile(v[:,4],[16,50,84]).tolist()))
    np.savez_compressed(args.out/(tile+'_population_trials.npz'),
        columns=np.array(['teff_lo_K','distance_pc','true_F444W_snr','F814W_total_njy','F814W_aperture_true_snr','recovered_aperture_njy','measured_SNR_F814W','pass_signed_lt_2','gaussian_probability']),trials=a)
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--calibration',type=Path,default=Path('/home/bahareh/results/ali_calibration_final'))
    p.add_argument('--out',type=Path,default=Path('/home/bahareh/results/ali_f814w_population_validation'))
    p.add_argument('--watch',action='store_true')
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True);pending={f'{t}{i}' for t in 'AB' for i in range(1,11)};allrows=[]
    while pending:
        for tile in sorted(pending):
            marker=a.out/(tile+'.json')
            if not (a.calibration/(tile+'.json')).exists():continue
            if marker.exists():rows=json.loads(marker.read_text())
            else:
                print('Population-flux F814W validation',tile,flush=True);rows=validate(tile,a);marker.write_text(json.dumps(rows,indent=2)+'\n')
            allrows.extend(rows);pending.remove(tile)
            (a.out/'population_validation_summary.json').write_text(json.dumps(dict(records=allrows,
                definition='Existing thirteen-bin atmosphere medians; F444W true SNR 5, 8, 15 at local NIRCam depth; real held-out F814W backgrounds; source shot noise included.'),indent=2)+'\n')
        if not a.watch:break
        if pending:time.sleep(15)
