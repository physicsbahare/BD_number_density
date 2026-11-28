"""Validate the radial domain using cached photometry, without new injections."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import quad_vec
from scipy.special import ndtr
from numpy.polynomial.hermite import hermgauss
from numpy.polynomial.legendre import leggauss
from src.ali_count_run import depth_groups
from src.ali_spatial_counts import density_components
from src.cosmosweb_counts import TEFF_DENSITY_2024
from src.f814w_targeted_recovery import wilson


def halo_tail_power(n=2.77):
    """Cumulative count scales as R**(3-n) for a nonzero noise floor."""
    return 3-n


def zero_source_probability(noise,denominator,background):
    noise=np.asarray(noise);denominator=np.asarray(denominator);background=np.asarray(background)
    return np.prod(ndtr((2*denominator[:,:-1]-background[:,:-1])/noise[:,:-1]),axis=1)*ndtr((background[:,-1]-5*denominator[:,-1])/noise[:,-1])


def empirical_blank_check(calibration,cache):
    records=[];total=0;k4=0;k5=0
    for f in sorted(calibration.glob('*.npz')):
        report=json.loads(f.with_suffix('.json').read_text());pos=np.load(cache/f.stem/'positions.npz');n=len(pos['blank_xy']);band_snr={}
        with np.load(f) as d:
            for b in ['F115W','F150W','F277W','F814W','F444W']:
                with np.load(cache/f.stem/(b+'.npz')) as raw:
                    xy=raw['native_stamp_xy'][:n] if 'native_stamp_xy' in raw else pos['blank_xy']
                    lookup={tuple(v):i for i,v in enumerate(xy)};ix=np.array([lookup[tuple(v)] for v in d[b+'_test_xy']]);wt=raw['wht_stamps'][ix,7,7]
                flux=d[b+'_test_flux_njy'];r=report[b]
                snr=(flux*np.sqrt(wt)-r['catalog_noise_equalized_offset'])/r['catalog_noise_equalized_scale']
                band_snr[b]={int(i):float(s) for i,s in zip(ix,snr)}
        common=sorted(set.intersection(*(set(v) for v in band_snr.values())));s=np.array([[band_snr[b][i] for b in ['F115W','F150W','F277W','F814W','F444W']] for i in common])
        p4=(s[:,4]>5)&np.all(s[:,:3]<2,axis=1);p5=p4&(s[:,3]<2);nn=len(common);a=int(p4.sum());b=int(p5.sum())
        total+=nn;k4+=a;k5+=b;records.append(dict(tile=f.stem,n=nn,four_band_blank_pass=a,five_band_blank_pass=b))
    return dict(n=total,four_band_pass=k4,five_band_pass=k5,five_band_probability=k5/total,five_band_wilson95=wilson(k5,total),tiles=records,
        interpretation='No source injected: chance satisfaction of signed native catalog-SNR cuts on real held-out backgrounds. This is not true-source catalog-detection completeness.')


def run(a):
    a.out.mkdir(parents=True,exist_ok=True)
    blankfile=a.out/'native_zero_source_selection.json'
    blank=json.loads(blankfile.read_text()) if blankfile.exists() else empirical_blank_check(a.calibration,a.cache)
    blankfile.write_text(json.dumps(blank,indent=2)+'\n');print('Measured zero-source rate',blank['five_band_probability'],flush=True)
    with np.load(a.maps) as z:
        use=z['five_area_sr']>0;v=z['properties'][use].astype(float);areas=z['five_area_sr'][use];tid=z['tile_id'][use];ftid=z['f814w_tile_id'][use];tiles=z['tiles']
    group=depth_groups(v,areas,tid*len(tiles)+ftid,list(range(10)),a.cells);ng=group.max()+1
    ar=np.bincount(group,weights=areas);means=np.column_stack([np.bincount(group,weights=areas*v[:,k])/ar for k in range(20)]);q=np.bincount(group,weights=areas*v[:,25])/ar
    tg=np.rint(np.bincount(group,weights=areas*ftid)/ar).astype(int);cdfgroups=[]
    for t in np.unique(tg):
        name=str(tiles[t]);report=json.loads((a.calibration/(name+'.json')).read_text())['F814W']
        with np.load(a.calibration/(name+'.npz')) as z:w=z['F814W_test_wht'][::2];res=z['F814W_test_z'][::2]
        edges=np.quantile(w,[.25,.5,.75]);samplestr=np.searchsorted(edges,w);mapstr=np.searchsorted(edges,(report['catalog_noise_equalized_scale']/means[:,8])**2)
        for j in range(4):
            ix=np.flatnonzero((tg==t)&(mapstr==j));sample=np.sort(res[samplestr==j]);cdfgroups.append((ix,sample,np.arange(1,len(sample)+1)/len(sample)))
    hn,hw=hermgauss(15);hw/=np.sqrt(np.pi)
    def optical_cdf(z,shot=None):
        p=np.zeros(len(z))
        for ix,s,prob in cdfgroups:
            if not len(ix):continue
            if shot is None:p[ix]=np.interp(z[ix],s,prob,left=0,right=1)
            else:p[ix]=sum(wh*np.interp(z[ix]-np.sqrt(2)*h*shot[ix],s,prob,left=0,right=1) for h,wh in zip(hn,hw))
        return p
    noise,den,mu,ee=means[:,:5],means[:,5:10],means[:,10:15],means[:,15:20]
    def probability(signal,source_noise=True):
        flux=signal+mu;p=np.prod(ndtr((2*den[:,:3]-flux[:,:3])/noise[:,:3]),axis=1)*ndtr((flux[:,4]-5*den[:,4])/noise[:,4])
        shot=np.sqrt(np.maximum(signal[:,3]*q,0))/noise[:,3] if source_noise else None
        return p*optical_cdf((2*den[:,3]-flux[:,3])/noise[:,3],shot)
    p0=probability(np.zeros_like(noise),False)
    summary=dict(area_deg2=float(ar.sum()*(180/np.pi)**2),cells_per_tile=a.cells,zero_source_probability_area_mean=float(np.average(p0,weights=ar)),zero_source_probability_area_p16_p50_p84=np.percentile(p0,[16,50,84]).tolist(),halo_density_exponent=2.77,cumulative_noise_floor_power=halo_tail_power(),untruncated_integral_converges=False)
    (a.out/'domain_status.json').write_text(json.dumps(summary,indent=2)+'\n')
    table=pd.read_csv(Path(__file__).resolve().parents[1]/'data/elfowl_fiveband_bin_reach.csv');rows=[];f25=10**((31.4-25)/2.5)
    if a.relative_floor_only:
        data=pd.read_csv(a.out/'native_radial_shell_counts.csv');differences=[]
        for i,row in data[data.r_high_pc.isna()].iterrows():
            model=table[table.teff_lo_K==row.teff_lo_K].iloc[0];rho=TEFF_DENSITY_2024[(int(row.teff_lo_K),int(row.teff_hi_K))][0]
            amp444=f25*model.d25_F444W_pc**2*ee[:,4];reach=np.sqrt(amp444/np.maximum(5*den[:,4]-mu[:,4],1e-30))
            results=[];lo=np.log(row.relative_r_low);hi=np.log(row.relative_r_high)
            for nodes in [64,128]:
                x,w=leggauss(nodes);value=np.zeros(3)
                for point,weight in zip((lo+hi)/2+x*(hi-lo)/2,w*(hi-lo)/2):
                    r=np.exp(point)*reach;value+=np.sum(density_components(r,rho)*(r**3*p0*ar*weight)[None,:],axis=1)
                results.append(value)
            for key,val in zip(['noise_floor_thin','noise_floor_thick','noise_floor_halo'],results[1]):data.at[i,key]=val
            data.at[i,'noise_floor_total']=sum(results[1]);data.at[i,'excess_above_no_source_floor']=row.total-sum(results[1]);differences.append(float(np.max(abs(results[1]-results[0]))))
        data.to_csv(a.out/'native_radial_shell_counts.csv',index=False)
        (a.out/'relative_floor_convergence.json').write_text(json.dumps(dict(nodes=[64,128],max_component_difference=max(differences),sum_max_component_differences=sum(differences)),indent=2)+'\n')
        print('Added relative-boundary noise-floor decomposition',max(differences),flush=True);return
    # Absolute shells quantify the cutoff dependence, separated from the
    # no-source forced-photometry floor. No shell cutoff is adopted as physical.
    edges=[.1,2000.,5000.,10000.,20000.,50000.,100000.,200000.,1000000.]
    for _,binrow in table.iterrows():
        lo,hi=int(binrow.teff_lo_K),int(binrow.teff_hi_K);rho,err,status=TEFF_DENSITY_2024[(lo,hi)]
        amp=f25*np.array([binrow['d25_'+b+'_pc']**2 for b in ['F115W','F150W','F277W','F814W','F444W']])*ee
        def integrand(logr):
            r=np.exp(logr);p=probability(amp/r**2);dens=density_components(r,rho)*r**3
            return np.concatenate([dens*float(np.dot(ar,p)),dens*float(np.dot(ar,p0))])
        for low,high in zip(edges[:-1],edges[1:]):
            value,qerr=quad_vec(integrand,np.log(low),np.log(high),epsabs=1e-7,epsrel=a.epsrel,limit=500)
            rows.append(dict(teff_lo_K=lo,teff_hi_K=hi,r_low_pc=low,r_high_pc=high,thin=value[0],thick=value[1],halo=value[2],total=sum(value[:3]),noise_floor_thin=value[3],noise_floor_thick=value[4],noise_floor_halo=value[5],noise_floor_total=sum(value[3:]),excess_above_no_source_floor=sum(value[:3]-value[3:]),quadrature_error=qerr,density_status=status,density_error=None if err is None else sum(value[:3])*err/rho))
        # Same radial integral in scaled intrinsic native-SNR reach coordinates.
        reach=np.sqrt(amp[:,4]/np.maximum(5*den[:,4]-mu[:,4],1e-30))
        def scaled_integrand(logt):
            t=np.exp(logt);r=t*reach;p=probability(amp/r[:,None]**2);d=density_components(r,rho)
            return np.sum(d*(r**3*p*ar)[None,:],axis=1)
        for low,high in [(.001,1.),(1.,1.5),(1.5,2.),(2.,5.),(5.,10.)]:
            value,qerr=quad_vec(scaled_integrand,np.log(low),np.log(high),epsabs=1e-7,epsrel=a.epsrel,limit=500)
            rows.append(dict(teff_lo_K=lo,teff_hi_K=hi,r_low_pc=None,r_high_pc=None,relative_r_low=low,relative_r_high=high,thin=value[0],thick=value[1],halo=value[2],total=sum(value),quadrature_error=qerr,density_status=status,density_error=None if err is None else sum(value)*err/rho))
        pd.DataFrame(rows).to_csv(a.out/'native_radial_shell_counts.csv',index=False);print('Radial shell bin',lo,hi,flush=True)
    (a.out/'domain_status.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--calibration',type=Path,default=Path('/home/bahareh/results/ali_calibration_final'));p.add_argument('--cache',type=Path,default=Path('/home/bahareh/results/ali_fiveband_reconstruction'));p.add_argument('--maps',type=Path,default=Path('/home/bahareh/results/ali_final_counts/spatial_maps.npz'));p.add_argument('--out',type=Path,default=Path('/home/bahareh/results/ali_domain_validation'));p.add_argument('--cells',type=int,default=256);p.add_argument('--epsrel',type=float,default=3e-5);p.add_argument('--relative-floor-only',action='store_true');run(p.parse_args())
