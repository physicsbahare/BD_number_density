"""Execute depth-boundary sensitivities using completed spatial maps only."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from numpy.polynomial.hermite import hermgauss
from src.ali_count_run import depth_groups
from src.ali_spatial_counts import integrate_cells
from src.cosmosweb_counts import TEFF_DENSITY_2024


def run(a):
    a.out.mkdir(parents=True,exist_ok=True)
    with np.load(a.maps) as z:
        valid=z['five_area_sr']>0;v=z['properties'][valid].astype(float)
        area=z['five_area_sr'][valid];tid=z['tile_id'][valid];ftid=z['f814w_tile_id'][valid];tiles=z['tiles']
    group=depth_groups(v,area,tid*len(tiles)+ftid,list(range(10)),a.cells);ar=np.bincount(group,weights=area)
    means=np.column_stack([np.bincount(group,weights=area*v[:,k])/ar for k in range(20)])
    q=np.bincount(group,weights=area*v[:,25])/ar
    tg=np.rint(np.bincount(group,weights=area*ftid)/ar).astype(int);strata=[]
    for t in np.unique(tg):
        name=str(tiles[t]);report=json.loads((a.calibration/(name+'.json')).read_text())['F814W']
        with np.load(a.calibration/(name+'.npz')) as z:w=z['F814W_test_wht'][::2];res=z['F814W_test_z'][::2]
        edges=np.quantile(w,[.25,.5,.75]);which=np.searchsorted(edges,w)
        depth=np.searchsorted(edges,(report['catalog_noise_equalized_scale']/means[:,8])**2)
        for j in range(4):
            ix=np.flatnonzero((tg==t)&(depth==j));sample=np.sort(res[which==j])
            strata.append((ix,sample,np.arange(1,len(sample)+1)/len(sample)))
    hn,hw=hermgauss(a.shot_nodes);hw/=np.sqrt(np.pi)
    def cdf(z,shot=None):
        p=np.zeros(len(z))
        for ix,s,pr in strata:
            if shot is None:p[ix]=np.interp(z[ix],s,pr,left=0,right=1)
            else:p[ix]=sum(wh*np.interp(z[ix]-np.sqrt(2)*h*shot[ix],s,pr,left=0,right=1) for h,wh in zip(hn,hw))
        return p
    configs={
        'boundary4':dict(depth_sigma=4.),'boundary6':dict(depth_sigma=6.),
        'noise095':dict(noise_scale=.95),'noise105':dict(noise_scale=1.05),
        'threshold_minus05':dict(threshold_shift=-.5),'threshold_plus05':dict(threshold_shift=.5),
        'native_r50':dict(native_r50=True),'convergence':dict(),
        'correlation':dict(nircam_correlation=.224419),
        'native_correlation':dict(nircam_correlation=.1671),
    }
    table=pd.read_csv(Path(__file__).resolve().parents[1]/'data/elfowl_fiveband_bin_reach.csv');rows=[];reachrows=[]
    for scenario in a.scenarios:
        cfg=configs[scenario];noise=means[:,:5].copy();ref=means[:,5:10].copy();mu=means[:,10:15].copy();ee=means[:,15:20]
        noise[:,3]*=cfg.get('noise_scale',1);mu[:,3]-=cfg.get('threshold_shift',0)*ref[:,3]
        for _,row in table.iterrows():
            lo,hi=int(row.teff_lo_K),int(row.teff_hi_K);rho,err,status=TEFF_DENSITY_2024[(lo,hi)]
            d25=row[['d25_'+b+'_pc' for b in ['F115W','F150W','F277W','F814W','F444W']]].to_numpy(float)
            # Only the integration boundary changes; probability denominators
            # and empirical noise retain their calibrated values.
            limit=5*noise[:,4] if not cfg.get('native_r50') else 5*ref[:,4]-mu[:,4]
            depth=cfg.get('depth_sigma',5.)
            custom=limit/(5*noise[:,4]) if cfg.get('native_r50') else 1.
            val,qe=integrate_cells(d25,rho,noise,ref,ee,ar,background_mean=mu,
                empirical_cdf=cdf,f814w_source_coefficient=q,depth_sigma=depth,
                boundary_flux_scale=custom,nircam_correlation=cfg.get('nircam_correlation',0),epsrel=1e-6)
            rows.append(dict(scenario=scenario,teff_lo_K=lo,teff_hi_K=hi,thin=val[0],thick=val[1],halo=val[2],total=sum(val),
                measured_density_sigma=None if err is None else sum(val)*err/rho,density_status=status,quadrature_error=qe))
            reach=np.sqrt(10**((31.4-25)/2.5)*d25[4]**2*ee[:,4]/(depth*noise[:,4]*custom))
            order=np.argsort(reach);cum=np.cumsum(ar[order])/sum(ar)
            quant=np.interp([.16,.5,.84],cum,reach[order])
            reachrows.append(dict(scenario=scenario,teff_lo_K=lo,teff_hi_K=hi,r16_pc=quant[0],r50_pc=quant[1],r84_pc=quant[2]))
        pd.DataFrame(rows).to_csv(a.out/'sensitivity_bins.csv',index=False)
        pd.DataFrame(reachrows).to_csv(a.out/'distance_boundaries.csv',index=False)
        print(scenario,sum(r['total'] for r in rows if r['scenario']==scenario),flush=True)
    (a.out/'execution.json').write_text(json.dumps(dict(maps=str(a.maps),calibration=str(a.calibration),cells=a.cells,shot_nodes=a.shot_nodes,scenarios=a.scenarios,area_deg2=float(sum(ar)*(180/np.pi)**2)),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--maps',type=Path,default=Path('/home/bahareh/results/ali_final_counts/spatial_maps.npz'))
    p.add_argument('--calibration',type=Path,default=Path('/home/bahareh/results/ali_calibration_final'))
    p.add_argument('--out',type=Path,required=True);p.add_argument('--cells',type=int,default=512)
    p.add_argument('--shot-nodes',type=int,default=15)
    p.add_argument('--scenarios',nargs='+',default=['boundary4','boundary6','noise095','noise105','threshold_minus05','threshold_plus05','native_r50','native_correlation'])
    run(p.parse_args())
