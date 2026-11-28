"""Integrate reconstructed spatial depth maps using existing population inputs."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from src.ali_spatial_counts import integrate_cells
from src.cosmosweb_counts import TEFF_DENSITY_2024
from src.ali_depth_analysis import BANDS,psf_basis,encircled_map
from numpy.polynomial.hermite import hermgauss


def fill_missing_f814w(data,xy,tid,area5,max_distance):
    """Interpolate only inside real ACS coverage, retaining donor provenance."""
    good=(data[:,23]>0)&np.isfinite(data[:,3])&(data[:,3]>0)
    need=(area5>0)&~good
    provenance=tid.copy()
    if np.any(need):
        dist,j=cKDTree(xy[good]).query(xy[need])
        if np.max(dist)>max_distance:raise ValueError('F814W interpolation exceeds allowed distance')
        base=data[good][j]
        for col in [3,8,13,18,23]+([25] if data.shape[1]>25 else []):data[need,col]=base[:,col]
        provenance[need]=tid[good][j]
    return provenance,int(need.sum())


def depth_groups(values,weights,tile_ids,indices,cells_per_tile):
    """Deterministic adaptive quadrature in joint log-depth space.

    Every native-footprint cell keeps its exact solid angle. Representatives
    are area-weighted means; doubling leaves provides a convergence check.
    """
    group=np.empty(len(values),int);number=0
    for tile in np.unique(tile_ids):
        leaves=[np.flatnonzero(tile_ids==tile)]
        while len(leaves)<cells_per_tile:
            new=[];changed=False
            for ids in leaves:
                if len(ids)<4:new.append(ids);continue
                x=np.log(values[ids][:,indices]);axis=np.argmax(np.var(x,axis=0))
                order=ids[np.argsort(x[:,axis])];cum=np.cumsum(weights[order])
                split=int(np.clip(np.searchsorted(cum,cum[-1]/2),1,len(order)-1))
                new.extend([order[:split],order[split:]]);changed=True
            leaves=new
            if not changed:break
        for ids in leaves:group[ids]=number;number+=1
    return group


def run(args):
    args.out.mkdir(parents=True,exist_ok=True)
    files=sorted(args.calibration.glob('*.npz'))
    if len(files)!=20 and not args.allow_partial:raise ValueError('All 20 calibrated tiles are required')
    records=[];allxy=[];tileids=[];testz={};testweight={};snrscale={}
    for tid,f in enumerate(files):
        with np.load(f) as d:
            xy=d['global_xy'];allxy.append(xy);tileids.append(np.full(len(xy),tid))
            # sigmas, SNR denominators, means, and encircled energies, plus validity.
            report=json.loads(f.with_suffix('.json').read_text())['F814W']
            ph,pb=psf_basis(Path(report['psf_path']))
            moment2=encircled_map(ph,pb,d['grid_xy'],report['aperture_radius_arcsec'],kernel_power=2)
            wht=(report['catalog_noise_equalized_scale']/d['F814W_sigma_snr_njy'])**2
            source_coefficient=moment2/d['F814W_ee']*float(d['F814W_shot_coefficient'])/wht
            records.append(np.column_stack([d[b+'_'+suffix] for suffix in
                 ['sigma_noise_njy','sigma_snr_njy','background_njy','ee','valid'] for b in BANDS]+[source_coefficient]))
            if args.nircam_calibration is not None:
                with np.load(args.nircam_calibration/f.name) as legacy:
                    if not np.array_equal(legacy['global_xy'],xy):raise ValueError('NIRCam calibration grids differ')
                    for j,suffix in enumerate(['sigma_noise_njy','sigma_snr_njy','background_njy','ee','valid']):
                        for k in [0,1,2,4]:records[-1][:,j*5+k]=legacy[BANDS[k]+'_'+suffix]
            testz[tid]=np.sort(d['F814W_test_z'])
            # Disjoint halves of held-out blanks: even indices estimate the
            # CDF, odd indices independently validate it with injections.
            testweight[tid]=(np.array(d['F814W_test_wht'][::2]),np.array(d['F814W_test_z'][::2]))
            snrscale[tid]=json.loads(f.with_suffix('.json').read_text())['F814W']['catalog_noise_equalized_scale']
    xy=np.vstack(allxy);data=np.vstack(records);tid=np.concatenate(tileids)
    nircols=[i+5*j for j in range(4) for i in [0,1,2,4]]
    good=np.all(data[:,20:23]>0,axis=1)&(data[:,24]>0)&np.isfinite(data[:,nircols]).all(axis=1)
    # Choose the smallest F444W noise among duplicate overlapping tile samples.
    order=np.lexsort((data[:,4],xy[:,1],xy[:,0]));xy,data,tid,good=xy[order],data[order],tid[order],good[order]
    # Prefer valid common four-band samples before taking unique coordinates.
    order=np.lexsort((data[:,4],~good,xy[:,1],xy[:,0]));xy,data,tid,good=xy[order],data[order],tid[order],good[order]
    _,ix=np.unique(xy,axis=0,return_index=True);xy,data,tid,good=xy[ix],data[ix],tid[ix],good[ix]
    if args.coverage and args.coverage.exists():
        with np.load(args.coverage) as a:
            a4,a5=a['four_area_sr'],a['five_area_sr'];origin=a['global_origin_xy'];step=int(a['step'])
        if args.original_coverage is not None:
            with np.load(args.original_coverage) as a:
                if not np.array_equal(a['global_origin_xy'],origin) or int(a['step'])!=step:raise ValueError('Original and five-band coverage grids differ')
                a4=a['original_four_area_sr']
        if np.any(a5>a4+1e-20):raise ValueError('Five-band footprint exceeds original four-band coverage')
        yy,xx=np.nonzero(a4>0);target=np.column_stack([origin[0]+step//2+xx*step,origin[1]+step//2+yy*step])
        query=cKDTree(xy[good]);dist,j=query.query(target);good_indices=np.flatnonzero(good)[j]
        if np.max(dist)>5*step:raise ValueError('Depth interpolation spans more than five cells')
        boundary_fraction=float(np.sum(a4[yy,xx][dist>0])/np.sum(a4))
        data,tid=data[good_indices],tid[good_indices];xy=target
        area4,area5=a4[yy,xx],a5[yy,xx]
    else:
        step=100
        data,tid,xy=data[good],tid[good],xy[good]
        rad=np.deg2rad(.03/3600);area4=(step*rad)**2/(1+np.sum((xy*rad)**2,axis=1))**1.5
        area5=area4*(data[:,23]>0);boundary_fraction=None
    # Missing F814W values may only be interpolated INSIDE a native valid five-band cell.
    ftid,n_f814w_interpolated=fill_missing_f814w(data,xy,tid,area5,5*step)
    d25=pd.read_csv(Path(__file__).resolve().parents[1]/'data/elfowl_fiveband_bin_reach.csv')
    np.savez_compressed(args.out/'spatial_maps.npz',global_xy=xy,properties=data.astype('float32'),
        four_area_sr=area4,five_area_sr=area5,tile_id=tid,f814w_tile_id=ftid,tiles=np.array([f.stem for f in files]))
    output=[]
    for selection,area,indices in [('four_original',area4,[0,1,2,4]),('four_matched',area5,[0,1,2,4]),('five_gaussian',area5,[0,1,2,3,4]),('five_empirical',area5,[0,1,2,3,4])]:
        if selection not in args.selections:continue
        use=area>0;d=data[use];t=tid[use];a=area[use]
        features=indices+[i+5 for i in indices] if selection=='four_original' else list(range(10))
        ft=ftid[use]
        group=depth_groups(d,a,t*len(files)+ft,features,args.cells_per_tile);ng=group.max()+1
        ar=np.bincount(group,weights=a,minlength=ng)
        means=np.column_stack([np.bincount(group,weights=a*d[:,k],minlength=ng)/ar for k in range(20)])
        tg=np.rint(np.bincount(group,weights=a*ft,minlength=ng)/ar).astype(int)
        noise=means[:,indices]
        ref=noise.copy() if args.snr_mode=='empirical' else means[:,np.array(indices)+5]
        ee=means[:,np.array(indices)+15]
        mu=np.zeros_like(noise) if args.zero_background else means[:,np.array(indices)+10]
        if args.snr_mode=='preserved_nircam':
            for k,i in enumerate(indices):
                if i!=3:ref[:,k]=noise[:,k];mu[:,k]=0.
        if 'F814W' in [BANDS[i] for i in indices]:
            ik=indices.index(3)
            noise[:,ik]*=args.f814w_noise_scale
            ref[:,ik]*=args.f814w_snr_scale
            ee[:,ik]*=args.f814w_flux_scale
            mu[:,ik]-=args.f814w_threshold_shift*ref[:,ik]
        cdf_groups=[]
        for tile in np.unique(tg):
            weights,zblank=testweight[tile];edges=np.quantile(weights,[.25,.5,.75]);blankstratum=np.searchsorted(edges,weights)
            mapweight=(snrscale[tile]/means[:,8])**2;stratum=np.searchsorted(edges,mapweight)
            for depth in range(4):
                sub=(tg==tile)&(stratum==depth);sample=np.sort(zblank[blankstratum==depth])
                if np.any(sub):cdf_groups.append((np.flatnonzero(sub),sample,np.arange(1,len(sample)+1)/len(sample)))
        hn,hw=hermgauss(args.shot_nodes);hw/=np.sqrt(np.pi)
        def cdf(z,shot_sigma=None):
            p=np.zeros(len(z))
            for sub,sample,probability in cdf_groups:
                if shot_sigma is None:p[sub]=np.interp(z[sub],sample,probability,left=0.,right=1.)
                else:p[sub]=sum(w*np.interp(z[sub]-h*np.sqrt(2)*shot_sigma[sub],sample,probability,left=0.,right=1.) for h,w in zip(hn,hw))
            return p
        shot_coefficient=np.bincount(group,weights=a*d[:,25])/ar if len(indices)==5 and not args.no_source_noise else None
        print(selection,'cells',len(a),'grouped',ng,flush=True)
        for _,row in d25.iterrows():
            lo,hi=int(row.teff_lo_K),int(row.teff_hi_K);rho,err,status=TEFF_DENSITY_2024[(lo,hi)]
            reach=[row['d25_'+BANDS[k]+'_pc'] for k in indices]
            val,qerr=integrate_cells(reach,rho,noise,ref,ee,ar,rmax=args.rmax,
                background_mean=mu,epsrel=args.epsrel,empirical_cdf=cdf if selection=='five_empirical' else None,radial_domain=args.radial_domain,f814w_source_coefficient=shot_coefficient,nircam_correlation=args.nircam_correlation)
            output.append(dict(selection=selection,teff_lo_K=lo,teff_hi_K=hi,
                 thin=float(val[0]),thick=float(val[1]),halo=float(val[2]),total=float(sum(val)),
                 density_uncertainty=None if err is None else float(sum(val)*err/rho),
                 density_status=status,quadrature_error_bound=qerr))
        print(selection,'total',sum(r['total'] for r in output if r['selection']==selection),flush=True)
    pd.DataFrame(output).to_csv(args.out/'temperature_component_counts.csv',index=False)
    summary=dict(status='COMPLETE_NUMERICAL_RUN' if len(files)==20 and args.coverage else 'PARTIAL_TILE_PREVIEW',
        tiles=[f.stem for f in files],snr_mode=args.snr_mode,cells_per_tile=args.cells_per_tile,rmax_pc=args.rmax,radial_domain=args.radial_domain,
        f814w_noise_scale=args.f814w_noise_scale,f814w_snr_scale=args.f814w_snr_scale,f814w_threshold_shift=args.f814w_threshold_shift,
        four_area_deg2=float(area4.sum()*(180/np.pi)**2),five_area_deg2=float(area5.sum()*(180/np.pi)**2),
        boundary_interpolated_area_fraction=boundary_fraction,results={})
    summary['nircam_calibration']=str(args.nircam_calibration) if args.nircam_calibration else str(args.calibration)
    summary['f814w_interpolated_cells']=n_f814w_interpolated
    summary['f814w_source_noise']=not args.no_source_noise
    summary['source_noise_quadrature_nodes']=args.shot_nodes
    summary['f814w_flux_scale']=args.f814w_flux_scale
    summary['nircam_F277W_F444W_correlation_sensitivity']=args.nircam_correlation
    for key in sorted(set(r['selection'] for r in output)):
        rows=[r for r in output if r['selection']==key];n=sum(r['total'] for r in rows)
        summary['results'][key]=dict(total=n,thin=sum(r['thin'] for r in rows),thick=sum(r['thick'] for r in rows),halo=sum(r['halo'] for r in rows),
            measured_density_uncertainty=float(np.sqrt(sum((r['density_uncertainty'] or 0)**2 for r in rows))),
            incomplete_density_bin_contribution=sum(r['total'] for r in rows if r['density_status']=='lower_limit'),
            quadrature_error_bound=sum(r['quadrature_error_bound'] for r in rows))
    (args.out/'count_summary.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--calibration',type=Path,required=True)
    p.add_argument('--coverage',type=Path)
    p.add_argument('--original-coverage',type=Path)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--snr-mode',choices=['catalog','empirical','preserved_nircam'],default='catalog')
    p.add_argument('--nircam-calibration',type=Path)
    p.add_argument('--zero-background',action='store_true')
    p.add_argument('--allow-partial',action='store_true')
    p.add_argument('--log-step',type=float,default=.05)
    p.add_argument('--cells-per-tile',type=int,default=512)
    p.add_argument('--rmax',type=float,default=50000.)
    p.add_argument('--radial-domain',choices=['depth','catalog_depth','extended'],default='depth')
    p.add_argument('--epsrel',type=float,default=1e-6)
    p.add_argument('--f814w-noise-scale',type=float,default=1.)
    p.add_argument('--f814w-snr-scale',type=float,default=1.)
    p.add_argument('--f814w-flux-scale',type=float,default=1.)
    p.add_argument('--f814w-threshold-shift',type=float,default=0.)
    p.add_argument('--no-source-noise',action='store_true')
    p.add_argument('--shot-nodes',type=int,default=15)
    p.add_argument('--nircam-correlation',type=float,default=0.)
    p.add_argument('--selections',nargs='+',default=['four_original','four_matched','five_gaussian','five_empirical'])
    run(p.parse_args())
