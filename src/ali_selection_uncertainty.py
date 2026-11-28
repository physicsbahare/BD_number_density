"""Conditional finite-blank-CDF uncertainty via linear selection weights.

This is the sampling uncertainty of the empirical CDF, conditional on the
noise/PSF/SNR calibration and the explicit depth-limited radial domain. Shared
blank scenes across temperature bins are kept correlated. Spatially separate
blank apertures are approximated as independent; no model systematics included.
"""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from numpy.polynomial.legendre import leggauss
from numpy.polynomial.hermite import hermgauss
from src.ali_count_run import depth_groups
from src.ali_spatial_counts import density_components,signed_probability
from src.cosmosweb_counts import TEFF_DENSITY_2024


def run(a):
    a.out.parent.mkdir(parents=True,exist_ok=True)
    with np.load(a.maps) as d:
        use=d['five_area_sr']>0;values=d['properties'][use].astype(float);areas=d['five_area_sr'][use];tids=d['tile_id'][use];ftids=d['f814w_tile_id'][use];tiles=d['tiles']
    group=depth_groups(values,areas,tids*len(tiles)+ftids,list(range(10)),a.cells_per_tile);ng=group.max()+1
    area=np.bincount(group,weights=areas);v=np.column_stack([np.bincount(group,weights=areas*values[:,k])/area for k in range(20)])
    if a.preserved_nircam:
        for k in [0,1,2,4]:v[:,k+5]=v[:,k];v[:,k+10]=0.
    tid=np.rint(np.bincount(group,weights=areas*ftids)/area).astype(int)
    source_q=np.bincount(group,weights=areas*values[:,25])/area
    strata=[]
    for t in np.unique(tid):
        name=str(tiles[t]);r=json.loads((a.calibration/(name+'.json')).read_text())['F814W']
        with np.load(a.calibration/(name+'.npz')) as d:w=d['F814W_test_wht'][::2];z=d['F814W_test_z'][::2]
        edges=np.quantile(w,[.25,.5,.75]);which=np.searchsorted(edges,w)
        mapw=(r['catalog_noise_equalized_scale']/v[:,8])**2;depth=np.searchsorted(edges,mapw)
        for j in range(4):
            zsort=np.sort(z[which==j]);strata.append(dict(tile=name,depth=j,indices=np.flatnonzero((tid==t)&(depth==j)),z=zsort,weights=np.zeros(len(zsort))))
    table=pd.read_csv(Path(__file__).resolve().parents[1]/'data/elfowl_fiveband_bin_reach.csv')
    nodes,quadw=leggauss(a.nodes);nodes=(nodes+1)/2;quadw/=2;f25=10**((31.4-25)/2.5)
    hn,hw=hermgauss(a.shot_nodes);hw/=np.sqrt(np.pi)
    for _,row in table.iterrows():
        rho=TEFF_DENSITY_2024[(int(row.teff_lo_K),int(row.teff_hi_K))][0]
        amp=f25*np.array([row['d25_'+b+'_pc']**2 for b in ['F115W','F150W','F277W','F814W','F444W']])*v[:,15:20]
        reach=np.sqrt(amp[:,-1]/(5*v[:,4]))
        if a.radial_domain=='catalog_depth':reach=np.sqrt(amp[:,-1]/(5*v[:,9]))
        for s in strata:
            ids=s['indices'];z=s['z'];masses=np.zeros(len(z)+1)
            for node,qw in zip(nodes,quadw):
                r=node*reach[ids];flux=amp[ids]/r[:,None]**2+v[ids,10:15]
                p=np.ones(len(ids))
                for k in [0,1,2,4]:p*=signed_probability(flux[:,k],v[ids,k],v[ids,k+5],5 if k==4 else 2,k==4)
                threshold=(2*v[ids,8]-flux[:,3])/v[ids,3]
                shot=np.sqrt(np.maximum(amp[ids,3]/r**2*source_q[ids],0))/v[ids,3]
                weights=density_components(r,rho).sum(axis=0)*r*r*reach[ids]*qw*area[ids]*p
                for h,wh in zip(hn,hw):
                    masses+=np.bincount(np.searchsorted(z,threshold-h*np.sqrt(2)*shot,side='right'),weights=weights*wh,minlength=len(z)+1)
            # Each sample contributes whenever its residual is below threshold.
            s['weights']+=np.cumsum(masses[:0:-1])[::-1]
    total=sum(s['weights'].mean() for s in strata)
    variance=sum(np.var(s['weights'],ddof=1)/len(s['weights']) for s in strata)
    result=dict(method=__doc__,nodes=a.nodes,shot_nodes=a.shot_nodes,preserved_nircam=a.preserved_nircam,cells_per_tile=a.cells_per_tile,radial_domain=a.radial_domain,
        expected_count_ecdf_step_quadrature=float(total),conditional_cdf_sampling_standard_error=float(np.sqrt(variance)),
        strata=[dict(tile=s['tile'],depth_quartile=s['depth'],n=len(s['z']),count=float(s['weights'].mean()),variance_of_mean=float(np.var(s['weights'],ddof=1)/len(s['weights']))) for s in strata])
    a.out.write_text(json.dumps(result,indent=2)+'\n');print(total,np.sqrt(variance),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--maps',type=Path,required=True);p.add_argument('--calibration',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--nodes',type=int,default=256);p.add_argument('--shot-nodes',type=int,default=15);p.add_argument('--cells-per-tile',type=int,default=512);p.add_argument('--preserved-nircam',action='store_true');p.add_argument('--radial-domain',choices=['depth','catalog_depth'],default='depth');run(p.parse_args())
