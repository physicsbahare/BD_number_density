"""Plot actual spatial sensitivity, coverage and conditional population counts."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def run(a):
    a.out.mkdir(parents=True,exist_ok=True);plt.rcParams.update({'font.size':10,'svg.fonttype':'none'})
    def save(fig,name):
        for ext in ['svg','png']:fig.savefig(a.out/(name+'.'+ext),dpi=80 if ext=='svg' else 160)
        plt.close(fig)
    with np.load(a.counts/'spatial_maps.npz') as d:
        xy=d['global_xy']*.03/60;v=d['properties'];area=d['five_area_sr'];a4=d['four_area_sr']
    keep=area>0
    fig,axs=plt.subplots(1,2,figsize=(11,4.8),layout='constrained')
    for ax,col,title in zip(axs,[3,8],['Empirical 2σ aperture sensitivity','Nominal catalog SNR=2 flux scale']):
        mag=31.4-2.5*np.log10(2*v[keep,col]/v[keep,18]);lo,hi=np.percentile(mag,[1,99])
        points=ax.scatter(xy[keep,0],xy[keep,1],c=mag,s=.35,cmap='viridis',vmin=lo,vmax=hi,rasterized=True)
        ax.set(aspect='equal',xlabel='TAN x (arcmin)',ylabel='TAN y (arcmin)',title=title);fig.colorbar(points,ax=ax,label='F814W total point-source AB magnitude')
    save(fig,'f814w_depth_maps')
    with np.load(a.coverage) as d:four=d['four_area_sr'];five=d['five_area_sr'];origin=d['global_origin_xy'];step=int(d['step'])
    mask=np.full(four.shape,np.nan);mask[four>0]=0;mask[five>0]=1
    extent=[origin[0]*.03/60,(origin[0]+four.shape[1]*step)*.03/60,origin[1]*.03/60,(origin[1]+four.shape[0]*step)*.03/60]
    fig,ax=plt.subplots(figsize=(6.5,6),layout='constrained');im=ax.imshow(mask,origin='lower',extent=extent,vmin=0,vmax=1,cmap='coolwarm',interpolation='nearest');ax.set(xlabel='TAN x (arcmin)',ylabel='TAN y (arcmin)',title='Native valid-aperture footprint intersection');c=fig.colorbar(im,ax=ax,ticks=[0,1]);c.ax.set_yticklabels(['NIRCam only','Five bands']);save(fig,'f814w_coverage')
    table=pd.read_csv(a.counts/'temperature_component_counts.csv');selection='five_empirical';d=table[table.selection==selection]
    fig,axs=plt.subplots(1,2,figsize=(11,4.3),layout='constrained');x=(d.teff_lo_K+d.teff_hi_K)/2;bottom=np.zeros(len(d))
    for key,color in zip(['thin','thick','halo'],['#2878b5','#f39c34','#7b52a1']):axs[0].bar(x,d[key],bottom=bottom,width=135,label=key,color=color);bottom+=d[key].to_numpy()
    axs[0].set(xlabel='Effective temperature (K)',ylabel='Conditional expected count',title='Empirical F814W selection');axs[0].legend()
    for key,label in [('four_original','4 bands, reconstructed full area'),('four_matched','4 bands, five-band area'),('five_gaussian','5 bands, Gaussian ACS'),('five_empirical','5 bands, empirical ACS')]:
        t=table[table.selection==key];axs[1].plot((t.teff_lo_K+t.teff_hi_K)/2,t.total,'o-',ms=3,label=label)
    axs[1].set(xlabel='Effective temperature (K)',ylabel='Conditional expected count');axs[1].legend(fontsize=8);save(fig,'temperature_component_counts')
    depth=dict(empirical_F814W_2sigma_AB_p16_p50_p84=np.percentile(31.4-2.5*np.log10(2*v[keep,3]/v[keep,18]),[16,50,84]).tolist(),nominal_catalog_2SNR_AB_p16_p50_p84=np.percentile(31.4-2.5*np.log10(2*v[keep,8]/v[keep,18]),[16,50,84]).tolist(),note='Spatial-cell percentiles, not area-weighted. Nominal catalog threshold omits the fitted background offset; it is an estimator-scale diagnostic.')
    (a.out/'depth_diagnostics.json').write_text(json.dumps(depth,indent=2)+'\n');print(depth)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--counts',type=Path,required=True);p.add_argument('--coverage',type=Path,required=True);p.add_argument('--out',type=Path,required=True);run(p.parse_args())
