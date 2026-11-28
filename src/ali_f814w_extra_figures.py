"""Noise-distribution and population-flux diagnostics from executed caches."""
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import ndtr,ndtri


def run(a):
    strata=[[] for _ in range(4)]
    for p in sorted(a.calibration.glob('*.npz')):
        with np.load(p) as d:
            w=d['F814W_test_wht'];z=d['F814W_test_z'];k=np.searchsorted(np.quantile(w,[.25,.5,.75]),w)
            for i in range(4):strata[i].extend(z[k==i].tolist())
    plt.rcParams.update({'font.size':10,'svg.fonttype':'none'})
    fig,axs=plt.subplots(1,2,figsize=(10,4),layout='constrained');q=np.linspace(.005,.995,150)
    for i,s in enumerate(strata):
        z=np.sort(s);axs[0].plot(z,np.arange(1,len(z)+1)/len(z),label=f'WHT quartile {i+1}')
        axs[1].plot(ndtri(q),np.quantile(z,q),label=f'Quartile {i+1}')
    x=np.linspace(-4,6,400);axs[0].plot(x,ndtr(x),'k--',label='Unit Gaussian');axs[0].set(xlim=(-4,6),xlabel='Blank aperture residual / empirical noise',ylabel='Cumulative probability');axs[0].legend(fontsize=8)
    axs[1].plot([-3,3],[-3,3],'k--');axs[1].set(xlabel='Unit-Gaussian quantile',ylabel='Measured residual quantile')
    for ext in ['png','svg']:fig.savefig(a.out/('f814w_noise_distribution.'+ext),dpi=160)
    plt.close(fig)
    p=json.loads((a.out/'population_flux_validation.json').read_text());rows=p['records'];x=np.array([r['teff_lo_K']+75 for r in rows]);y=np.array([r['empirical_probability'] for r in rows]);lo=np.array([r['wilson95'][0] for r in rows]);hi=np.array([r['wilson95'][1] for r in rows])
    fig,ax=plt.subplots(figsize=(7,4.2),layout='constrained');ax.errorbar(x,y,yerr=[y-lo,hi-y],fmt='o',label='Real-background injections (95% CI)');ax.plot(x,[r['gaussian_probability'] for r in rows],'--',label='Gaussian prediction');ax.set(xlabel='Effective temperature (K)',ylabel='Fraction with signed catalog SNR < 2',title='Representative population fluxes; balanced test, not count weighting');ax.legend(fontsize=8)
    for ext in ['png','svg']:fig.savefig(a.out/('f814w_population_flux_validation.'+ext),dpi=160)
    plt.close(fig)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--calibration',type=Path,required=True);p.add_argument('--out',type=Path,required=True);run(p.parse_args())
