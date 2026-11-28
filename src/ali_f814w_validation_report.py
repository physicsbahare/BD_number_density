"""Aggregate executed F814W trials; save compact tables and scientific figures."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import ndtr
from src.f814w_targeted_recovery import wilson


def run(a):
    a.out.mkdir(parents=True,exist_ok=True)
    s=json.loads((a.recovery/'recovery_summary.json').read_text())
    if len(s['tiles'])!=20:raise ValueError('All twenty recovery fields required')
    pd.DataFrame(s['records']).to_csv(a.out/'injection_by_tile.csv',index=False)
    stars=pd.DataFrame(s['star_controls']);stars.to_csv(a.out/'psf_controls.csv',index=False)
    rows=[];curves=[];bydepth=[];trialcount=0;nb=0
    for tile in s['tiles']:
        with np.load(a.recovery/(tile+'_recovery.npz')) as z:
            levels=z['levels'];snr=z['measured_snr'];gauss=z['gaussian_probability']
            expected=z['expected_aperture_njy'];rec=z['recovered_aperture_njy'];sig=z['blank_sigma_njy'];mu=z['blank_mean_njy'];flux=z['injected_total_njy']
            if not curves: allp=[[] for _ in levels];allg=[[] for _ in levels];bias=[[] for _ in levels];inj=[[] for _ in levels]
            curves.append(tile);nb+=len(sig);trialcount+=snr.size
            edges=np.quantile(sig,[.25,.5,.75]);depth=np.searchsorted(edges,sig)
            for j,level in enumerate(levels):
                allp[j].extend((snr[j]<2).tolist());allg[j].extend(gauss[j].tolist());bias[j].extend((rec[j]-mu-expected[j]).tolist());inj[j].extend(flux[j].tolist())
                for k in range(4):
                    use=depth==k;n=int(use.sum());npass=int((snr[j,use]<2).sum())
                    bydepth.append(dict(tile=tile,noise_quartile=k,true_aperture_snr=float(level),n=n,n_pass=npass,gaussian_probability=float(gauss[j,use].mean())))
    for j,level in enumerate(levels):
        p=np.array(allp[j]);g=np.array(allg[j]);n=len(p);v=p.mean();se=np.sqrt(v*(1-v)/n)
        rows.append(dict(true_aperture_snr=float(level),n=n,n_pass=int(p.sum()),empirical_probability=v,
            standard_error=se,wilson95=wilson(p.sum(),n),gaussian_probability=g.mean(),difference=v-g.mean(),
            difference_over_binomial_se=float((v-g.mean())/se) if se else None,
            aperture_bias_njy=float(np.mean(bias[j])),injected_flux_njy_p16_p50_p84=np.percentile(inj[j],[16,50,84]).tolist(),
            convergence={str(k):float(p[:k].mean()) for k in [1000,10000,40000,n] if k<=n}))
    summary=dict(status='EXECUTED_TARGETED_APERTURE_VALIDATION',n_tiles=20,distinct_background_scenes=nb,
        injection_trials=trialcount,probabilities=rows,psf_controls_n=len(stars),
        psf_core_observed_over_model_p16_p50_p84=np.percentile(stars.ratio,[16,50,84]).tolist(),
        limitations=['Scenes reused across flux levels: trials are not all independent; distinct random background apertures are treated as approximately independent.',
        'Backgrounds exclude catalog neighbors; no detection, deblending, morphology, visual or quality-flag completeness.',
        'PSF controls validate the core ratio at 0.10/0.18 arcsec, not the absolute full-PSF aperture correction.',
        'Binomial uncertainties condition on estimated noise, registration, catalog-SNR mapping and atmosphere fluxes.'])
    (a.out/'injection_validation.json').write_text(json.dumps(summary,indent=2)+'\n')
    pd.DataFrame(rows).to_csv(a.out/'injection_probabilities.csv',index=False)
    pd.DataFrame(bydepth).to_csv(a.out/'injection_depth_strata.csv',index=False)
    pop=json.loads((a.population/'population_validation_summary.json').read_text())
    pd.DataFrame(pop['records']).to_csv(a.out/'population_flux_validation_by_tile.csv',index=False)
    pr=pd.DataFrame(pop['records']);combined=[]
    for lo,v in pr.groupby('teff_lo_K'):
        n=int(v.n.sum());k=int(v.n_pass.sum());combined.append(dict(teff_lo_K=int(lo),n=n,n_pass=k,empirical_probability=k/n,gaussian_probability=float(np.average(v.gaussian_probability,weights=v.n)),wilson95=wilson(k,n)))
    (a.out/'population_flux_validation.json').write_text(json.dumps(dict(records=combined,n_trials=int(pr.n.sum())),indent=2)+'\n')
    plt.rcParams.update({'font.size':10,'svg.fonttype':'none'})
    fig,ax=plt.subplots(1,2,figsize=(10,4),layout='constrained');x=np.array(levels);p=np.array([r['empirical_probability'] for r in rows]);g=np.array([r['gaussian_probability'] for r in rows]);se=np.array([r['standard_error'] for r in rows])
    ax[0].errorbar(x,p,yerr=1.96*se,fmt='o',label='Injected sources (95% binomial CI)');ax[0].plot(x,g,'--',label='Gaussian prediction');ax[0].set(xlabel='Injected aperture flux / empirical noise',ylabel='P(signed catalog SNR < 2)',ylim=(-.02,1.02));ax[0].legend(fontsize=8)
    ax[1].errorbar(x,p-g,yerr=1.96*se,fmt='o');ax[1].axhline(0,color='k',lw=.7);ax[1].set(xlabel='Injected aperture flux / empirical noise',ylabel='Empirical − Gaussian probability')
    for ext in ['svg','png']:fig.savefig(a.out/('f814w_injection_probability.'+ext),dpi=160)
    plt.close(fig)
    fig,ax=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    ax[0].plot(levels,[np.mean(v) for v in bias],'o-');ax[0].axhline(0,color='k',lw=.7);ax[0].set(xlabel='Injected aperture flux / empirical noise',ylabel='Recovered − background − input (nJy)')
    ax[1].hist(stars.ratio,bins=18);ax[1].axvline(1,color='k');ax[1].set(xlabel='Observed / PSF core aperture ratio',ylabel='Unmasked stellar controls')
    for ext in ['svg','png']:fig.savefig(a.out/('f814w_flux_and_psf_validation.'+ext),dpi=160)
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(6,4),layout='constrained')
    for j in [0,5,6]:
        p=np.array(allp[j]);ns=np.unique(np.geomspace(250,len(p),30).astype(int));cs=np.cumsum(p)
        ax.plot(ns,cs[ns-1]/ns,label=f'True aperture SNR {levels[j]:g}')
    ax.set(xscale='log',xlabel='Distinct background scenes',ylabel='Cumulative P(catalog SNR < 2)');ax.legend()
    for ext in ['svg','png']:fig.savefig(a.out/('f814w_recovery_convergence.'+ext),dpi=160)
    plt.close(fig)
    print(json.dumps({k:v for k,v in summary.items() if k not in ['probabilities','limitations']},indent=2))
    print(pd.DataFrame(rows)[['true_aperture_snr','n','empirical_probability','gaussian_probability','difference_over_binomial_se']].to_string(index=False))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--recovery',type=Path,required=True);p.add_argument('--population',type=Path,required=True);p.add_argument('--out',type=Path,required=True);run(p.parse_args())
