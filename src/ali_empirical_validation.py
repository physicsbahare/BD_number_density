"""Independent odd-scene validation of the even-scene empirical ACS CDF."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import ndtri
from numpy.polynomial.hermite import hermgauss
from src.f814w_targeted_recovery import wilson


def run(a):
    rows=[]
    hn,hw=hermgauss(21);hw/=np.sqrt(np.pi)
    for f in sorted(a.calibration.glob('*.npz')):
        with np.load(f) as z:
            w=z['F814W_test_wht'];res=z['F814W_test_z'];den=z['F814W_test_snr_denom_njy'];off=z['F814W_test_snr_offset_njy'];mu=z['F814W_test_background_njy'];sig=z['F814W_test_sigma_njy']
        trainw=w[::2];trainz=res[::2];edges=np.quantile(trainw,[.25,.5,.75]);train=np.searchsorted(edges,trainw);test=np.searchsorted(edges,w[1::2])
        with np.load(a.recovery/(f.stem+'_recovery.npz')) as r:
            levels=r['levels'];expected=r['expected_aperture_njy'][:,1::2];passed=r['measured_snr'][:,1::2]<2;gauss=r['gaussian_probability'][:,1::2]
        for j,level in enumerate(levels):
            numerator=2*den[1::2]+off[1::2]-mu[1::2]-expected[j]
            threshold=numerator/sig[1::2];pred=np.zeros(len(threshold))
            normal_z=ndtri(gauss[j]);ok=np.isfinite(normal_z)&(np.abs(normal_z)>1e-8)
            combined_sd=sig[1::2].copy();combined_sd[ok]=np.abs(numerator[ok]/normal_z[ok])
            shot_sd=np.sqrt(np.maximum(combined_sd**2-sig[1::2]**2,0))/sig[1::2]
            for k in range(4):
                samples=np.sort(trainz[train==k]);use=test==k
                probs=np.arange(1,len(samples)+1)/len(samples)
                pred[use]=sum(weight*np.interp(threshold[use]-node*np.sqrt(2)*shot_sd[use],samples,probs,left=0,right=1) for node,weight in zip(hn,hw))
            rows.append(dict(tile=f.stem,level=float(level),n=len(pred),n_pass=int(passed[j].sum()),empirical_cdf_probability=float(pred.mean()),gaussian_probability=float(gauss[j].mean())))
    d=pd.DataFrame(rows);d.to_csv(a.out/'independent_empirical_validation_by_tile.csv',index=False);combined=[]
    for level,v in d.groupby('level'):
        n=int(v.n.sum());k=int(v.n_pass.sum());p=k/n;pred=float(np.average(v.empirical_cdf_probability,weights=v.n));se=np.sqrt(p*(1-p)/n)
        combined.append(dict(level=float(level),n=n,n_pass=k,recovered_probability=p,empirical_cdf_probability=pred,gaussian_probability=float(np.average(v.gaussian_probability,weights=v.n)),wilson95=wilson(k,n),difference=p-pred,difference_over_binomial_se=(p-pred)/se if se else None))
    (a.out/'independent_empirical_validation.json').write_text(json.dumps(dict(records=combined,
        definition='Noise parameters from original even blank apertures. ECDF from even-indexed held-out scenes; validation uses disjoint odd-indexed held-out scenes. Four WHT strata per tile. ECDF convolved with source photon noise using 21-node Gaussian quadrature; the exact injected variance is recovered from the saved Gaussian probability.',
        limitations='Binomial SE excludes finite ECDF-estimation error; no detection or catalog morphology completeness.'),indent=2)+'\n')
    print(pd.DataFrame(combined)[['level','n','recovered_probability','empirical_cdf_probability','difference_over_binomial_se']].to_string(index=False))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['calibration','recovery','out']:p.add_argument('--'+name,type=Path,required=True)
    run(p.parse_args())
