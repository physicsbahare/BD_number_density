"""Verify original atmosphere normalization and regenerate thirteen-bin reaches."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd


def run(a):
    data=pd.read_csv(a.models,sep='\t');rows=[]
    for lo in range(300,2250,150):
        subset=data[(data.Teff>=lo)&(data.Teff<lo+150)]
        if len(subset)==0:raise ValueError('Missing atmosphere bin')
        row=dict(teff_lo_K=lo,teff_hi_K=lo+150,n_models=len(subset))
        for band in ['F115W','F150W','F277W','F444W','F814W']:
            row['d25_'+band+'_pc']=float(np.median(np.sqrt(subset['f_'+band])*subset.Radius_pc))
        rows.append(row)
    d=pd.DataFrame(rows);old=pd.read_csv(a.archived);checks=[]
    for i,row in d.iterrows():
        for b in ['F115W','F444W']:
            difference=abs(row['d25_'+b+'_pc']-old.iloc[i]['d25_'+b+'_pc']);checks.append(difference<=.0051)
    if not all(checks):raise ValueError('Original archived reaches not reproduced')
    d.to_csv(a.out,index=False);print(json.dumps(dict(n_models=len(data),n_archived_checks=len(checks),all_within_rounding=True)))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['models','archived','out']:p.add_argument('--'+n,type=Path,required=True)
    run(p.parse_args())
