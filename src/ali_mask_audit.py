"""Inspect existing JWST star masks without applying an unrequested flag cut."""
import argparse,json,tarfile
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


def run(a):
    meta=json.loads((a.coverage/'geometry.json').read_text());origin=np.array(meta['global_origin_xy']);step=meta['step']
    with np.load(a.coverage/'area_maps.npz') as d:area=d['five_area_sr']
    yy,xx=np.nonzero(area>0);globalxy=origin+step//2+np.column_stack([xx,yy])*step
    masked=np.zeros(len(xx),bool);rows=[]
    with tarfile.open(a.archive,'r|gz') as tar:
        for member in tar:
            if not member.isfile() or not member.name.endswith('.fits'):continue
            tile=member.name.split('_')[-1].split('.')[0];f=tar.extractfile(member);h=fits.Header.fromfile(f)
            if h['BITPIX']!=8:raise ValueError('Unexpected mask dtype')
            rec=next(v for v in meta['inputs'] if v['band']=='F444W' and v['tile']==tile)
            offset=np.array(rec['global_offset_xy']);local=globalxy-offset;sx=h['NAXIS1'];sy=h['NAXIS2']
            inside=np.flatnonzero((local[:,0]>=0)&(local[:,0]<sx)&(local[:,1]>=0)&(local[:,1]<sy));counts=np.zeros(256,dtype=np.int64)
            for y in range(0,sy,128):
                n=min(128,sy-y);raw=f.read(n*sx)
                if len(raw)!=n*sx:raise IOError('Truncated star mask')
                z=np.frombuffer(raw,'uint8').reshape(n,sx);counts+=np.bincount(z.ravel(),minlength=256)
                ix=inside[(local[inside,1]>=y)&(local[inside,1]<y+n)]
                masked[ix]|=z[local[ix,1]-y,local[ix,0]]!=0
            rows.append(dict(tile=tile,header=str(h),pixel_value_counts={str(i):int(v) for i,v in enumerate(counts) if v}))
            print('Audited star mask',tile,flush=True)
    covered=area[yy,xx];result=dict(archive=str(a.archive),tiles=rows,nonzero_mask_center_area_deg2=float(covered[masked].sum()*(180/np.pi)**2),
        five_band_area_without_nonzero_mask_centers_deg2=float(covered[~masked].sum()*(180/np.pi)**2),
        definition='3 arcsec-cell center diagnostic of nonzero JWST star mask. No quality or bright-star cut is silently added to the five specified photometric criteria; this is not an exact native masked-footprint calculation.')
    a.out.write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--archive',type=Path,required=True);p.add_argument('--coverage',type=Path,required=True);p.add_argument('--out',type=Path,required=True);run(p.parse_args())
