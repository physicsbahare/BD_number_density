"""Bound fractional-phase footprint sensitivity by one extra pixel erosion."""
import argparse,json
from pathlib import Path
import numpy as np


def run(a):
    meta=json.loads((a.coverage/'geometry.json').read_text());ny,nx=meta['native_shape_yx'];nb=(nx+7)//8
    maps=[np.memmap(a.coverage/(b+'.bits'),mode='r',dtype='uint8',shape=(ny,nb)) for b in ['F115W','F150W','F277W','F444W','F814W']]
    lut=np.array([i.bit_count() for i in range(256)],dtype='uint8');n0=0;n1=0
    for y in range(0,ny,256):
        start=max(0,y-1);end=min(ny,y+257);joint=maps[0][start:end].copy()
        for m in maps[1:]:joint&=m[start:end]
        left=(joint<<1)|(np.roll(joint,-1,axis=1)>>7);right=(joint>>1)|(np.roll(joint,1,axis=1)<<7)
        row=joint&left&right;row[:,0]&=127;row[:,-1]&=254
        eroded=row.copy();eroded[1:]&=row[:-1];eroded[:-1]&=row[1:]
        if start==0:eroded[0]=0
        if end==ny:eroded[-1]=0
        lo=y-start;hi=min(y+256,ny)-start;n0+=int(lut[joint[lo:hi]].sum());n1+=int(lut[eroded[lo:hi]].sum())
    result=dict(five_native_pixels=n0,extra_one_pixel_eroded_pixels=n1,lost_pixels=n0-n1,
        area_upper_bound_deg2=(n0-n1)*(.03/3600)**2,
        relative_area_upper_bound=(n0-n1)/n0,
        definition='Conservative phase-support sensitivity proxy: one extra native-pixel erosion of the global five-band valid-centered-aperture intersection. Pixel area at tangent point bounds the lost solid angle.')
    a.out.write_text(json.dumps(result,indent=2)+'\n');print(result)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--coverage',type=Path,required=True);p.add_argument('--out',type=Path,required=True);run(p.parse_args())
