"""Estimate JWST star-mask sky union from user-supplied COSMOS2025 30mas FITS.
No WHT download or MAST queries. Requires only numpy.
Use --sample-out points.csv --samples 10000 to generate reproducible weighted
sky positions for the epoch-geometry intersection script.
Usage: python -m pm_footprint.jwst_mask_union_repro --archive FILE.tar.gz
0.3"/0.6" results are raster approximations, not a calibrated dropout effective area.
"""
import argparse, tarfile, numpy as np, time, json, hashlib
from .mask_wcs import (pixel_to_lambert, cd_from_header,
                       tan_pixel_area_deg2, legacy_column_scaling_error_arcsec)
from pathlib import Path
parser=argparse.ArgumentParser()
parser.add_argument('--archive',required=True)
parser.add_argument('--output',default='jwst_mask_union_conservative.json')
parser.add_argument('--sample-out',default=None,help='Optional output CSV of weighted masked-sky Lambert points (radians)')
parser.add_argument('--samples',type=int,default=10000)
parser.add_argument('--seed',type=int,default=20261031)
parser.add_argument('--shift-x-cells',type=float,default=0.0,help='Shift grid origin in x by this fraction of one cell')
parser.add_argument('--shift-y-cells',type=float,default=0.0,help='Shift grid origin in y by this fraction of one cell')
args=parser.parse_args()
F=args.archive;R=np.pi/180
archive_hasher=hashlib.sha256()
with open(F,'rb') as archive_stream:
 for block in iter(lambda: archive_stream.read(1024*1024), b''):
  archive_hasher.update(block)
archive_sha256=archive_hasher.hexdigest()
r0=150.12*R; d0=2.2*R; s0=np.sin(d0);c0=np.cos(d0)
GRIDS=[0.3,0.6]; ORIGIN=-0.7; SIZE=1.4
T=time.time()
def hdr(f):
 h={}
 while True:
  b=f.read(2880)
  if len(b)!=2880:raise RuntimeError('truncated header')
  end=False
  for p in range(0,2880,80):
   s=b[p:p+80].decode('ascii','replace')
   if s[:8].strip()=='END':end=True;break
   if s[8:10]=='= ':
    v=s[10:].split('/',1)[0].strip().strip("'")
    try:v=float(v.replace('D','E')) if any(x in v for x in '.ED') else int(v)
    except:pass
    h[s[:8].strip()]=v
  if end:return h

def sky(x,y,h):
 # FITS convention: CD[i,j] = CDELT[i] * PC[i,j] (row scaling).
 # The former column-scaled implementation misplaced rotated sky pixels.
 return pixel_to_lambert(x,y,h)


def depos(xe,ye,frac,g):
 gs=g/3600.;N=round(SIZE/gs)
 fx=(xe-ORIGIN)/gs-0.5-args.shift_x_cells;fy=(ye-ORIGIN)/gs-0.5-args.shift_y_cells
 x=np.floor(fx).astype(np.int32);y=np.floor(fy).astype(np.int32)
 a=fx-x;b=fy-y
 assert (x>=0).all() and (x+1<N).all() and (y>=0).all() and (y+1<N).all()
 keys=np.concatenate([y.astype('int64')*N+x,y.astype('int64')*N+x+1,(y+1).astype('int64')*N+x,(y+1).astype('int64')*N+x+1])
 values=np.concatenate([frac*(1-a)*(1-b),frac*a*(1-b),frac*(1-a)*b,frac*a*b]).astype(np.float32)*(0.3/g)**2
 uniq,inv=np.unique(keys,return_inverse=True)
 values=np.bincount(inv,weights=values,minlength=len(uniq)).astype(np.float32)
 excess=np.maximum(values-1,0).sum()*(g/3600)**2
 values=np.minimum(values,1)
 return uniq,values,float(excess)

out={str(g):[] for g in GRIDS};tiles=[]
with tarfile.open(F,'r|gz') as tar:
 for n,m in enumerate(tar):
  if not m.isfile():continue
  f=tar.extractfile(m);h=hdr(f);nx=h['NAXIS1'];ny=h['NAXIS2']
  assert nx==19200 and ny==24910 and h['BITPIX']==8
  corners_x=np.array([0,nx-1,0,nx-1,(nx-1)/2.])
  corners_y=np.array([0,0,ny-1,ny-1,(ny-1)/2.])
  legacy_shift=legacy_column_scaling_error_arcsec(h,corners_x,corners_y)
  nominal=(.03/3600)**2
  center_area=float(tan_pixel_area_deg2((nx-1)/2.,(ny-1)/2.,h))
  cnt=0; xx=[];yy=[];ff=[]
  for i in range(0,ny,1000):
   rows=min(1000,ny-i)
   dat=f.read(rows*nx)
   a=np.frombuffer(dat,dtype=np.uint8).reshape(rows//10,10,nx//10,10)
   if np.any(a>1): raise ValueError('Nonbinary JWST mask values in '+m.name)
   b=a.sum(axis=(1,3),dtype=np.uint16)
   j,k=np.nonzero(b)
   if not len(j):continue
   v=b[j,k];cnt+=int(v.sum())
   xe,ye=sky(k*10+4.5,(j+i//10)*10+4.5,h)
   xx.append(xe.astype(np.float32)); yy.append(ye.astype(np.float32));ff.append((v/100).astype(np.float32))
  x=np.concatenate(xx);y=np.concatenate(yy);frac=np.concatenate(ff)
  native=cnt*(.03/3600)**2
  t={'tile':m.name,'pixels_code1':cnt,'native_area_deg2':native,
     'wcs_pc_off_diagonal':[float(h['PC1_2']),float(h['PC2_1'])],
     'wcs_cd_det_deg2_per_pixel':float(np.linalg.det(cd_from_header(h))),
     'wcs_center_tan_pixel_area_deg2':center_area,
     'wcs_center_area_ratio_to_nominal_30mas':center_area/nominal,
     'legacy_column_scaled_pc_shift_arcsec_max_corners':float(np.max(legacy_shift))}
  for g in GRIDS:
   key,v,excess=depos(x,y,frac,g)
   out[str(g)].append((key,v))
   dep=v.sum(dtype=np.float64)*(g/3600)**2
   t[f'grid_{g}_area_deg2']=dep
   t[f'grid_{g}_area_ratio_to_native']=dep/native
   t[f'grid_{g}_clipped_excess_deg2']=excess
  tiles.append(t)
  print(f'{n+1}/20 {m.name} native={native:.9f} grid0.3/native={t["grid_0.3_area_ratio_to_native"]:.5f} grid0.6/native={t["grid_0.6_area_ratio_to_native"]:.5f} elapsed={time.time()-T:.1f}s',flush=True)
res={'source_archive_sha256':archive_sha256, 'source_archive_size_bytes':Path(F).stat().st_size, 'grid_origin_shift_in_cells':[args.shift_x_cells,args.shift_y_cells], 'source':'user supplied JWST binary starmask archive; no WHT or MAST used','wcs_convention':'FITS CD[i,j]=CDELT[i]*PC[i,j]; earlier mask union/MC files require WCS validation before scientific use','method':'native 30mas binary mask aggregated in 10x10 blocks, bilinear mass-conserving fractional deposition onto Lambert equal-area grid; per sky cell max occupancy is an APPROXIMATE union lower bound; sum clipped at 1 an APPROXIMATE upper bound; grid resolution sensitivity 0.3/0.6 arcsec','tiles':tiles,'sum_tile_native_masked_area_deg2_NOT_UNION':sum(q['native_area_deg2'] for q in tiles),'grid_results':{}}
for g in GRIDS:
 keys=np.concatenate([kv[0] for kv in out[str(g)]]);values=np.concatenate([kv[1] for kv in out[str(g)]])
 uniq,inv=np.unique(keys,return_inverse=True)
 upper=np.bincount(inv,weights=values,minlength=len(uniq));lower=np.zeros(len(uniq),dtype=np.float32)
 np.maximum.at(lower,inv,values)
 lo=lower.sum(dtype=np.float64)*(g/3600)**2
 hi=np.minimum(upper,1).sum()*(g/3600)**2
 res['grid_results'][str(g)]={'union_lower_bound_deg2':lo,'union_upper_bound_deg2':hi,'occupied_cells':len(uniq)}
 if g == 0.6 and args.sample_out:
  if args.samples < 100:
   raise ValueError('At least 100 Monte Carlo samples required')
  weights=np.minimum(upper,1)
  rng=np.random.default_rng(args.seed)
  picked=uniq[rng.choice(len(uniq),size=args.samples,replace=True,p=weights/weights.sum())]
  gs=g/3600.
  x=(ORIGIN+(picked%round(SIZE/gs)+0.5+args.shift_x_cells)*gs)*R
  y=(ORIGIN+(picked//round(SIZE/gs)+0.5+args.shift_y_cells)*gs)*R
  np.savetxt(args.sample_out,np.column_stack((x,y)),delimiter=',',
             header='x_lambert_rad,y_lambert_rad',comments='',fmt='%.13g')
  res['weighted_mc_sample']={'file':str(args.sample_out),'n':args.samples,
                              'sha256':hashlib.sha256(Path(args.sample_out).read_bytes()).hexdigest(),
                              'seed':args.seed,'source_grid_arcsec':g,
                              'sampling_area_deg2':hi,
                              'weight':'clipped sum across tiles (upper raster approximation)'}

 print('UNION',g,'lo',lo,'hi',hi,'elapsed',time.time()-T,flush=True)
Path(args.output).write_text(json.dumps(res,indent=2)+'\n')
print('DONE',time.time()-T,flush=True)
