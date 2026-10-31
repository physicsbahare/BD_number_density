"""Estimate JWST star-mask sky union from user-supplied COSMOS2025 30mas FITS.
No WHT download or MAST queries. Requires only numpy.
Usage: python -m pm_footprint.jwst_mask_union_repro --archive FILE.tar.gz
0.3"/0.6" results are raster approximations, not a calibrated dropout effective area.
"""
import argparse, tarfile, numpy as np, time, json
from pathlib import Path
parser=argparse.ArgumentParser()
parser.add_argument('--archive',required=True)
parser.add_argument('--output',default='jwst_mask_union_conservative.json')
args=parser.parse_args()
F=args.archive;R=np.pi/180
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
 u=(x+1-h['CRPIX1'])*h['CDELT1']*R;v=(y+1-h['CRPIX2'])*h['CDELT2']*R
 xi=h['PC1_1']*u+h['PC1_2']*v; eta=h['PC2_1']*u+h['PC2_2']*v
 a=h['CRVAL1']*R;d=h['CRVAL2']*R;sd=np.sin(d);cd=np.cos(d)
 ra=a+np.arctan2(xi,cd-eta*sd)
 dec=np.arctan2(sd+eta*cd,np.sqrt((cd-eta*sd)**2+xi**2))
 dra=ra-r0;k=np.sqrt(2/(1+s0*np.sin(dec)+c0*np.cos(dec)*np.cos(dra)))
 return k*np.cos(dec)*np.sin(dra)/R, k*(c0*np.sin(dec)-s0*np.cos(dec)*np.cos(dra))/R

def depos(xe,ye,frac,g):
 gs=g/3600.;N=round(SIZE/gs)
 fx=(xe-ORIGIN)/gs-0.5;fy=(ye-ORIGIN)/gs-0.5
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
  cnt=0; xx=[];yy=[];ff=[]
  for i in range(0,ny,1000):
   rows=min(1000,ny-i)
   dat=f.read(rows*nx)
   a=np.frombuffer(dat,dtype=np.uint8).reshape(rows//10,10,nx//10,10)
   b=a.sum(axis=(1,3),dtype=np.uint16)
   j,k=np.nonzero(b)
   if not len(j):continue
   v=b[j,k];cnt+=int(v.sum())
   xe,ye=sky(k*10+4.5,(j+i//10)*10+4.5,h)
   xx.append(xe.astype(np.float32)); yy.append(ye.astype(np.float32));ff.append((v/100).astype(np.float32))
  x=np.concatenate(xx);y=np.concatenate(yy);frac=np.concatenate(ff)
  native=cnt*(.03/3600)**2
  t={'tile':m.name,'pixels_code1':cnt,'native_area_deg2':native}
  for g in GRIDS:
   key,v,excess=depos(x,y,frac,g)
   out[str(g)].append((key,v))
   dep=v.sum(dtype=np.float64)*(g/3600)**2
   t[f'grid_{g}_area_deg2']=dep
   t[f'grid_{g}_area_ratio_to_native']=dep/native
   t[f'grid_{g}_clipped_excess_deg2']=excess
  tiles.append(t)
  print(f'{n+1}/20 {m.name} native={native:.9f} grid0.3/native={t["grid_0.3_area_ratio_to_native"]:.5f} grid0.6/native={t["grid_0.6_area_ratio_to_native"]:.5f} elapsed={time.time()-T:.1f}s',flush=True)
res={'source':'user supplied JWST binary starmask archive; no WHT or MAST used','method':'native 30mas binary mask aggregated in 10x10 blocks, bilinear mass-conserving fractional deposition onto Lambert equal-area grid; per sky cell max occupancy is an APPROXIMATE union lower bound; sum clipped at 1 an APPROXIMATE upper bound; grid resolution sensitivity 0.3/0.6 arcsec','tiles':tiles,'sum_tile_native_masked_area_deg2_NOT_UNION':sum(q['native_area_deg2'] for q in tiles),'grid_results':{}}
for g in GRIDS:
 keys=np.concatenate([kv[0] for kv in out[str(g)]]);values=np.concatenate([kv[1] for kv in out[str(g)]])
 uniq,inv=np.unique(keys,return_inverse=True)
 upper=np.bincount(inv,weights=values,minlength=len(uniq));lower=np.zeros(len(uniq),dtype=np.float32)
 np.maximum.at(lower,inv,values)
 lo=lower.sum(dtype=np.float64)*(g/3600)**2
 hi=np.minimum(upper,1).sum()*(g/3600)**2
 res['grid_results'][str(g)]={'union_lower_bound_deg2':lo,'union_upper_bound_deg2':hi,'occupied_cells':len(uniq)}
 print('UNION',g,'lo',lo,'hi',hi,'elapsed',time.time()-T,flush=True)
Path(args.output).write_text(json.dumps(res,indent=2)+'\n')
print('DONE',time.time()-T,flush=True)
