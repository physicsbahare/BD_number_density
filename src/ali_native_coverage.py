"""Native-pixel, overlap-deduplicated positive-weight footprints, local files only."""
from __future__ import annotations
import argparse,gzip,json,time
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import minimum_filter
from src.depth_reconstruction import read_header,aperture_kernel

BANDS=['F115W','F150W','F277W','F444W','F814W']


def aperture_validity(mask):
    mask=np.asarray(mask,dtype=bool)
    # Coverage is dominated by connected covered/uncovered regions. Only bad
    # pixels adjacent to coverage can affect erosion of covered pixels.
    adjacent=np.zeros(mask.shape,bool)
    adjacent[1:]|=mask[:-1];adjacent[:-1]|=mask[1:]
    adjacent[:,1:]|=mask[:,:-1];adjacent[:,:-1]|=mask[:,1:]
    boundary=np.flatnonzero(adjacent&~mask)
    if len(boundary)>2000000:
        return minimum_filter(mask,size=(5,7),mode='constant',cval=0)&minimum_filter(mask,size=(7,5),mode='constant',cval=0)
    out=mask.copy();ny,nx=mask.shape;yy,xx=np.divmod(boundary,nx)
    for dy in range(-3,4):
        for dx in range(-3,4):
            if abs(dx)==abs(dy)==3:continue
            y,x=yy+dy,xx+dx;ok=(y>=0)&(y<ny)&(x>=0)&(x<nx)
            out[y[ok],x[ok]]=False
    out[:3]=False;out[-3:]=False;out[:,:3]=False;out[:,-3:]=False
    return out


def run(args):
    out=args.out;out.mkdir(parents=True,exist_ok=True)
    nir=args.root/'COSMOS-Web/NIRCam_mosaics'
    hst=args.root/'COSMOS_Web/HST_F814W_v0.5.2026jul/tiles'
    tiles=[f'{a}{i}' for a in 'AB' for i in range(1,11)]
    saved=json.loads((out/'geometry.json').read_text()) if (out/'geometry.json').exists() else None
    headers={} if saved else {t:read_header(next(nir.glob('*f444w*30mas_'+t+'_*_sci.fits.gz'))) for t in tiles}
    gh=fits.Header.fromstring(saved['global_wcs'],sep='\n') if saved else headers['A1'].copy()
    gh['CRPIX1']=gh['CRPIX2']=1.;gw=WCS(gh).celestial
    records=saved['inputs'] if saved else []
    for b in ([] if saved else BANDS):
        for t in tiles:
            sci=next((hst if b=='F814W' else nir).glob('*_30mas_tile_'+t+'_*_drz.fits.gz' if b=='F814W' else '*'+b.lower()+'_*30mas_'+t+'_*_sci.fits.gz'))
            sh=read_header(sci);nx,ny=sh['NAXIS1'],sh['NAXIS2'];sw=WCS(sh).celestial
            xy=np.array([[0,0],[nx-1,0],[(nx-1)/2,(ny-1)/2],[0,ny-1],[nx-1,ny-1]])
            xx,yy=gw.world_to_pixel(sw.pixel_to_world(xy[:,0],xy[:,1]))
            offsets=np.column_stack([xx,yy])-xy
            med=np.median(offsets,axis=0);res=np.max(np.abs(offsets-med))
            if res>.05:raise ValueError('Non-translational WCS mapping needs full resampling')
            offset=np.rint(med).astype(int)
            wht=Path(str(sci).replace('_drz.fits.gz','_wht.fits.gz').replace('_sci.fits.gz','_wht.fits.gz'))
            records.append(dict(band=b,tile=t,weight_path=str(wht),shape=[ny,nx],
                global_offset_xy=offset.tolist(),wcs_translation_residual_native_pixel=float(res),
                fractional_translation_xy=(med-offset).tolist()))
    step=args.step
    xmin=int(np.floor((min(r['global_offset_xy'][0] for r in records)-step)/step)*step-step//2)
    ymin=int(np.floor((min(r['global_offset_xy'][1] for r in records)-step)/step)*step-step//2)
    xmax=max(r['global_offset_xy'][0]+r['shape'][1] for r in records)+step
    ymax=max(r['global_offset_xy'][1]+r['shape'][0] for r in records)+step
    nx=int(np.ceil((xmax-xmin)/step))*step;ny=int(np.ceil((ymax-ymin)/step))*step
    if saved:
        xmin,ymin=saved['global_origin_xy'];ny,nx=saved['native_shape_yx']
    nbyte=(nx+7)//8
    meta=dict(global_origin_xy=[xmin,ymin],native_shape_yx=[ny,nx],step=step,
              global_wcs=gw.to_header().tostring(sep='\n',padding=False),inputs=records,
              definition='Finite positive WHT; native 0.10 arcsec-radius aperture support eroded per tile before union. Header-WCS translations rounded to nearest native pixel; residual and fractional offsets retained.')
    (out/'geometry.json').write_text(json.dumps(meta,indent=2)+'\n')
    struct=aperture_kernel(.1,.03)>0
    for band in BANDS:
        marker=out/(band+'.complete')
        if marker.exists():print('Reusing native coverage',band,flush=True);continue
        bitpath=out/(band+'.bits');progress_path=out/(band+'.tiles.json')
        completed=json.loads(progress_path.read_text()) if progress_path.exists() else []
        bitmap=np.memmap(bitpath,dtype='uint8',mode='r+' if completed else 'w+',shape=(ny,nbyte))
        if not completed:bitmap[:]=0
        for r in records:
            if r['band']!=band:continue
            if r['tile'] in completed:continue
            if band=='F814W' and args.registration_dir is not None:
                rp=args.registration_dir/(r['tile']+'.json')
                while not rp.exists():time.sleep(15)
                delta=np.array(json.loads(rp.read_text())['F814W']['native_static_registration_offset_xy_pixels'])
                header_offset=r.get('header_global_offset_xy',r['global_offset_xy'])
                actual=np.array(header_offset)+np.array(r['fractional_translation_xy'])-delta
                r['header_global_offset_xy']=header_offset
                r['applied_static_registration_xy']=delta.tolist()
                r['global_offset_xy']=np.rint(actual).astype(int).tolist()
            print('Native coverage',band,r['tile'],flush=True)
            sy,sx=r['shape'];valid=np.empty((sy,sx),bool)
            cache=args.cache/r['tile']/(band+'_wht_coverage.npz')
            cache.parent.mkdir(parents=True,exist_ok=True)
            if cache.exists():
                with np.load(cache) as z:packed=z['packed']
                for kind in ['sci','err']:
                    cp=args.cache/r['tile']/(band+'_'+kind+'_coverage.npz')
                    if cp.exists():
                        with np.load(cp) as z:packed&=z['packed']
                valid[:]=np.unpackbits(packed,axis=1)[:,:sx]
                print('Reused packed validity',band,r['tile'],flush=True)
            else:
                with gzip.open(r['weight_path'],'rb') as f:
                    h=fits.Header.fromfile(f)
                    if h['BITPIX']!=-32:raise ValueError('Unexpected WHT BITPIX')
                    for y in range(0,sy,128):
                        n=min(128,sy-y);raw=f.read(n*sx*4)
                        if len(raw)!=n*sx*4:raise IOError('Truncated WHT')
                        a=np.frombuffer(raw,'>f4').reshape(n,sx)
                        valid[y:y+n]=np.isfinite(a)&(a>0)
                np.savez_compressed(cache,packed=np.packbits(valid,axis=1),shape_yx=[sy,sx])
            # Exact decomposition of the 45-pixel fractional-aperture support:
            # union of 5x7 and 7x5 rectangles, avoiding a slow 45-neighbor loop.
            valid=aperture_validity(valid)
            xo,yo=np.array(r['global_offset_xy'])-[xmin,ymin]
            byte0,pad=divmod(int(xo),8)
            for y in range(0,sy,128):
                a=valid[y:y+128]
                packed=np.packbits(np.pad(a,((0,0),(pad,0))),axis=1)
                target=bitmap[yo+y:yo+y+len(a),byte0:byte0+packed.shape[1]]
                target|=packed
            bitmap.flush()
            del valid
            completed.append(r['tile']);progress_path.write_text(json.dumps(completed)+'\n')
        marker.write_text('complete\n');del bitmap
    (out/'geometry.json').write_text(json.dumps(meta,indent=2)+'\n')
    maps=[np.memmap(out/(b+'.bits'),dtype='uint8',mode='r',shape=(ny,nbyte)) for b in BANDS]
    four=np.zeros((ny//step,nx//step));five=np.zeros_like(four)
    # TAN solid-angle Jacobian; rotation is orthogonal, so radius is invariant.
    pixrad=np.deg2rad(.03/3600)
    xangle=(np.arange(nx)+xmin)*pixrad
    for y in range(0,ny,step):
        joint=maps[0][y:y+step].copy()
        for m in maps[1:4]:joint&=m[y:y+step]
        allfive=joint&maps[4][y:y+step]
        j4=np.unpackbits(joint,axis=1)[:,:nx]
        j5=np.unpackbits(allfive,axis=1)[:,:nx]
        yangle=(np.arange(y,y+step)+ymin)*pixrad
        jac=(1+yangle[:,None]**2+xangle[None,:]**2)**-1.5
        four[y//step]=(j4*jac).reshape(step,nx//step,step).sum(axis=(0,2))*pixrad**2
        five[y//step]=(j5*jac).reshape(step,nx//step,step).sum(axis=(0,2))*pixrad**2
    np.savez_compressed(out/'area_maps.npz',four_area_sr=four,five_area_sr=five,
                        global_origin_xy=[xmin,ymin],step=step)
    summary=dict(four_area_deg2=float(four.sum()*(180/np.pi)**2),five_area_deg2=float(five.sum()*(180/np.pi)**2),
                 geometry=meta,common_area_never_exceeds_four=bool(np.all(five<=four+1e-20)))
    (out/'coverage_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('Native areas:',summary['four_area_deg2'],summary['five_area_deg2'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path('/mnt/E_202/JWST'))
    p.add_argument('--out',type=Path,default=Path('/home/bahareh/results/ali_native_coverage'))
    p.add_argument('--step',type=int,default=100)
    p.add_argument('--cache',type=Path,default=Path('/home/bahareh/results/ali_fiveband_reconstruction'))
    p.add_argument('--registration-dir',type=Path)
    run(p.parse_args())
