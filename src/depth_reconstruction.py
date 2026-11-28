"""Bounded-memory extraction from local gzip mosaics for Ali depth rebuilding."""
from __future__ import annotations

import argparse
import gzip
import json
import time
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from photutils.aperture import CircularAperture
from scipy.spatial import cKDTree


def aperture_kernel(radius_arcsec, scale_arcsec):
    r = radius_arcsec/scale_arcsec
    half = int(np.ceil(r))+1
    return CircularAperture((half, half), r).to_mask(method='exact').to_image((2*half+1,2*half+1))


def streamed_samples(path, center_xy, stamp_xy, half, chunk_rows=128,coverage_path=None,positive_coverage=False):
    """Extract central pixels and postage stamps in a single forward pass.

    Coordinates are zero based integer pixels; stamps crossing edges become
    NaN. FITS BSCALE/BZERO are honored. No full mosaic array is allocated.
    """
    opener = gzip.open if str(path).endswith('.gz') else open
    center_xy = np.asarray(center_xy, dtype=int).reshape(-1,2)
    stamp_xy = np.asarray(stamp_xy, dtype=int).reshape(-1,2)
    centers = np.full(len(center_xy), np.nan, dtype='float32')
    stamps = np.full((len(stamp_xy),2*half+1,2*half+1), np.nan, dtype='float32')
    with opener(path,'rb') as f:
        h = fits.Header.fromfile(f)
        if h['NAXIS'] != 2:
            raise ValueError('Expected primary 2D FITS image')
        nx,ny = h['NAXIS1'],h['NAXIS2']
        packed=np.empty((ny,(nx+7)//8),dtype='uint8') if coverage_path is not None else None
        dt = np.dtype({-32:'>f4',-64:'>f8',32:'>i4',16:'>i2',8:'u1'}[h['BITPIX']])
        for y0 in range(0,ny,chunk_rows):
            y1 = min(ny,y0+chunk_rows)
            raw = f.read((y1-y0)*nx*dt.itemsize)
            if len(raw) != (y1-y0)*nx*dt.itemsize:
                raise IOError('Truncated FITS image: '+str(path))
            a = np.frombuffer(raw,dt).reshape(y1-y0,nx)
            bs,bz = h.get('BSCALE',1),h.get('BZERO',0)
            if packed is not None:
                flag=np.isfinite(a)
                if positive_coverage:flag&=(a.astype(float)*bs+bz)>0
                packed[y0:y1]=np.packbits(flag,axis=1)
            use = np.flatnonzero((center_xy[:,1]>=y0)&(center_xy[:,1]<y1)&(center_xy[:,0]>=0)&(center_xy[:,0]<nx))
            centers[use] = a[center_xy[use,1]-y0,center_xy[use,0]].astype('float64')*bs+bz
            use = np.flatnonzero((stamp_xy[:,1]+half>=y0)&(stamp_xy[:,1]-half<y1)&(stamp_xy[:,0]>=half)&(stamp_xy[:,0]<nx-half)&(stamp_xy[:,1]>=half)&(stamp_xy[:,1]<ny-half))
            for i in use:
                x,y = stamp_xy[i]
                low,high = max(y-half,y0),min(y+half+1,y1)
                stamps[i,low-y+half:high-y+half,:] = a[low-y0:high-y0,x-half:x+half+1].astype('float64')*bs+bz
    if packed is not None:np.savez_compressed(coverage_path,packed=packed,shape_yx=[ny,nx])
    return h,centers,stamps


def read_header(path):
    with gzip.open(path,'rb') as f:
        return fits.Header.fromfile(f)


def geometry(header,step):
    # COSMOS-Web tiles share CRVAL/PC/CDELT; CRPIX encodes integer offsets.
    cr = np.array([header['CRPIX1'],header['CRPIX2']])
    if not np.allclose(cr,np.rint(cr),rtol=0,atol=1e-6):
        raise ValueError('Noninteger tile offset needs general sky resampling')
    nx,ny = header['NAXIS1'],header['NAXIS2']
    gx = np.arange(int(np.ceil((1-cr[0])/step)),int(np.floor((nx-cr[0])/step))+1)*step
    gy = np.arange(int(np.ceil((1-cr[1])/step)),int(np.floor((ny-cr[1])/step))+1)*step
    xx,yy = np.meshgrid(gx,gy)
    global_xy = np.column_stack([xx.ravel(),yy.ravel()]).astype('int32')
    return global_xy, (global_xy+cr-1).astype(int)


def run_tile(tile,args):
    out = args.out/tile
    done = out/'COMPLETE.json'
    if done.exists():
        print('Reusing',tile,flush=True)
        return
    out.mkdir(parents=True,exist_ok=True)
    nir = args.root/'COSMOS-Web/NIRCam_mosaics'
    hst = args.root/'COSMOS_Web/HST_F814W_v0.5.2026jul/tiles'
    ref = next(nir.glob('mosaic_nircam_f444w_*30mas_'+tile+'_*_sci.fits.gz'))
    header = read_header(ref)
    global_xy,grid_xy = geometry(header,args.step)
    nx,ny = header['NAXIS1'],header['NAXIS2']
    with fits.open(args.catalog,memmap=True) as h:
        d = h[1].data
        use = np.flatnonzero(np.char.strip(d['tile'].astype(str))==tile)
        ra,dec = np.array(d['ra'][use]),np.array(d['dec'][use])
        cx,cy = WCS(header).celestial.world_to_pixel_values(ra,dec)
        cat_xy = np.column_stack([cx,cy])
        finite = np.isfinite(cat_xy).all(axis=1)
        tree = cKDTree(cat_xy[finite])
        rng = np.random.default_rng(20261126+int(tile[1:])+(100 if tile[0]=='B' else 0))
        candidates = rng.integers([30,30],[nx-30,ny-30],size=(args.blanks*3,2))
        # Supplement segmentation masks with a 0.75 arcsec catalog-source veto.
        distance,_ = tree.query(candidates)
        blank_xy = candidates[distance>25][:args.blanks]
        choose = np.linspace(0,len(use)-1,min(1000,len(use)),dtype=int)
        probes = np.rint(cat_xy[choose]).astype(int)
        cat_ids = np.array(d['id'][use[choose]])
        cat_snrs = np.column_stack([d['snr_'+b][use[choose]] for b in ['hst-f814w','f115w','f150w','f277w','f444w']])
        star = np.asarray(d['flag_star'][use])>0
        bright = (np.asarray(d['snr_hst-f814w'][use])>30)&(np.asarray(d['snr_hst-f814w'][use])<500)
        distance2,_ = tree.query(cat_xy[finite],k=2)
        isolated = np.zeros(len(use),bool);isolated[finite]=distance2[:,1]>80
        star_pick = np.flatnonzero(star&bright&isolated)[:150]
        star_xy = np.rint(cat_xy[star_pick]).astype(int)
    stamp_xy = np.vstack([blank_xy,probes])
    np.savez_compressed(out/'positions.npz',global_xy=global_xy,grid_xy=grid_xy,blank_xy=blank_xy,
                        catalog_xy=probes,catalog_id=cat_ids,catalog_snr=cat_snrs,
                        star_xy=star_xy,star_radec=np.column_stack([ra[star_pick],dec[star_pick]]))
    report = dict(tile=tile,step_native_pixels=args.step,blank_requested=args.blanks,
                  n_blank=len(blank_xy),n_catalog=len(probes),n_stars=len(star_xy),
                  science_wcs=header.tostring(sep='\n',padding=False),bands={})
    # Segmentation-map stamps are checked rather than assuming blank sky.
    seg = next((args.root/'COSMOS-Web/auxiliary/segmentation_maps').glob('*_'+tile+'_segmap_*.fits.gz'))
    sh,_,ss = streamed_samples(seg,[],blank_xy,7)
    if (sh['NAXIS1'],sh['NAXIS2']) != (nx,ny):
        raise ValueError('Segmentation grid dimensions differ')
    np.save(out/'blank_seg_clear.npy',np.isfinite(ss).all(axis=(1,2))&np.all(ss==0,axis=(1,2)))
    for band in ['F115W','F150W','F277W','F444W','F814W']:
        cache = out/(band+'.npz')
        if cache.exists():
            print('Reusing',tile,band,flush=True);continue
        paths = {}
        for kind in ['sci','err','wht']:
            paths[kind] = next((hst if band=='F814W' else nir).glob(
                '*_30mas_tile_'+tile+'_*_'+('drz' if kind=='sci' else kind)+'.fits.gz'
                if band=='F814W' else 'mosaic_nircam_'+band.lower()+'_*30mas_'+tile+'_*_'+kind+'.fits.gz'))
        native_header=read_header(paths['sci'])
        native_wcs=WCS(native_header).celestial
        reference_wcs=WCS(header).celestial
        # ACS headers declare FK5 while NIRCam declares ICRS. Honor the
        # frame transformation; identical PC/CRPIX is not identical sky WCS.
        if band=='F814W':
            gsky=reference_wcs.pixel_to_world(grid_xy[:,0],grid_xy[:,1])
            ssky=reference_wcs.pixel_to_world(stamp_xy[:,0],stamp_xy[:,1])
            native_grid=np.rint(np.column_stack(native_wcs.world_to_pixel(gsky))).astype(int)
            native_stamps=np.rint(np.column_stack(native_wcs.world_to_pixel(ssky))).astype(int)
            if len(star_xy):
                stsky=reference_wcs.pixel_to_world(star_xy[:,0],star_xy[:,1])
                native_stars=np.rint(np.column_stack(native_wcs.world_to_pixel(stsky))).astype(int)
            else:native_stars=star_xy
        else:
            native_grid,native_stamps,native_stars=grid_xy,stamp_xy,star_xy
        vals = {}
        t0 = time.time()
        for kind,path in paths.items():
            print(tile,band,kind,path.name,flush=True)
            bh,center,stamps = streamed_samples(path,native_grid,native_stamps,7,
                coverage_path=out/(band+'_'+kind+'_coverage.npz'),positive_coverage=kind!='sci')
            if (bh['NAXIS1'],bh['NAXIS2']) != (nx,ny):
                raise ValueError('Image dimension mismatch')
            if kind=='sci':
                wref,wband = WCS(header).celestial,WCS(bh).celestial
                p = np.array([[0,0],[nx/2,ny/2],[nx-1,ny-1]])
                diff = wref.pixel_to_world(p[:,0],p[:,1]).separation(wband.pixel_to_world(p[:,0],p[:,1])).arcsec
                if diff.max()>1e-4 and band!='F814W':
                    raise ValueError('Science grid differs; cannot align by pixel')
                vals['same_pixel_wcs_separation_arcsec']=diff
                vals['native_grid_xy']=native_grid
                vals['native_stamp_xy']=native_stamps
                conversion = 1. if band=='F814W' else float(bh['PIXAR_SR'])*1e15
                if band=='F814W' and bh.get('BUNIT')!='nJy':
                    raise ValueError('HST science must be nJy per pixel')
                if band!='F814W' and bh.get('BUNIT')!='MJy/sr':
                    raise ValueError('NIRCam science unit mismatch')
                vals['science_to_njy_per_pixel'] = conversion
            if kind in ['sci','err']:
                center = center*vals['science_to_njy_per_pixel']
                stamps = stamps*vals['science_to_njy_per_pixel']
            vals[kind+'_grid'] = center
            vals[kind+'_stamps'] = stamps
        # PSF validation stars need larger native F814W scenes only.
        if band=='F814W' and len(star_xy):
            _,_,stars = streamed_samples(paths['sci'],[],native_stars,35)
            vals['star_stamps'] = stars
            vals['native_star_xy']=native_stars
        np.savez_compressed(cache,**vals)
        report['bands'][band]=dict(paths={k:str(v) for k,v in paths.items()},seconds=time.time()-t0)
        (out/'progress.json').write_text(json.dumps(report,indent=2)+'\n')
        print('Saved',cache,'elapsed_s',round(time.time()-t0),flush=True)
    done.write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path('/mnt/E_202/JWST'))
    p.add_argument('--catalog',type=Path,default=Path('/home/bahareh/Data/Cosmos_Web/COSMOSWeb_mastercatalog_v1.1.fits'))
    p.add_argument('--out',type=Path,default=Path('/home/bahareh/results/ali_fiveband_reconstruction'))
    p.add_argument('--tiles',nargs='+',default=[f'{a}{i}' for a in 'AB' for i in range(1,11)])
    p.add_argument('--step',type=int,default=100)
    p.add_argument('--blanks',type=int,default=10000)
    args=p.parse_args()
    for tile in args.tiles:run_tile(tile,args)
