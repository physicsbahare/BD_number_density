"""Audit detection-image zero pixels against already computed WHT coverage.

This reads detection products only. It neither reconstructs survey depths nor
injects sources. The packed-mask union and cached solid angles are reused.
"""
import argparse,json,tarfile,time
from pathlib import Path
import numpy as np
from astropy.io import fits
from src.depth_reconstruction import aperture_kernel


def packed_count_cells(bits,nx,step):
    a=np.unpackbits(bits,axis=1)[:,:nx]
    return a.reshape(len(a),nx//step,step).sum(axis=(0,2))


def run(a):
    a.out.mkdir(parents=True,exist_ok=True);meta=json.loads((a.coverage/'geometry.json').read_text())
    ny,nx=meta['native_shape_yx'];nb=(nx+7)//8;origin=np.array(meta['global_origin_xy']);step=meta['step']
    bp=a.out/'detection_positive.bits';marker=a.out/'archive_complete.json'
    if marker.exists():records=json.loads(marker.read_text())['tiles']
    elif a.reuse_completed_grid_cache:
        expected=[f'{b}{i}' for b in 'AB' for i in range(1,11)]
        if bp.stat().st_size!=ny*nb:raise ValueError('Packed bitmap size mismatch')
        bitmap=np.memmap(bp,mode='r',dtype='uint8',shape=(ny,nb));records=[]
        pilot=json.loads(Path('/home/bahareh/results/ali_detection_zero_pilot.json').read_text())
        for tile in expected:
            with np.load(a.out/(tile+'_detection_grid.npz')) as z:
                xy=z['global_xy'].astype(int)-origin;chi=z['chi']
                known=((bitmap[xy[:,1],xy[:,0]//8]>>(7-xy[:,0]%8))&1)>0
                if not np.all(known[chi>0]):raise ValueError('Completed positive-mask/grid check failed')
            records.append(dict(tile=tile,member='detection_images/detection_chi2pos_SWLW_'+tile+'.fits',
                NSCI_MU=pilot['native_headers']['NSCI_MU'] if tile=='A1' else None,
                NSCI_SIG=pilot['native_headers']['NSCI_SIG'] if tile=='A1' else None,
                recovery='All FITS members completed; bitmap and twenty grids verified after an unhandled trailing non-FITS archive member. Unpersisted per-tile header diagnostics are not fabricated.'))
        del bitmap;marker.write_text(json.dumps(dict(tiles=records),indent=2)+'\n')
    else:
        bitmap=np.memmap(bp,mode='w+',dtype='uint8',shape=(ny,nb));bitmap[:]=0;records=[]
        with tarfile.open(a.archive,'r|gz') as tar:
            for member in tar:
                if not member.isfile() or not member.name.endswith('.fits'):continue
                tile=member.name.split('_')[-1].split('.')[0];print('Detection image',tile,flush=True)
                f=tar.extractfile(member);h=fits.Header.fromfile(f);sx,sy=h['NAXIS1'],h['NAXIS2'];t0=time.time()
                rec=next(x for x in meta['inputs'] if x['band']=='F444W' and x['tile']==tile);xo,yo=np.array(rec['global_offset_xy'])-origin
                byte0,pad=divmod(int(xo),8)
                with np.load(a.calibration/(tile+'.npz')) as d:grid=d['grid_xy'].astype(int);gxy=d['global_xy']
                chi=np.full(len(grid),np.nan,dtype='float32');counts=np.zeros(3,dtype='int64')
                for y in range(0,sy,128):
                    n=min(128,sy-y);raw=f.read(n*sx*4)
                    if len(raw)!=n*sx*4:raise IOError('Truncated detection image')
                    image=np.frombuffer(raw,'>f4').reshape(n,sx)
                    counts+=np.array([np.count_nonzero(image==0),np.count_nonzero(image>0),np.count_nonzero(~np.isfinite(image))])
                    p=np.packbits(np.pad(np.isfinite(image)&(image>0),((0,0),(pad,0))),axis=1)
                    bitmap[yo+y:yo+y+n,byte0:byte0+p.shape[1]]|=p
                    use=np.flatnonzero((grid[:,1]>=y)&(grid[:,1]<y+n));chi[use]=image[grid[use,1]-y,grid[use,0]]
                np.savez_compressed(a.out/(tile+'_detection_grid.npz'),global_xy=gxy,chi=chi)
                records.append(dict(tile=tile,member=member.name,native_zero_positive_nonfinite=counts.tolist(),NSCI_MU=h.get('NSCI_MU'),NSCI_SIG=h.get('NSCI_SIG'),BUNIT_inherited=h.get('BUNIT'),seconds=time.time()-t0))
                bitmap.flush()
                (a.out/'archive_progress.json').write_text(json.dumps(dict(tiles=records),indent=2)+'\n')
        del bitmap;marker.write_text(json.dumps(dict(tiles=records),indent=2)+'\n')
    pos=np.memmap(bp,mode='r',dtype='uint8',shape=(ny,nb))
    original=[np.memmap(a.original/(b+'.bits'),mode='r',dtype='uint8',shape=(ny,nb)) for b in ['F115W','F150W','F277W','F444W']]
    native=[np.memmap(a.coverage/(b+'.bits'),mode='r',dtype='uint8',shape=(ny,nb)) for b in ['F115W','F150W','F277W','F444W','F814W']]
    with np.load(a.original/'original_four_area_map.npz') as d:origarea=d['original_four_area_sr']
    with np.load(a.coverage/'area_maps.npz') as d:fourarea=d['four_area_sr'];fivearea=d['five_area_sr']
    wrong=[np.zeros_like(origarea) for _ in range(3)]
    for y in range(0,ny,step):
        masks=[]
        for items in [original,native[:4],native]:
            v=items[0][y:y+step].copy()
            for x in items[1:]:v&=x[y:y+step]
            masks.append(v)
        for j,(mask,areas) in enumerate(zip(masks,[origarea,fourarea,fivearea])):
            n=packed_count_cells(mask,nx,step);npositive=packed_count_cells(mask&pos[y:y+step],nx,step)
            wrong[j][y//step]=areas[y//step]*np.divide(npositive,n,out=np.zeros_like(n,dtype=float),where=n>0)
    conv=(180/np.pi)**2
    summary=dict(original_positive_wht_area_deg2=float(origarea.sum()*conv),native_four_area_deg2=float(fourarea.sum()*conv),native_five_area_deg2=float(fivearea.sum()*conv),
        erroneous_chi_positive_original_four_area_deg2=float(wrong[0].sum()*conv),erroneous_chi_positive_native_four_area_deg2=float(wrong[1].sum()*conv),erroneous_chi_positive_five_area_deg2=float(wrong[2].sum()*conv),
        historical_area_deg2=.496306,tiles=records,
        method='Coverage from cached native WHT masks. Detection positivity is an intentionally incorrect alternative. Per-cell exact cached solid angle is scaled by native positive-pixel fraction; within-cell TAN Jacobian variation is negligible relative to mask phase uncertainty.')
    k=aperture_kernel(.1,.03);reports=json.loads((a.calibration/'calibration.json').read_text());diag=[]
    for rec in records:
        if rec.get('NSCI_SIG') is None:continue
        tile=rec['tile'];r=reports[tile]['F444W'];factor=1e15*2.11539874851881e-14
        diag.append(dict(tile=tile,catalog_fitted_scale=r['catalog_noise_equalized_scale'],white_pixel_scale_area=factor*rec['NSCI_SIG']*np.sqrt(k.sum()),white_pixel_scale_fractional_variance=factor*rec['NSCI_SIG']*np.sqrt((k*k).sum()),
            fitted_offset=r['catalog_noise_equalized_offset'],pixel_background_offset=factor*rec['NSCI_MU']*k.sum(),empirical_to_catalog_sigma=1/r['catalog_sigma_over_empirical_sigma_median']))
    summary['native_snr_pixel_noise_diagnostics']=diag
    (a.out/'detection_footprint_validation.json').write_text(json.dumps(summary,indent=2)+'\n')
    np.savez_compressed(a.out/'erroneous_chi_positive_area_maps.npz',original_four_area_sr=wrong[0],four_area_sr=wrong[1],five_area_sr=wrong[2],global_origin_xy=origin,step=step)
    print({k:v for k,v in summary.items() if not isinstance(v,list)},flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--archive',type=Path,default=Path('/mnt/E_202/JWST/COSMOS-Web/auxiliary/detection_images.tar.gz'))
    p.add_argument('--coverage',type=Path,default=Path('/home/bahareh/results/ali_native_coverage'))
    p.add_argument('--original',type=Path,default=Path('/home/bahareh/results/ali_original_coverage'))
    p.add_argument('--calibration',type=Path,default=Path('/home/bahareh/results/ali_calibration_final'))
    p.add_argument('--out',type=Path,default=Path('/home/bahareh/results/ali_detection_validation'))
    p.add_argument('--reuse-completed-grid-cache',action='store_true')
    run(p.parse_args())
