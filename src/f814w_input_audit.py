"""Read-only F814W input audit; does not calculate counts or completeness.

Read gzip FITS headers without expanding image arrays. Optionally stream an
ERR/WHT pair to check their numerical relationship with bounded memory.
"""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import astropy
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales


def read_header(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rb') as stream:
        return fits.Header.fromfile(stream)


def header_record(path):
    h = read_header(path)
    keys = ['NAXIS1', 'NAXIS2', 'BITPIX', 'BUNIT', 'FILETYPE', 'FILTER',
            'FILTER1', 'FILTER2', 'PHOTFLAM', 'PHOTPLAM', 'PHOTZPT',
            'EXPTIME', 'WHTTYPE', 'DRIZCORR', 'D001SCAL', 'D001PIXF']
    record = dict(path=str(path), size_bytes=path.stat().st_size,
                  header={k: h.get(k) for k in keys})
    if h.get('NAXIS') == 2:
        w = WCS(h).celestial
        nx, ny = h['NAXIS1'], h['NAXIS2']
        xy = np.array([(x, y) for y in [0, (ny-1)/2, ny-1]
                       for x in [0, (nx-1)/2, nx-1]])
        sky = w.pixel_to_world(xy[:, 0], xy[:, 1])
        record.update(pixel_scale_arcsec=(proj_plane_pixel_scales(w)*3600).tolist(),
                      test_pixels_xy=xy.tolist(),
                      test_sky_radec_deg=np.column_stack([sky.ra.deg, sky.dec.deg]).tolist(),
                      wcs_header=w.to_header().tostring(sep='\n', padding=False))
    return record


def snr_comparison(path, nrows):
    output = {}
    with fits.open(path, memmap=True) as hdus:
        table = hdus[1].data
        # Deterministic sampling throughout the row range, rather than a tile prefix.
        indices = np.linspace(0, len(table)-1, min(nrows, len(table)), dtype=int)
        output['catalog_rows'] = len(table)
        output['sample_rows'] = len(indices)
        output['sampling'] = 'evenly spaced row indices including both endpoints'
        output['hdu_comments'] = list(hdus[1].header.get('COMMENT', []))
        for band in ['hst-f814w', 'f115w', 'f150w', 'f277w', 'f444w']:
            snr = np.array(table['snr_'+band][indices])
            rows = []
            for mode in ['auto', 'aper']:
                flux = np.array(table['flux_'+mode+'_'+band][indices])
                err = np.array(table['flux_err_'+mode+'_'+band][indices])
                if flux.ndim == 1:
                    flux, err = flux[:, None], err[:, None]
                for j in range(flux.shape[1]):
                    valid = (np.isfinite(snr) & (snr != 0) & (snr != -99)
                             & np.isfinite(flux[:, j]) & (flux[:, j] != -99)
                             & np.isfinite(err[:, j]) & (err[:, j] > 0))
                    ratio = flux[valid, j]/err[valid, j]
                    s = snr[valid]
                    rows.append(dict(mode=mode, aperture_index=j,
                                     n_valid=int(valid.sum()),
                                     equal_fraction=float(np.mean(np.isclose(ratio, s, rtol=1e-6, atol=1e-8))),
                                     flux_over_error_divided_by_snr_p16_p50_p84=np.percentile(ratio/s, [16, 50, 84]).tolist()))
            output[band] = rows
    return output


def stream_weight_check(err_path, wht_path):
    """Check ERR^2 WHT on a deterministic pixel subsample, no noise assumption."""
    values = []
    nvalid = 0
    with gzip.open(err_path, 'rb') as es, gzip.open(wht_path, 'rb') as ws:
        eh, wh = fits.Header.fromfile(es), fits.Header.fromfile(ws)
        assert (eh['NAXIS1'], eh['NAXIS2']) == (wh['NAXIS1'], wh['NAXIS2'])
        assert eh['BITPIX'] == wh['BITPIX'] == -32
        assert eh.get('BSCALE', 1) == wh.get('BSCALE', 1) == 1
        assert eh.get('BZERO', 0) == wh.get('BZERO', 0) == 0
        nx, ny = eh['NAXIS1'], eh['NAXIS2']
        for y in range(0, ny, 64):
            count = nx*min(64, ny-y)
            eb, wb = es.read(4*count), ws.read(4*count)
            if len(eb) != 4*count or len(wb) != 4*count:
                raise IOError('Truncated FITS pixel stream')
            e, w = np.frombuffer(eb, dtype='>f4'), np.frombuffer(wb, dtype='>f4')
            valid = np.isfinite(e) & (e > 0) & np.isfinite(w) & (w > 0)
            nvalid += int(valid.sum())
            sampled = np.arange(0, count, 97)
            sampled = sampled[valid[sampled]]
            values.append(e[sampled].astype(float)**2 * w[sampled])
    v = np.concatenate(values)
    return dict(err_path=str(err_path), wht_path=str(wht_path),
                joint_positive_finite_pixels=nvalid, sampled_pixels=len(v),
                err_squared_times_wht_p0_p16_p50_p84_p100=np.percentile(v, [0, 16, 50, 84, 100]).tolist(),
                note='Pixel-level relationship only; no aperture correlated-noise or catalog SNR calibration.')


def run(root, catalog, out, sample_rows, check_weights):
    out.mkdir(parents=True, exist_ok=True)
    inventory = []
    # Deliberately confined to E_202/JWST; no downloads or unrelated drives.
    for dirname in ['COSMOS-Web', 'COSMOS_Web']:
        inventory.extend(p for p in (root/dirname).rglob('*')
                         if p.is_file() and 'f814' in p.name.lower())
    records = []
    for path in sorted(inventory):
        if path.name.endswith('.fits.gz') and path.stat().st_size:
            print('Header:', path.name, flush=True)
            records.append(header_record(path))
    by_path = {x['path']: x for x in records}
    parity = []
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    for r in records:
        if not r['path'].endswith('_drz.fits.gz'):
            continue
        for kind in ['err', 'wht']:
            partner = by_path.get(r['path'].replace('_drz.fits.gz', '_'+kind+'.fits.gz'))
            if partner is None:
                continue
            a = SkyCoord(np.array(r['test_sky_radec_deg'])*u.deg)
            b = SkyCoord(np.array(partner['test_sky_radec_deg'])*u.deg)
            parity.append(dict(science_path=r['path'], other_path=partner['path'],
                               max_nine_point_wcs_difference_arcsec=float(a.separation(b).arcsec.max())))
    psfs = []
    for p in sorted(inventory):
        if p.suffix != '.psf':
            continue
        with fits.open(p, memmap=True) as h:
            psfs.append(dict(path=str(p), header={k:h[1].header.get(k) for k in
                        ['FILTER', 'PSF_SAMP', 'PSF_FWHM', 'POLDEG1', 'PSFAXIS1', 'PSFAXIS2', 'PSFAXIS3', 'ACCEPTED']},
                        note='PSFEx spatial basis, not a directly injectable normalized 2D PSF; matching to July 2026 images remains unvalidated.'))
    report = dict(status='INPUT_AUDIT_ONLY_COUNTS_AND_RECOVERY_NOT_EXECUTED',
                  numpy_version=np.__version__, astropy_version=astropy.__version__,
                  root=str(root), catalog_path=str(catalog),
                  inventory=[dict(path=str(p),size_bytes=p.stat().st_size) for p in sorted(inventory)],
                  fits_headers=records, wcs_pair_checks=parity, psfs=psfs,
                  catalog_snr_comparison=snr_comparison(catalog, sample_rows),
                  final_five_band_expected_count=None, common_five_band_area_deg2=None,
                  warning='Header rectangles are not valid survey coverage. No footprint, depth, injection, or count result is asserted.')
    if check_weights:
        tile_dir = root/'COSMOS_Web/HST_F814W_v0.5.2026jul/tiles'
        e = next(tile_dir.glob('*_60mas_tile_A1_*_err.fits.gz'))
        w = next(tile_dir.glob('*_60mas_tile_A1_*_wht.fits.gz'))
        print('Streaming ERR/WHT pixels:', e.name, flush=True)
        report['weight_pixel_check'] = stream_weight_check(e, w)
    (out/'input_audit.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print('Saved', out/'input_audit.json', flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=Path('/mnt/E_202/JWST'))
    p.add_argument('--catalog', type=Path, default=Path('/home/bahareh/Data/Cosmos_Web/COSMOSWeb_mastercatalog_v1.1.fits'))
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--sample-rows', type=int, default=60000)
    p.add_argument('--check-weights', action='store_true')
    args = p.parse_args()
    run(args.root, args.catalog, args.out, args.sample_rows, args.check_weights)
