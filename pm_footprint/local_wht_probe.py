"""Local, read-only FITS WHT inventory; never loads a full 50 GB image.

Run from repository root:
 python -m pm_footprint.local_wht_probe --v13 "/path/F444W_psf_matched_F444W_WHT_v13.fits" --cosmos-dir "/path/NIRCam_mosaics" --out pm_footprint/results/local_wht_probe.json

For compressed COSMOS-Web tiles this reads headers only; do not mistake them
for calibrated depth maps. No image bytes or coordinates leave the machine.
"""
import argparse
import json
import hashlib
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

FILTERS = ("f115w", "f150w", "f277w", "f444w")


def probe(path, sample=True):
    item = {"file": path.name, "bytes": path.stat().st_size,
            "compressed": path.name.endswith(".gz"), "hdus": []}
    with fits.open(path, memmap=not item["compressed"], lazy_load_hdus=True) as hdul:
        for i, hdu in enumerate(hdul):
            hdr = hdu.header
            shape = [int(hdr.get("NAXIS"+str(j), 0))
                     for j in range(int(hdr.get("NAXIS", 0)), 0, -1)]
            d = {"index": i, "extname": hdr.get("EXTNAME"),
                 "shape": shape, "bitpix": hdr.get("BITPIX"),
                 "bunit": hdr.get("BUNIT"), "wht_type": hdr.get("WHTTYPE"),
                 "photmjsr": hdr.get("PHOTMJSR"),
                 "has_celestial_wcs": False}
            try:
                w = WCS(hdr)
                d["has_celestial_wcs"] = bool(w.has_celestial)
                if w.has_celestial:
                    d["pixel_scale_matrix_deg"] = w.celestial.pixel_scale_matrix.tolist()
                    d["crval"] = [float(x) for x in w.celestial.wcs.crval]
            except Exception as exc:
                d["wcs_error"] = str(exc)[:150]
            if sample and not item["compressed"] and len(shape) == 2 and min(shape) > 10:
                # Memmap + scalar reads: no whole-image array materialization.
                data = hdu.data
                ys = np.linspace(0, shape[0]-1, 11, dtype=int)
                xs = np.linspace(0, shape[1]-1, 11, dtype=int)
                vals = np.array([float(data[y, x]) for y in ys for x in xs])
                good = np.isfinite(vals)
                d["sparse_121_samples"] = {
                    "finite": int(good.sum()),
                    "positive": int(np.count_nonzero(vals[good] > 0)),
                    "min_finite": float(vals[good].min()) if good.any() else None,
                    "median_positive": float(np.median(vals[good & (vals > 0)]))
                    if np.any(good & (vals > 0)) else None,
                }
            item["hdus"].append(d)
    return item


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--v13", type=Path, required=True)
    p.add_argument("--cosmos-dir", type=Path, required=True)
    p.add_argument("--out", type=Path, default=Path("pm_footprint/results/local_wht_probe.json"))
    args = p.parse_args()
    if not args.v13.is_file():
        p.error("v13 WHT file does not exist")
    tiles = sorted(x for x in args.cosmos_dir.glob("*wht.fits*")
                   if any(f in x.name.lower() for f in FILTERS))
    if not tiles:
        p.error("No COSMOS-Web WHT tiles found; check --cosmos-dir")
    result = {"purpose": "metadata and sparse QC only, NOT effective area",
              "v13": probe(args.v13), "cosmos_web": [probe(t, sample=False) for t in tiles],
              "tile_count": len(tiles),
              "warning": "WHT>0 is geometric coverage, not 26.2 AB completeness. "
                         "Do not compare WHT values across reductions without unit/noise calibration."}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print("Wrote", args.out, "with", len(tiles), "COSMOS-Web tile headers")


if __name__ == "__main__":
    main()
