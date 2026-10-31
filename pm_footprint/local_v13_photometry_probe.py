"""Read-only, low-I/O preflight for Zhen-Kai v13 photometry and Dave125.

Only headers of giant FITS images and catalogue are read. Four ERR maps
are sampled in 3x3 native-pixel neighborhoods at frozen source centers by
os.pread; no giant ndarray allocation, decompression, or image downloads.
Pixel ERR is *not* an aperture uncertainty or 5-sigma magnitude depth.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

FILTERS = ("F115W", "F150W", "F277W", "F444W")
COMMON_HEADER_KEYS = (
    "BUNIT", "WHTTYPE", "ERRTYPE", "PHOTMJSR", "PIXAR_SR",
    "EXPTIME", "EFFEXPTM", "GAIN", "PHOTFLAM", "PHOTPLAM",
    "ZP", "MAGZP", "ABZP", "NCOMBINE", "DRIZCORR",
)
NOTEBOOK_TERMS = (
    "photmjsr", "flux", "aper", "err", "weight", "noise",
    "psf", "encircled", "ratio", "f277w", "flag", "sigma",
)


def hash_file_small(path, max_bytes=2_000_000):
    if path is None or not path.exists() or path.stat().st_size > max_bytes:
        return None
    return hashlib.sha256(path.read_bytes()).hexdigest()


def primary_header(path):
    with path.open("rb") as stream:
        h = fits.Header.fromfile(stream, endcard=True, padding=True)
        offset = stream.tell()
    shape = [int(h.get("NAXIS" + str(j), 0))
             for j in range(int(h.get("NAXIS", 0)), 0, -1)]
    record = {
        "name": path.name, "bytes": path.stat().st_size,
        "bitpix": h.get("BITPIX"), "shape_yx": shape,
        "data_offset_bytes": offset,
        "keywords": {k: h.get(k) for k in COMMON_HEADER_KEYS},
    }
    wcs = None
    try:
        wcs = WCS(h).celestial
        if not wcs.has_celestial:
            wcs = None
    except Exception as exc:
        record["wcs_error"] = str(exc)[:160]
    record["has_celestial_wcs"] = wcs is not None
    if wcs is not None:
        record["ctype"] = list(wcs.wcs.ctype)
        record["crval_deg"] = [float(v) for v in wcs.wcs.crval]
        record["pixel_scale_matrix_deg"] = (
            wcs.pixel_scale_matrix.astype(float).tolist()
        )
    return h, offset, record, wcs


def read_catalog_schema(path):
    # Read only FITS header blocks, never table data. A standard FITS
    # table has NAXIS=0 primary and BINTABLE immediately after.
    with path.open("rb") as stream:
        p = fits.Header.fromfile(stream, endcard=True, padding=True)
        if p.get("NAXIS", 0) != 0:
            raise ValueError("Catalog primary has image data; stop before scanning 6GB")
        t = fits.Header.fromfile(stream, endcard=True, padding=True)
    if t.get("XTENSION", "").strip().upper() not in ("BINTABLE", "TABLE"):
        raise ValueError("Expected catalog table as first FITS extension")
    ncol = int(t.get("TFIELDS", 0))
    return {
        "file": path.name, "file_size_bytes": path.stat().st_size,
        "rows": int(t.get("NAXIS2", 0)), "row_size_bytes": int(t.get("NAXIS1", 0)),
        "extension": t.get("EXTNAME"), "column_count": ncol,
        "columns": [
            {"name": str(t.get(f"TTYPE{i}", "")),
             "format": str(t.get(f"TFORM{i}", "")),
             "unit": t.get(f"TUNIT{i}")}
            for i in range(1, ncol+1)
        ],
        "note": "Schema/header only; no 6GB table records loaded",
    }


def read_native_patch(fd, offset, h, x, y, radius=1):
    nx, ny = int(h["NAXIS1"]), int(h["NAXIS2"])
    if h.get("BITPIX") != -32 or h.get("BSCALE", 1) != 1 or h.get("BZERO", 0) != 0:
        raise ValueError("ERR image must be unscaled FITS big-endian float32")
    ix, iy = int(np.rint(x)), int(np.rint(y))
    if not (radius <= ix < nx-radius and radius <= iy < ny-radius):
        return None, False
    chunks = []
    for row in range(iy-radius, iy+radius+1):
        nvals = 2*radius+1
        byte_offset = offset + 4*(row*nx + ix-radius)
        data = os.pread(fd, nvals*4, int(byte_offset))
        if len(data) != 4*nvals:
            raise IOError("Truncated ERR FITS read")
        chunks.append(np.frombuffer(data, dtype=">f4").astype(float))
    return np.concatenate(chunks), True


def summarize_patch(patch):
    if patch is None:
        return {"center_err": None, "median_positive_err_3x3": None,
                "finite_count": 0, "positive_count": 0}
    f = np.isfinite(patch)
    positive = f & (patch > 0)
    return {
        "center_err": float(patch[len(patch)//2])
        if f[len(patch)//2] else None,
        "median_positive_err_3x3": float(np.median(patch[positive]))
        if np.any(positive) else None,
        "finite_count": int(f.sum()), "positive_count": int(positive.sum())
    }


def notebook_summary(path):
    if path is None or not path.is_file():
        return None
    nb = json.loads(path.read_text(encoding="utf-8"))
    cells = nb.get("cells", [])
    matches = []
    for i, c in enumerate(cells):
        if c.get("cell_type") != "code":
            continue
        s = "".join(c.get("source", [])).lower()
        terms = [term for term in NOTEBOOK_TERMS if term in s]
        if terms:
            matches.append({"cell_index": i, "keywords": terms})
    return {
        "name": path.name, "sha256": hash_file_small(path),
        "total_cells": len(cells),
        "code_cells": sum(c.get("cell_type") == "code" for c in cells),
        "keyword_hits": matches,
        "note": "No notebook executed and no source code or private paths exported."
    }


def run(root, catalog, roster, out_prefix, notebook=None, ee=None):
    root, catalog, roster = map(Path, (root, catalog, roster))
    if not root.is_dir() or not catalog.is_file() or not roster.is_file():
        raise FileNotFoundError("Root directory, uncompressed catalog and roster must exist")
    with roster.open(newline="") as stream:
        sources = list(csv.DictReader(stream))
    if len(sources) != 125 or len({c["candidate_id"] for c in sources}) != 125:
        raise ValueError("Expected the frozen 125 distinct candidate IDs")
    roster_sha = hashlib.sha256(roster.read_bytes()).hexdigest()
    if roster_sha != "d43ce32e33b898f6bf57551affe41f209e0b455a4e4e63b4e99a8eea03236bbd":
        raise ValueError("Frozen candidate roster SHA256 mismatch")
    schema = read_catalog_schema(catalog)
    positions = np.array([[float(s["ra_deg"]),float(s["dec_deg"])]
                          for s in sources], dtype=float)
    files = {}
    for f in FILTERS:
        for typ, name in (
            ("ERR", f"{f}_psf_matched_F444W_ERR_v13.fits"),
            ("SCI", f"{f}_psf_matched_F444W_SCI_v13.fits"
             if f != "F444W" else
             "F444W_psf_matched_F444W_SCI_v13_bkgsub.fits"),
        ):
            p = root / name
            if not p.is_file():
                raise FileNotFoundError(f"Missing image: {p}")
            header, offset, item, wcs = primary_header(p)
            files[f"{f}_{typ}"] = (p,header,offset,item,wcs)
            if len(item["shape_yx"]) != 2:
                raise ValueError(f"Expected 2D image: {p}")
    pw = root / "F444W_psf_matched_F444W_WHT_v13.fits"
    _, _, wht_item, wht_wcs = primary_header(pw)
    if wht_wcs is None:
        raise ValueError("F444W WHT has no celestial WCS")
    sampling_rows = [
        {"candidate_id": s["candidate_id"], "ra_deg":s["ra_deg"],
         "dec_deg":s["dec_deg"]} for s in sources
    ]
    report = {"status":"LOCAL_V13_4BAND_PREFLIGHT_NO_SELECTION_AREA",
              "catalog":schema,
              "roster_sha256":roster_sha,
              "image_headers":{"F444W_WHT":wht_item},
              "err_map_source_center_qc":{},
              "notebook":notebook_summary(Path(notebook)) if notebook else None,
              "encircled_energy_table":{
                  "file":Path(ee).name,
                  "sha256":hash_file_small(Path(ee)),
                  "note":"PSF encircled energy information available; not yet interpreted"
              } if ee and Path(ee).is_file() else None}
    for key, (path,h,offset,item,w) in files.items():
        report["image_headers"][key] = item
        if w is None:
            raise ValueError(f"No celestial WCS in {path.name}")
        x,y = w.all_world2pix(positions[:,0], positions[:,1], 0)
        xref,yref = wht_wcs.all_world2pix(positions[:,0],positions[:,1],0)
        report["image_headers"][key]["max_center_pixel_offset_from_F444W_WHT"] = (
            float(np.max(np.hypot(x-xref,y-yref))))
        if not key.endswith("_ERR"):
            continue
        outside, positive_centers = 0, 0
        fd = os.open(path, os.O_RDONLY)
        try:
            for i, row in enumerate(sampling_rows):
                arr, inside = read_native_patch(fd,offset,h,x[i],y[i])
                if not inside:
                    outside += 1
                qc = summarize_patch(arr)
                if qc["center_err"] is not None and qc["center_err"] > 0:
                    positive_centers += 1
                for field, val in qc.items():
                    row[f"{key}_{field}"] = val
        finally:
            os.close(fd)
        report["err_map_source_center_qc"][key] = {
            "out_of_bounds_3x3":outside,
            "center_err_positive":positive_centers,
            "warning":"Native pixel ERR values are not catalog aperture/auto-flux uncertainties."
        }
        print(f"Checked {key}: {positive_centers}/125 positive ERR centers",flush=True)
    report["limitations"] = [
        "Only image headers and 3x3 ERR pixel samples at 125 preselected catalog positions.",
        "Raw ERR pixel values cannot be converted to a 26.2-AB 5sigma map without units, aperture definition, empirical correlated noise and photometric calibration.",
        "This does not establish original Dave negative F277W flux handling, catalog flags, blending or dropout completeness.",
        "SCI full images were not loaded; F356W second epoch images were not provided.",
        "Source-centered ERR samples are selection-biased and cannot measure all-sky usable area."
    ]
    out_prefix = Path(out_prefix)
    out_prefix.parent.mkdir(parents=True,exist_ok=True)
    out_prefix.with_suffix(".json").write_text(
        json.dumps(report,indent=2,allow_nan=False)+"\n")
    csv_path=out_prefix.parent / (out_prefix.name+"_err_candidates.csv")
    with csv_path.open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(sampling_rows[0]))
        writer.writeheader();writer.writerows(sampling_rows)
    print("Wrote",out_prefix.with_suffix(".json"),"and",csv_path)
    return report


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--root",type=Path,required=True,
                   help="Zhenkai folder containing uncompressed v13 ERR/SCI/WHT maps")
    p.add_argument("--catalog",type=Path,required=True,
                   help="Uncompressed cosmos_F444W_v13_catalog_v1.fits")
    p.add_argument("--roster",type=Path,default=Path("pm_footprint/dave125_positions.csv"))
    p.add_argument("--notebook",type=Path)
    p.add_argument("--ee",type=Path)
    p.add_argument("--out-prefix",type=Path,
                   default=Path("pm_footprint/results/local_v13_photometry_probe"))
    a=p.parse_args()
    run(a.root,a.catalog,a.roster,a.out_prefix,a.notebook,a.ee)


if __name__=="__main__":
    main()
