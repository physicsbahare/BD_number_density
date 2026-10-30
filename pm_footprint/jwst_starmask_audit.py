"""Audit user-supplied COSMOS2025 JWST starmask FITS tar.gz against Dave125.
Does not fetch WHT, MAST or any other external file.
Usage: python -m pm_footprint.jwst_starmask_audit --archive /path/to/cosmos_web_starmask_jwst.tar.gz
"""
import argparse
import csv
import hashlib
import json
import math
import tarfile
from pathlib import Path

import numpy as np


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_header(stream):
    header = {}
    while True:
        block = stream.read(2880)
        if len(block) != 2880:
            raise ValueError("Truncated FITS header")
        finished = False
        for pos in range(0, 2880, 80):
            card = block[pos:pos + 80].decode("ascii", "replace")
            if card[:8].strip() == "END":
                finished = True
                break
            if card[8:10] == "= ":
                value = card[10:].split("/", 1)[0].strip().strip("'")
                try:
                    value = float(value) if any(q in value for q in ".Ee") else int(value)
                except ValueError:
                    pass
                header[card[:8].strip()] = value
        if finished:
            return header


def position_pixels(h, ras_deg, decs_deg):
    ra = np.deg2rad(ras_deg)
    dec = np.deg2rad(decs_deg)
    ra0, dec0 = np.deg2rad([h["CRVAL1"], h["CRVAL2"]])
    dra = (ra - ra0 + np.pi) % (2 * np.pi) - np.pi
    cosc = np.sin(dec0) * np.sin(dec) + np.cos(dec0) * np.cos(dec) * np.cos(dra)
    xi = np.cos(dec) * np.sin(dra) / cosc
    eta = (np.cos(dec0) * np.sin(dec) - np.sin(dec0) * np.cos(dec) * np.cos(dra)) / cosc
    cd = np.array([
        [h["PC1_1"] * h["CDELT1"], h["PC1_2"] * h["CDELT1"]],
        [h["PC2_1"] * h["CDELT2"], h["PC2_2"] * h["CDELT2"]],
    ])
    xy = np.linalg.solve(cd, np.rad2deg([xi, eta]))
    return (
        np.rint(xy[0] + h["CRPIX1"] - 1).astype(int),
        np.rint(xy[1] + h["CRPIX2"] - 1).astype(int),
    )


def audit(archive, positions):
    with open(positions, newline="") as f:
        candidates = list(csv.DictReader(f))
    ids = np.array([int(c["candidate_id"]) for c in candidates])
    assert len(ids) == 125 and len(set(ids)) == 125, "Expected 125 unique Dave candidates"
    ras = np.array([float(c["ra_deg"]) for c in candidates])
    decs = np.array([float(c["dec_deg"]) for c in candidates])
    candidate = {str(int(i)): {"tile_codes": [], "min_mask1_distance_arcsec_in_1p02arcsec_square": None} for i in ids}
    tiles = []
    with tarfile.open(archive, "r|gz") as tar:
        for member in tar:
            if not member.isfile() or not member.name.endswith(".fits"):
                continue
            if not member.name.startswith("cosmos_web_starmask_jwst_"):
                raise ValueError("Unexpected tar member: " + member.name)
            tile = member.name.split("jwst_")[-1].split(".fits")[0]
            f = tar.extractfile(member)
            h = read_header(f)
            if h["BITPIX"] != 8 or h["NAXIS"] != 2:
                raise ValueError("Unexpected FITS data type/shape")
            width, height = int(h["NAXIS1"]), int(h["NAXIS2"])
            x, y = position_pixels(h, ras, decs)
            inside = (x >= 0) & (x < width) & (y >= 0) & (y < height)
            counts = np.zeros(256, dtype=np.int64)
            candidate_vals = {}
            nearest = {str(int(ids[j])): float("inf") for j in np.flatnonzero(inside)}
            # A square radius > 1 arcsec encloses a circular 1-arcsec neighbourhood.
            pad = math.ceil(1 / (abs(h["CDELT1"]) * 3600))
            row = 0
            while row < height:
                take = min(384, height - row)
                data = f.read(take * width)
                if len(data) != take * width:
                    raise ValueError("Truncated image data: " + tile)
                array = np.frombuffer(data, dtype=np.uint8).reshape(take, width)
                counts += np.bincount(array.ravel(), minlength=256)
                center_here = inside & (y >= row) & (y < row + take)
                for j in np.flatnonzero(center_here):
                    candidate_vals[str(int(ids[j]))] = int(array[y[j] - row, x[j]])
                near_here = inside & (y + pad >= row) & (y - pad < row + take)
                for j in np.flatnonzero(near_here):
                    xx, yy = int(x[j]), int(y[j])
                    y0, y1 = max(0, yy - pad, row), min(height, yy + pad + 1, row + take)
                    x0, x1 = max(0, xx - pad), min(width, xx + pad + 1)
                    if y0 >= y1:
                        continue
                    ys, xs = np.nonzero(array[y0 - row:y1 - row, x0:x1] == 1)
                    if len(ys):
                        sep = np.min((ys + y0 - yy) ** 2 + (xs + x0 - xx) ** 2)
                        cid = str(int(ids[j]))
                        nearest[cid] = min(nearest[cid], math.sqrt(float(sep)) * abs(h["CDELT1"]) * 3600)
                row += take
            codes = {str(v): int(counts[v]) for v in np.flatnonzero(counts)}
            if set(codes) != {"0", "1"}:
                raise ValueError("Unexpected mask codes: " + str(codes))
            for cid, value in candidate_vals.items():
                candidate[cid]["tile_codes"].append({"tile": tile, "code": value})
                if math.isfinite(nearest[cid]):
                    old = candidate[cid]["min_mask1_distance_arcsec_in_1p02arcsec_square"]
                    candidate[cid]["min_mask1_distance_arcsec_in_1p02arcsec_square"] = min(old, nearest[cid]) if old is not None else nearest[cid]
            tiles.append({"tile": tile, "shape": [width, height], "pixel_scale_arcsec": abs(h["CDELT1"]) * 3600, "codes": codes})
    assert len(tiles) == 20, "Expected all twenty COSMOS-Web JWST masks"
    missing = [int(cid) for cid, d in candidate.items() if not d["tile_codes"]]
    masked = [int(cid) for cid, d in candidate.items() if any(z["code"] == 1 for z in d["tile_codes"])]
    within_one = [int(cid) for cid, d in candidate.items() if d["min_mask1_distance_arcsec_in_1p02arcsec_square"] is not None and d["min_mask1_distance_arcsec_in_1p02arcsec_square"] <= 1.0]
    result = {
        "archive_sha256": sha256(archive),
        "candidate_positions_sha256": sha256(positions),
        "n_tiles": len(tiles),
        "n_candidate_entries": len(candidates),
        "center_in_mask_ids": sorted(masked),
        "candidate_missing_tile_ids": sorted(missing),
        "candidate_within_1arcsec_of_mask_ids": sorted(within_one),
        "tile_masked_area_sum_deg2_NOT_UNION": sum(d["codes"]["1"] * (d["pixel_scale_arcsec"] / 3600) ** 2 for d in tiles),
        "tiles": tiles,
        "candidates": candidate,
        "cautions": [
            "Pixel value 1 is treated as excluded; compare to official mask semantics.",
            "Summing per-tile masked pixels double-counts overlapping sky areas and does not give the effective area.",
            "Center-pixel membership does not equal PSF/aperture/segmentation overlap.",
            "No F444W SNR, blue dropout completeness, or multi-epoch astrometric efficiency is inferred.",
            "Archive is a user-provided local attachment, not downloaded by GitHub Actions.",
        ],
    }
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", required=True)
    parser.add_argument("--positions", default=str(Path(__file__).with_name("dave125_positions.csv")))
    parser.add_argument("--output", default="jwst_starmask_audit.json")
    args = parser.parse_args()
    result = audit(args.archive, args.positions)
    Path(args.output).write_text(json.dumps(result, indent=2) + "\n")
    print(f"Tiles={result['n_tiles']}; masked centers={len(result['center_in_mask_ids'])}; "
          f"within 1 arcsec={result['candidate_within_1arcsec_of_mask_ids']}")


if __name__ == "__main__":
    main()
