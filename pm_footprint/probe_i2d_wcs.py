"""Compare corner-derived WHT astrometry against the original DR1 i2d SCI WCS.

Downloads only the first 4 MiB of the gzip stream for each official i2d
tile; NEVER downloads full multi-GB mosaics or image pixel data.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import html
from html.parser import HTMLParser
import json
from pathlib import Path
import zlib
from urllib.parse import urljoin
from urllib.request import Request, urlopen

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from .wht_depth_mask import official_corner_wcs, TILES

RELEASE = "https://cosmos2025.iap.fr/nircam.html"


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.hrefs = []

    def handle_starttag(self, tag, attrs):
        if tag.lower() == "a":
            url = dict(attrs).get("href")
            if url:
                self.hrefs.append(html.unescape(url))


def find_i2d_url(tile, page=None):
    if tile not in TILES:
        raise ValueError(f"Not a DR1 tile: {tile}")
    if page is None:
        with urlopen(Request(RELEASE, headers={"User-Agent": "COSMOS-Web-science-WCS-audit"}),
                     timeout=60) as response:
            page = response.read().decode("utf-8")
    links = Links()
    links.feed(page)
    stem = f"mosaic_nircam_f444w_COSMOS-Web_60mas_{tile}_v1.0_i2d"
    matches = [urljoin(RELEASE, href) for href in links.hrefs
               if stem.lower() in href.lower()]
    if len(matches) != 1:
        raise RuntimeError(f"Expected one official i2d link for {tile}; got {matches}")
    return matches[0]


def header_from_prefix(data, offset):
    """Read a FITS header without touching any image data."""
    for pos in range(offset, len(data) - 80, 80):
        if data[pos:pos + 8] == b"END     ":
            end = ((pos + 80 + 2879) // 2880) * 2880
            if end > len(data):
                break
            return fits.Header.fromstring(data[offset:end].decode("ascii")), end
    raise ValueError("FITS header END not present in compressed prefix")


def i2d_sci_header(url, prefix_bytes=4 * 1024 * 1024):
    request = Request(url, headers={
        "Range": f"bytes=0-{prefix_bytes-1}",
        "User-Agent": "COSMOS-Web-science-WCS-audit",
    })
    with urlopen(request, timeout=120) as response:
        compressed = response.read(prefix_bytes)
        code = response.status
    if len(compressed) < 1024 or compressed[:2] != b"\x1f\x8b":
        raise RuntimeError(f"Official i2d prefix is not gzip: HTTP {code}")
    try:
        raw = zlib.decompressobj(16 + zlib.MAX_WBITS).decompress(
            compressed, max_length=4 * 1024 * 1024)
    except zlib.error as exc:
        raise RuntimeError(f"Invalid i2d gzip prefix: {exc}") from exc
    primary, offset = header_from_prefix(raw, 0)
    if primary.get("NAXIS", 0) != 0:
        raise RuntimeError("Unexpected non-empty i2d primary: cannot skip image data")
    sci, _ = header_from_prefix(raw, offset)
    if sci.get("EXTNAME", "").upper() != "SCI":
        raise RuntimeError(f"Expected i2d SCI extension; got {sci.get('EXTNAME')}")
    wcs = WCS(sci, relax=True).celestial
    if not wcs.has_celestial:
        raise RuntimeError("Original i2d SCI header lacks celestial WCS")
    return sci, wcs, {"http_status": code, "compressed_prefix_bytes": len(compressed),
                      "decompressed_prefix_bytes": len(raw)}


def compare(tile, url):
    sci, true_wcs, metadata = i2d_sci_header(url)
    nx, ny = int(sci["NAXIS1"]), int(sci["NAXIS2"])
    inferred, fit = official_corner_wcs(tile, nx, ny)
    pixels = np.array([[0., 0.], [nx - 1., 0.],
                       [nx - 1., ny - 1.], [0., ny - 1.],
                       [(nx - 1.) / 2, (ny - 1.) / 2]])
    ra_true, dec_true = true_wcs.all_pix2world(pixels[:, 0], pixels[:, 1], 0)
    ra_fit, dec_fit = inferred.all_pix2world(pixels[:, 0], pixels[:, 1], 0)
    sep = np.hypot((ra_true - ra_fit) * np.cos(np.deg2rad(dec_true)),
                   dec_true - dec_fit) * 3600
    result = {
        "tile": tile, "official_i2d_url": url,
        "image_shape": [ny, nx],
        "sci_ctype": [sci.get("CTYPE1"), sci.get("CTYPE2")],
        "corner_vs_original_sci_wcs_arcsec": [float(x) for x in sep],
        "max_wcs_difference_arcsec": float(sep.max()),
        "corner_fit": fit, **metadata,
    }
    result["status"] = ("PASS_WITHIN_0P6ARCSEC" if sep.max() <= 0.6
                        else "FAIL_INDEPENDENT_WCS_MISMATCH")
    return result


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tiles", nargs="+", default=["A1", "B10"])
    p.add_argument("--out", default="pm_footprint/results/wcs_probe.json")
    a = p.parse_args()
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    results = []
    try:
        with urlopen(Request(RELEASE, headers={"User-Agent": "COSMOS-Web-science-WCS-audit"}),
                     timeout=60) as response:
            page = response.read().decode("utf-8")
        for tile in a.tiles:
            result = compare(tile, find_i2d_url(tile, page))
            results.append(result)
            print(json.dumps(result, indent=2), flush=True)
        output = {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "source": RELEASE,
            "status": ("PASS" if all(x["status"].startswith("PASS") for x in results)
                       else "FAIL_REQUIRES_REVIEW"),
            "results": results,
        }
        out.write_text(json.dumps(output, indent=2) + "\n")
        if output["status"] != "PASS":
            raise RuntimeError("Original SCI WCS differs from corner reconstruction")
    except Exception as exc:
        (out.parent / "wcs_probe_failure.json").write_text(json.dumps({
            "status": "FAILED_NO_ASTROMETRY_VALIDATION",
            "error": repr(exc), "completed_tiles": results,
        }, indent=2) + "\n")
        raise


if __name__ == "__main__":
    main()
