"""Small local FITS cutout pilot: image injection only, NOT completeness.

Usage:
    python -m injection_recovery.run_cutout_pilot \
        --manifest my_pilot.json --out my_pilot_output

The manifest must provide pre-calibrated SCI cutouts (nJy/pixel), VAR cutouts
((nJy/pixel)^2), same-pixel-grid PSFs and PER-BAND x,y positions. Never run on
whole JWST mosaics: the pilot refuses cutouts larger than 2048 x 2048.
Run the original source extraction separately on its SCI/VAR output.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np
from astropy.io import fits

from .core import JWST_BANDS, ALI_BANDS, inject_psf_njy


def _array(path: str, ext: int | str = 0):
    data, hdr = fits.getdata(path, ext=ext, header=True, memmap=True)
    return np.asarray(data, dtype=float), hdr


def run(manifest: dict, outdir: Path, seed: int = 20261109) -> dict:
    project = str(manifest.get("project", "")).lower()
    if project not in ("ali", "dave"):
        raise ValueError("project must be exactly ali or dave")
    bands = ALI_BANDS if project == "ali" else JWST_BANDS
    entries = manifest.get("bands", {})
    if not isinstance(entries, dict) or any(b not in entries for b in bands):
        raise ValueError("All required bands, including F814W for Ali, are needed")
    if not manifest.get("wcs_alignment_verified"):
        raise ValueError("Must explicitly verify per-band WCS alignment before injecting")
    if not manifest.get("image_units_verified"):
        raise ValueError("Must explicitly verify SCI nJy/pixel and VAR (nJy/pixel)^2")
    if outdir.exists():
        raise FileExistsError(f"Refusing to replace existing output directory: {outdir}")
    # Validate and prepare *all* inputs before writing any outputs.
    loaded = []
    for b in bands:
        d = entries[b]
        img, hdr = _array(d["science_path"], d.get("science_ext", 0))
        var, varhdr = _array(d["variance_path"], d.get("variance_ext", 0))
        psf, _ = _array(d["psf_path"], d.get("psf_ext", 0))
        if img.size > 2048**2:
            raise ValueError(f"{b}: full mosaics are not suitable for this cutout pilot")
        loaded.append((b, d, img, hdr, var, varhdr, psf))
    outdir.mkdir(parents=True)
    rng = np.random.default_rng(seed)
    report = {
        "status": "IMAGE_INJECTION_ONLY_NOT_RECOVERED_NOT_COMPLETENESS",
        "project": project,
        "seed": seed,
        "injection_id": str(manifest["injection_id"]),
        "bands": {},
        "warning": "Must rerun the original full photometry and selection pipeline.",
    }
    for b, d, img, hdr, var, varhdr, psf in loaded:
        science, variance = inject_psf_njy(
            img, var, psf, float(d["x"]), float(d["y"]),
            float(d["flux_njy"]), float(d["electrons_per_njy"]), rng,
        )
        path = outdir / f"{report['injection_id']}_{b}.fits"
        sci_hdu = fits.ImageHDU(science.astype("float32"), header=hdr, name="SCI")
        sci_hdu.header["BUNIT"] = "nJy/pixel"
        var_hdu = fits.ImageHDU(variance.astype("float32"), header=varhdr, name="VAR")
        var_hdu.header["BUNIT"] = "(nJy/pixel)^2"
        fits.HDUList([fits.PrimaryHDU(), sci_hdu, var_hdu]).writeto(path)
        report["bands"][b] = {
            "injected_path": str(path),
            "x_pixel": float(d["x"]),
            "y_pixel": float(d["y"]),
            "flux_njy": float(d["flux_njy"]),
        }
    (outdir / "injection_manifest.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--seed", default=20261109, type=int)
    opts = parser.parse_args()
    config = json.loads(opts.manifest.read_text())
    report = run(config, opts.out, opts.seed)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
