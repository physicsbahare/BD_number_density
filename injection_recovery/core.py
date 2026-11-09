"""Image-injection primitives and paper-specific photometric selections.

NOT an end-to-end completeness measurement. These functions must be integrated
with the actual COSMOS-Web forced-photometry, deblending, masking and catalogue
extraction pipeline before any recovered fractions are reported.

All science arrays must already be calibrated in nJy/pixel, variance in
(nJy/pixel)^2, and PSFs must be sampled on the same pixel grid as the image.
"""
from __future__ import annotations

from collections.abc import Mapping
import numpy as np
from scipy.ndimage import shift

JWST_BANDS = ("F115W", "F150W", "F277W", "F444W")
ALI_BANDS = ("F814W",) + JWST_BANDS


def _valid_photometry(flux: Mapping[str, float],
                      error: Mapping[str, float], bands: tuple[str, ...]) -> bool:
    """Do not treat missing/invalid measurements as nondetections."""
    try:
        return all(
            np.isfinite(float(flux[b])) and np.isfinite(float(error[b]))
            and float(error[b]) > 0
            for b in bands
        )
    except (KeyError, TypeError, ValueError):
        return False


def dave_photo_selection(flux: Mapping[str, float],
                         error: Mapping[str, float],
                         m_f444w_ab: float) -> bool:
    """Published Dave125 photometric cuts, not the full catalogue pipeline.

    Pass the *corrected* new_e_* errors. Negative S/N is signed, never abs().
    Positive F277W is used as measured; nonpositive F277W uses a 2-sigma floor.
    Does not assert completeness, source quality, or PM measurability.
    """
    if not _valid_photometry(flux, error, JWST_BANDS):
        return False
    if not np.isfinite(m_f444w_ab):
        return False
    sn = {b: float(flux[b]) / float(error[b]) for b in JWST_BANDS}
    e277 = float(error["F277W"])
    f277 = float(flux["F277W"])
    f277_eff = f277 if f277 > 0 else 2.0 * e277
    return bool(
        sn["F444W"] > 5.0
        and all(sn[b] < 2.0 for b in JWST_BANDS[:-1])
        and float(flux["F444W"]) / f277_eff > 30.0
        and m_f444w_ab < 26.2
    )


def ali_photo_selection(flux: Mapping[str, float],
                        error: Mapping[str, float],
                        f444w_auto: float,
                        warn_flag: int) -> bool:
    """Ali catalogue-level *photometric preselection* reconstructed from notes.

    F814W is required, in addition to the four NIRCam bands. The input
    F444W AUTO flux and warn_flag are catalogue-derived, NOT a replacement
    for the original catalogue extractor. Does not model visual rejection,
    morphology or any later paper revisions. Verify against Ali's frozen
    manuscript/candidate script before using as the official full selection.
    """
    if not _valid_photometry(flux, error, ALI_BANDS):
        return False
    try:
        if not np.isfinite(float(f444w_auto)):
            return False
        if int(warn_flag) != 0 or float(f444w_auto) <= 0:
            return False
    except (TypeError, ValueError):
        return False
    sn = {b: float(flux[b]) / float(error[b]) for b in ALI_BANDS}
    return bool(
        sn["F444W"] > 5.0
        and all(sn[b] < 2.0 for b in ALI_BANDS if b != "F444W")
    )


def inject_psf_njy(
    science: np.ndarray,
    variance: np.ndarray,
    psf: np.ndarray,
    x: float,
    y: float,
    flux_njy: float,
    electrons_per_njy: float,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Return noisy injected SCI, VAR copies; do NOT alter original arrays.

    (x, y) are zero-indexed *image pixel* coordinates. The normalized PSF
    is shifted subpixel using bilinear interpolation and placed at (x, y).
    Source shot noise uses a Gaussian approximation to Poisson noise. The
    pre-existing background noise and variance are preserved.

    Refuses injections crossing the image edge; callers must choose eligible
    unmasked positions and provide independently registered band-specific
    WCS coordinates. Flux-to-electron calibration is exposure-dependent.
    """
    img = np.asarray(science, dtype=float)
    var = np.asarray(variance, dtype=float)
    p = np.asarray(psf, dtype=float)
    if img.ndim != 2 or var.shape != img.shape:
        raise ValueError("SCI and VAR must be matched 2D arrays")
    if p.ndim != 2 or any(d % 2 == 0 for d in p.shape):
        raise ValueError("PSF must be an odd-sized 2D array")
    if not all(np.isfinite(z).all() for z in (img, var, p)):
        raise ValueError("SCI, VAR, PSF must be finite in the working image")
    if np.any(var < 0) or np.any(p < 0) or p.sum() <= 0:
        raise ValueError("VAR must be nonnegative and PSF nonnegative/positive sum")
    if not all(np.isfinite(z) for z in (x, y, flux_njy, electrons_per_njy)):
        raise ValueError("Injection parameters must be finite")
    if flux_njy < 0 or electrons_per_njy <= 0:
        raise ValueError("Flux >=0 and electrons_per_njy >0 are required")
    xc, yc = int(np.rint(x)), int(np.rint(y))
    hy, hx = p.shape[0] // 2, p.shape[1] // 2
    y0, y1, x0, x1 = yc-hy, yc+hy+1, xc-hx, xc+hx+1
    if y0 < 0 or x0 < 0 or y1 > img.shape[0] or x1 > img.shape[1]:
        raise ValueError("PSF footprint crosses image edge; reject this position")
    p = shift(p, shift=(y-yc, x-xc), order=1, mode="constant", cval=0.0)
    p = np.maximum(p, 0)
    p /= p.sum()
    added_flux = flux_njy * p
    added_var = added_flux / electrons_per_njy
    out, out_var = img.copy(), var.copy()
    out[y0:y1, x0:x1] += added_flux + rng.normal(
        loc=0.0, scale=np.sqrt(added_var), size=p.shape
    )
    out_var[y0:y1, x0:x1] += added_var
    return out, out_var
