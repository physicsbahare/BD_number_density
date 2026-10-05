"""Independent UNCOVER validation for the Galactic UCD count model.

This script applies the same thin+thick+halo line-of-sight integration used
for COSMOS-Web to the independent UNCOVER/A2744 field.

Two checks are intentionally kept separate:
1. F115W depth-only prediction using the repository's atmosphere-based d25
   values and Li et al. (2026) F115W=28.03 depth.
2. A late-T cross-check using the d_max values quoted by Li et al. (2026).

This is NOT yet a forward model of the Li et al. colour/SED selection.
"""

from __future__ import annotations

from pathlib import Path
import math

import numpy as np
import pandas as pd
import yaml
from scipy.integrate import quad
from scipy.stats import chi2
from astropy.coordinates import SkyCoord
import astropy.units as u

from src.cosmosweb_counts import (
    TEFF_DENSITY_2024,
    disk_rho,
    halo_rho,
    local_component_densities,
    distance_at_mag,
)

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "uncover_validation.yaml"
D25_TABLE = ROOT / "data" / "preliminary_cosmosweb_teff_counts.csv"


def field_count(
    dmax_pc: float,
    area_arcmin2: float,
    rho_local_total: float,
    ra_deg: float,
    dec_deg: float,
    h_thin_pc: float = 187.0,
    thick_fraction: float = 0.02,
    halo_fraction: float = 0.0025,
) -> float:
    """Integrate the Galactic density model in an arbitrary sightline."""
    c = SkyCoord(ra_deg * u.deg, dec_deg * u.deg, frame="icrs").galactic
    l = c.l.radian
    b = c.b.radian
    omega = (area_arcmin2 / 3600.0) * (np.pi / 180.0) ** 2

    rho_thin, rho_thick, rho_halo = local_component_densities(
        rho_local_total,
        thick_fraction=thick_fraction,
        halo_fraction=halo_fraction,
    )

    def integrand(d):
        return (
            disk_rho(d, l, b, rho_thin, 2600.0, h_thin_pc)
            + disk_rho(d, l, b, rho_thick, 3600.0, 900.0)
            + halo_rho(d, l, b, rho_halo)
        ) * d * d

    return omega * quad(integrand, 0.0, dmax_pc, epsabs=1e-9)[0]


def garwood_interval(n: int, confidence: float = 0.68) -> tuple[float, float]:
    """Exact central Poisson confidence interval for an observed count."""
    alpha = 1.0 - confidence
    lo = 0.0 if n == 0 else 0.5 * chi2.ppf(alpha / 2.0, 2 * n)
    hi = 0.5 * chi2.ppf(1.0 - alpha / 2.0, 2 * (n + 1))
    return float(lo), float(hi)


def parse_bin(label: str) -> tuple[int, int]:
    lo, hi = label.split("-")
    return int(lo), int(hi)


def main() -> None:
    cfg = yaml.safe_load(CONFIG.read_text())
    d25 = pd.read_csv(D25_TABLE)

    area = float(cfg["survey"]["effective_area_arcmin2"])
    depth = float(cfg["survey"]["f115w_limit_ab"])
    ra = float(cfg["field"]["ra_deg"])
    dec = float(cfg["field"]["dec_deg"])
    hthin = float(cfg["model"]["h_thin_pc"])
    thick = float(cfg["model"]["thick_to_thin_local"])
    halo = float(cfg["model"]["halo_to_thin_local"])

    rows = []
    for _, row in d25.iterrows():
        teff_bin = parse_bin(str(row["teff_bin_K"]))
        if teff_bin not in TEFF_DENSITY_2024:
            continue
        if teff_bin[0] < 450 or teff_bin[1] > 1500:
            continue

        rho = TEFF_DENSITY_2024[teff_bin][0]
        dmax_f115 = distance_at_mag(float(row["d25_F115W_pc"]), depth)
        n_f115 = field_count(
            dmax_f115, area, rho, ra, dec, hthin, thick, halo
        )

        label = f"{teff_bin[0]}-{teff_bin[1]}"
        paper_dmax = cfg["paper_detection_reach_pc"].get(label)
        n_paper = None
        if paper_dmax is not None:
            n_paper = field_count(
                float(paper_dmax), area, rho, ra, dec, hthin, thick, halo
            )

        rows.append(
            {
                "teff_bin_K": label,
                "rho_local_pc-3": rho,
                "dmax_F115W_pc": dmax_f115,
                "N_model_F115W": n_f115,
                "N_observed": int(cfg["observed"]["teff_counts"].get(label, 0)),
                "Li_dmax_pc": paper_dmax,
                "N_model_Li_dmax": n_paper,
            }
        )

    out = pd.DataFrame(rows)
    total_model = out["N_model_F115W"].sum()
    total_obs = int(cfg["observed"]["total_t_dwarfs"])
    total_lo, total_hi = garwood_interval(total_obs, 0.68)

    late = out[out["teff_bin_K"].isin(
        ["450-600", "600-750", "750-900", "900-1050"]
    )]
    late_model = late["N_model_F115W"].sum()
    late_model_paper = late["N_model_Li_dmax"].sum()
    late_obs = int(late["N_observed"].sum())
    late_lo95, late_hi95 = garwood_interval(late_obs, 0.95)

    print(out.to_string(index=False))
    print()
    print(f"UNCOVER area: {area:.1f} arcmin^2")
    print(f"Observed all-T count: {total_obs}")
    print(f"Depth-only model, 450-1500 K: {total_model:.3f}")
    print(f"Depth-only model surface density: {total_model/area:.4f} arcmin^-2")
    print(f"Observed surface density: {total_obs/area:.4f} arcmin^-2")
    print(
        f"Observed 68% exact Poisson count interval: "
        f"[{total_lo:.3f}, {total_hi:.3f}]"
    )
    print()
    print(f"Observed late-T (450-1050 K) count: {late_obs}")
    print(f"Depth-only late-T model: {late_model:.3f}")
    print(f"Li-dmax late-T model: {late_model_paper:.3f}")
    print(
        f"Observed late-T 95% exact Poisson interval: "
        f"[{late_lo95:.3f}, {late_hi95:.3f}]"
    )


if __name__ == "__main__":
    main()
