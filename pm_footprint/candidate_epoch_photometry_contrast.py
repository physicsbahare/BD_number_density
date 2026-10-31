"""F444W photometry versus >=365d F356W geometric coverage, 125 Dave candidates.

Reproducible from small committed CSV; no FITS, MAST, or network needed.
Statistical tests are *conditional diagnostics*, not survey completeness.
"""
import argparse
import csv
import json
import math
from pathlib import Path

DEFAULT = Path("pm_footprint/results/v13_epoch_candidate_photometry_fixture_20261031.csv")


def median(values):
    a = sorted(values)
    n = len(a)
    if not n:
        raise ValueError("Empty median")
    return a[n // 2] if n % 2 else (a[n // 2 - 1] + a[n // 2]) / 2


def mean(values):
    return sum(values) / len(values)


def groups(rows, radius_arcsec=0.5):
    parent = list(range(len(rows)))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i, a in enumerate(rows):
        for j in range(i + 1, len(rows)):
            b = rows[j]
            dec = math.radians((a["dec_deg"] + b["dec_deg"]) / 2)
            sep = 3600 * math.hypot(
                (a["ra_deg"] - b["ra_deg"]) * math.cos(dec),
                a["dec_deg"] - b["dec_deg"])
            if sep <= radius_arcsec:
                parent[root(j)] = root(i)
    out = {}
    for i, row in enumerate(rows):
        out.setdefault(root(i), []).append(row)
    return list(out.values())


def rng32(seed):
    x = seed & 0xffffffff
    while True:
        x ^= (x << 13) & 0xffffffff
        x ^= x >> 17
        x ^= (x << 5) & 0xffffffff
        x &= 0xffffffff
        yield x / 4294967296


def permutation(rows, key="mag_F444W", reps=20000, seed=20261031):
    a = [r[key] for r in rows if r["red365"]]
    b = [r[key] for r in rows if not r["red365"]]
    if not a or not b:
        raise ValueError("Both epoch groups required")
    dmed, dmean = median(a) - median(b), mean(a) - mean(b)
    values = [r[key] for r in rows]
    indices = list(range(len(values)))
    rng = rng32(seed)
    count_med = count_mean = 0
    for _ in range(reps):
        for j in range(len(indices) - 1, 0, -1):
            k = int(next(rng) * (j + 1))
            indices[j], indices[k] = indices[k], indices[j]
        x = [values[i] for i in indices[:len(a)]]
        y = [values[i] for i in indices[len(a):]]
        count_med += abs(median(x) - median(y)) >= abs(dmed) - 1e-12
        count_mean += abs(mean(x) - mean(y)) >= abs(dmean) - 1e-12
    return dict(red_n=len(a), other_n=len(b), red_median=median(a),
                other_median=median(b), delta_median=dmed, delta_mean=dmean,
                permutation_p_median=(count_med + 1) / (reps + 1),
                permutation_p_mean=(count_mean + 1) / (reps + 1),
                permutations=reps, seed=seed)


def passes(r, error_factor=1., f444_factor=1., f277_offset_njy=0.):
    f444 = r["f_F444W"] * f444_factor
    f277 = r["f_F277W"] + f277_offset_njy
    e277 = r["new_e_F277W"] * error_factor
    snr444 = f444 / (r["new_e_F444W"] * error_factor)
    blue = [r["f_F115W"] / (r["new_e_F115W"] * error_factor),
            r["f_F150W"] / (r["new_e_F150W"] * error_factor),
            f277 / e277]
    # Published rule: positive measured F277W uses its flux;
    # nonpositive measured F277W uses 2*new_e_F277W.
    ratio = f444 / (f277 if f277 > 0 else 2 * e277)
    mag = r["mag_F444W"] - 2.5 * math.log10(f444_factor)
    return (snr444 > 5 and all(s < 2 for s in blue)
            and ratio > 30 and mag < 26.2)


def analyze(rows, reps=20000):
    if len(rows) != 125 or len({r["candidate_id"] for r in rows}) != 125:
        raise ValueError("Expected frozen 125 unique IDs")
    if sum(r["red365"] for r in rows) != 62:
        raise ValueError("Expected 62 F356W >=365d geometric memberships")
    grouped = groups(rows)
    if any(len({r["red365"] for r in g}) != 1 for g in grouped):
        raise ValueError("Positional group crosses epoch coverage boundary")
    collapsed = [dict(red365=g[0]["red365"],
                      mag_F444W=mean([r["mag_F444W"] for r in g]),
                      snr444=mean([r["f_F444W"] / r["new_e_F444W"] for r in g]))
                 for g in grouped]
    scenarios = {
        "nominal": {},
        "new_e_plus_2percent": {"error_factor": 1.02},
        "new_e_plus_5percent": {"error_factor": 1.05},
        "F444W_flux_minus_2percent": {"f444_factor": 0.98},
        "F444W_flux_minus_5percent": {"f444_factor": 0.95},
        "F277W_flux_plus_1nJy": {"f277_offset_njy": 1.},
        "F277W_flux_plus_2nJy": {"f277_offset_njy": 2.},
    }
    cuts = {}
    for label, kwargs in scenarios.items():
        failed = [r["candidate_id"] for r in rows if not passes(r, **kwargs)]
        cuts[label] = {"passing": len(rows) - len(failed), "failing_ids": failed}
    return {
        "status": "SOURCE_CONDITIONED_DIAGNOSTIC_NOT_COMPLETENESS",
        "catalog_sources": len(rows),
        "positional_groups_0p5arcsec": len(grouped),
        "red_epoch_groups": sum(r["red365"] for r in collapsed),
        "source_mag": permutation(rows, reps=reps),
        "group_mag": permutation(collapsed, reps=reps),
        "group_f444w_snr": permutation(collapsed, "snr444", reps=reps),
        "selection_sensitivity": cuts,
        "input_origin": "User-supplied sky-verified 125 Zhen-Kai v13 rows joined to frozen F356W epoch list",
        "input_original_125_rows_sha256": "93400baee3daccfc5020b01f93d05e9b65ef238afa62b5b691f4aa90728e2a40",
        "caveats": [
            "125 sources were selected on F444W; this is not an independent depth measurement.",
            "Permutation assumes exchangeable locations; depth gradients and spatial clustering violate this null.",
            "0.5arcsec proximity groups are not verified unique astrophysical objects.",
            "Flux/error perturbations test stability of already-selected rows, not injection-recovery completeness.",
            "F356W coverage is geometric; neither F356W SNR nor PM precision has been measured.",
        ],
    }


def load(path):
    with path.open(newline="") as stream:
        data = list(csv.DictReader(stream))
    if not data:
        raise ValueError("Empty input")
    rows = []
    for r in data:
        x = {k: float(v) for k, v in r.items()
             if k not in ("candidate_id", "red365")}
        x["candidate_id"] = int(r["candidate_id"])
        x["red365"] = r["red365"] == "1"
        rows.append(x)
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, default=DEFAULT)
    p.add_argument("--out", type=Path, default=Path(
        "pm_footprint/results/candidate_epoch_photometry_contrast_20261031.json"))
    p.add_argument("--permutations", type=int, default=20000)
    args = p.parse_args()
    if args.permutations < 0:
        p.error("Permutations must be nonnegative")
    result = analyze(load(args.input), reps=args.permutations)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print("red/other median F444W:", result["source_mag"]["red_median"],
          result["source_mag"]["other_median"])


if __name__ == "__main__":
    main()
