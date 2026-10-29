"""Compare independent depth-14 and depth-15 MOC geometry runs.

This only validates HEALPix rasterization stability, NOT astrophysical PM
detectability or agreement with the COSMOS-Web usable-depth mask.
"""
import json
from pathlib import Path

ROOT = Path("pm_footprint/results")


def main():
    coarse = json.loads((ROOT / "summary.json").read_text())
    fine = json.loads((ROOT / "depth15/summary.json").read_text())
    if coarse["moc_max_depth"] != 14 or fine["moc_max_depth"] != 15:
        raise ValueError("Expected order-14 and order-15 summaries")
    if coarse["mast_rows"] != fine["mast_rows"]:
        raise ValueError("Different MAST row counts; comparison invalid")
    for field in ("compute_source_sha256", "mast_metadata_sha256",
                  "candidate_catalog_sha256", "baseline_days"):
        if field not in coarse or coarse[field] != fine.get(field):
            raise ValueError(f"Different run provenance: {field}")
    deltas = {}
    membership = {}
    for key in coarse["area_deg2"]:
        a = coarse["area_deg2"][key]
        b = fine["area_deg2"][key]
        deltas[key] = {
            "order14_deg2": a,
            "order15_deg2": b,
            "absolute_delta_deg2": b - a,
            "relative_delta": (b - a) / a if a else None,
        }
        ca = coarse["candidates_with_geometric_coverage"][key]
        cb = fine["candidates_with_geometric_coverage"][key]
        membership[key] = {"order14": ca, "order15": cb, "difference": cb - ca}
    unstable_area = [
        key for key, v in deltas.items()
        if abs(v["absolute_delta_deg2"]) > max(0.001, 0.02 * v["order15_deg2"])
    ]
    unstable_membership = [
        key for key, v in membership.items() if v["difference"] != 0
    ]
    report = {
        "compute_source_sha256": coarse["compute_source_sha256"],
        "mast_metadata_sha256": coarse["mast_metadata_sha256"],
        "candidate_catalog_sha256": coarse["candidate_catalog_sha256"],
        "baseline_days": coarse["baseline_days"],
        "status": ("RASTERIZATION_STABLE_REQUIRES_SCIENCE_QC"
                   if not unstable_area and not unstable_membership
                   else "BOUNDARY_OR_RASTERIZATION_REVIEW_NEEDED"),
        "area_comparison": deltas,
        "candidate_membership_comparison": membership,
        "unstable_area_masks": unstable_area,
        "unstable_membership_masks": unstable_membership,
        "notes": [
            "Stable HEALPix resolution does not establish survey usable-depth coverage.",
            "Same-filter F444W is the primary conservative imaging opportunity.",
            "Cross-filter geometrical coverage does not imply PM sensitivity.",
            "Check detector-level footprints and 30-day baseline before quoting."
        ],
    }
    (ROOT / "resolution_qc.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
