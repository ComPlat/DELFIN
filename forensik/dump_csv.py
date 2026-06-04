#!/usr/bin/env python3
"""Phase F deliverable: per-SMILES isocoverage CSV for paper supplementary."""
import csv
import json
from pathlib import Path

JSONL = Path("/home/qmchem_max/agent_workspace/quality_framework/xyz_archive/b00f9a0-full7-VOLLPOOL_isocoverage.jsonl")
OUT = Path("/home/qmchem_max/ComPlat/DELFIN/.claude/worktrees/agent-a6a0c001c3fc98a82/paper_data/isocov_per_smiles_b00f9a0.csv")
OUT.parent.mkdir(parents=True, exist_ok=True)

with JSONL.open() as f, OUT.open("w", newline="") as g:
    w = csv.writer(g)
    w.writerow([
        "smiles_id", "class", "metal", "metal_idx", "cn",
        "polyhedron_class", "theoretical_geoms",
        "n_theoretical_isomers", "n_named_observed",
        "coverage_pct", "n_frames",
        "n_frames_unknown_geometry", "n_frames_polyhedron_mismatch",
        "n_named_missing", "n_named_missing_cn_mismatch",
        "smiles",
    ])
    n = 0
    for line in f:
        r = json.loads(line)
        if r.get("error"):
            continue
        w.writerow([
            r.get("smiles_id", ""),
            r.get("class", ""),
            r.get("metal", ""),
            r.get("metal_idx", ""),
            r.get("cn", ""),
            r.get("polyhedron_class", ""),
            "|".join(r.get("theoretical_geoms", []) or []),
            len(r.get("theoretical_isomers", []) or []),
            len(r.get("named_observed", []) or []),
            r.get("coverage_pct", 0),
            r.get("n_frames", 0),
            r.get("n_frames_unknown_geometry", 0),
            r.get("n_frames_polyhedron_mismatch", 0),
            len(r.get("named_missing", []) or []),
            len(r.get("named_missing_cn_mismatch", []) or []),
            r.get("smiles", ""),
        ])
        n += 1
print(f"Wrote {n} rows to {OUT}")
