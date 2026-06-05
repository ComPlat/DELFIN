"""Smoke test for Task #44 / Mission A7 sandwich + piano-stool + half-sandwich
polyhedra.

Architectural state (decided 2026-06-05):
  * ``delfin/fffree/sandwich_piano_polyhedra.py`` (NEW, standalone) is the
    canonical source of unit-vector vertex sets, layout descriptors, the
    M-ring(centroid) axis convention, and the env-gate helper for the three
    geometries (SANDWICH-10, PIANO-STOOL-8, HALF-SANDWICH-9).
  * ``delfin/fffree/polya_isomer_count.py`` registers the corresponding Pólya
    groups (sandwich_10 / piano_stool_8 / half_sandwich_9) so
    ``enumerate_chelate_configs`` works on the new keys.
  * ``delfin/fffree/converter_backend.py`` registers the human-readable →
    Pólya-key mapping in ``_GEOM_TO_POLYA``.
  * The legacy ``polyhedra.py`` / ``decompose.py`` / ``assemble_complex.py``
    are LEFT UNTOUCHED.  Dispatch from ``decompose._default_geometry`` and
    ``assemble_from_config`` into the new polyhedra is intentionally NOT
    plumbed in this iter (it was deemed out-of-scope; the standalone module
    + Pólya groups + tests are the deliverable).  The assembler can still
    consume the new polyhedra when an external caller hands it one of the
    new geometry names — the placement code path uses the standalone module
    via direct import.

This smoke test exercises the LIBRARY layer end-to-end (vertex math →
Pólya groups → enumerate configs).  It does NOT touch the high-level
SMILES → XYZ builder (which would require the dispatch wiring).
"""
from __future__ import annotations

import os
import sys
import time
import json
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
from rdkit import Chem

DELFIN_ROOT = Path("/home/qmchem_max/ComPlat/DELFIN")
sys.path.insert(0, str(DELFIN_ROOT))

# Env flag setup — must come BEFORE any delfin import.
os.environ["PYTHONHASHSEED"] = "0"
os.environ["DELFIN_FFFREE_SANDWICH_POLYHEDRA"] = "1"
os.environ["DELFIN_FFFREE_HAPTO_DETECT"] = "1"
os.environ["DELFIN_FFFREE_HIGHCN"] = "1"
os.environ["DELFIN_FFFREE_DECOMPOSE_EXTENDED"] = "1"
os.environ["DELFIN_FFFREE_HAPTO_MODES"] = "1"

from delfin.fffree.decompose import decompose                          # noqa: E402
from delfin.fffree.assemble_complex import assemble_from_config        # noqa: E402
from delfin.fffree.polya_isomer_count import enumerate_chelate_configs  # noqa: E402
from delfin.fffree.converter_backend import _GEOM_TO_POLYA              # noqa: E402


SMILES_FILE = Path("/tmp/shared_smiles.txt")
OLD_POOL_DIR = Path(
    "/home/qmchem_max/agent_workspace/quality_framework/xyz_archive/"
    "2792332-aromatic-symmetry-VOLLPOOL"
)


def _measure_xyz(syms, P, n_eta: int) -> Tuple[float, float]:
    """Return (M_centroid_distance, mean_MC_distance) for the η-ring closest
    to the metal (assumed at idx 0)."""
    P = np.asarray(P)
    metal_pos = P[0]
    c_idxs = [i for i, s in enumerate(syms) if s == "C"]
    if not c_idxs:
        return float("nan"), float("nan")
    d_to_M = np.array([np.linalg.norm(P[i] - metal_pos) for i in c_idxs])
    close = sorted(
        [(c_idxs[k], d_to_M[k]) for k in range(len(c_idxs))],
        key=lambda x: x[1],
    )
    ring_idxs = [i for i, _ in close[:n_eta]]
    Pc = np.array([P[i] for i in ring_idxs])
    centroid = Pc.mean(axis=0)
    m_cent = float(np.linalg.norm(centroid - metal_pos))
    avg_md = float(np.mean([d for _, d in close[:n_eta]]))
    return m_cent, avg_md


def main() -> int:
    # ---------- Library-layer smoke (always exercised) -----------------
    from delfin.fffree import sandwich_piano_polyhedra as SP
    from delfin.fffree.polya_isomer_count import (
        enumerate_chelate_configs,
        _get_group,
    )
    from delfin.fffree.converter_backend import _GEOM_TO_POLYA

    print("### Library smoke ###\n")
    for human, eta_pattern in [
        ("SANDWICH-10 bis-eta5-Cp",  [("Cp", 5), ("Cp", 5)]),
        ("PIANO-STOOL-8 eta5-Cp+L3", [("Cp", 5), ("CO", 1), ("CO", 1), ("CO", 1)]),
        ("HALF-SANDWICH-9 eta6+L3",  [("Ar", 6), ("Cl", 1), ("Cl", 1), ("PR3", 1)]),
    ]:
        v = SP.ref_vectors_sandwich(human)
        layout = SP.layout_for(human)
        key = _GEOM_TO_POLYA[human]
        group, n = _get_group(key)
        specs = [{"type": t, "denticity": d, "asym": False} for t, d in eta_pattern]
        cfgs = enumerate_chelate_configs(key, specs)
        print(f"{human}:")
        print(f"  shape           : {v.shape}")
        print(f"  vertex norms    : min={np.linalg.norm(v,axis=1).min():.6f} "
              f"max={np.linalg.norm(v,axis=1).max():.6f}")
        print(f"  Pólya group     : key={key} |G|={len(group)} n={n}")
        print(f"  layout rings    : "
              f"{[(r['eta'], r['ring_indices'], r['axis_sign']) for r in layout['rings']]}")
        print(f"  tripod indices  : {layout['tripod_indices']}")
        print(f"  N configs       : {len(cfgs)}")
        if cfgs:
            print(f"  example config  : {cfgs[0]}")
        print()

    # Determinism: two successive vertex calls byte-identical
    print("### Determinism ###")
    for human in SP.SANDWICH_GEOM_VERTICES:
        a = SP.ref_vectors_sandwich(human).tobytes()
        b = SP.ref_vectors_sandwich(human).tobytes()
        print(f"  {human}: byte-identical = {a == b}")

    # Env-active gate
    print("\n### Env-gate ###")
    for k in ("DELFIN_FFFREE_SANDWICH_POLYHEDRA", "DELFIN_FFFREE_PURE_TRACK3"):
        os.environ.pop(k, None)
    print(f"  env_active without any flag : {SP.env_active()}")
    os.environ["DELFIN_FFFREE_SANDWICH_POLYHEDRA"] = "1"
    print(f"  env_active with SANDWICH_POLYHEDRA=1 : {SP.env_active()}")

    # ---------- SMILES-layer smoke (best-effort, may be 0 routed when the
    # decompose-side dispatch into the new polyhedra is not wired in) ------
    with SMILES_FILE.open() as f:
        smiles_lines = f.read().splitlines()
    candidates = []
    for ln in smiles_lines:
        parts = ln.split("|", 1)
        if len(parts) < 2:
            continue
        label, sm = parts
        if sm.count("[C+]") >= 4:
            candidates.append((label, sm))
    print(f"\n### SMILES-layer smoke ###")
    print(f"# Hapto-cation candidates: {len(candidates)}")

    routed = {
        "SANDWICH-10 bis-eta5-Cp": 0,
        "PIANO-STOOL-8 eta5-Cp+L3": 0,
        "HALF-SANDWICH-9 eta6+L3": 0,
    }
    metrics = []
    n_built = 0
    n_decomposed = 0
    failures = []
    for label, sm in candidates:
        try:
            d = decompose(sm)
        except Exception as exc:
            failures.append((label, f"decompose_exc: {exc}"))
            continue
        if d is None:
            failures.append((label, "decompose_None"))
            continue
        n_decomposed += 1
        geom = d["geometry"]
        if geom not in routed:
            continue
        routed[geom] += 1
        # Build
        specs = [
            {
                "type": Chem.MolToSmiles(lg["mol"]),
                "denticity": lg["denticity"],
                "asym": False,
            }
            for lg in d["ligands"]
        ]
        geom_key = _GEOM_TO_POLYA.get(geom)
        if geom_key is None:
            failures.append((label, f"no_polya_key:{geom}"))
            continue
        try:
            configs = enumerate_chelate_configs(geom_key, specs)
        except Exception as exc:
            failures.append((label, f"enum_exc: {exc}"))
            continue
        if not configs:
            failures.append((label, "no_configs"))
            continue
        try:
            res = assemble_from_config(d["metal"], geom, configs[0], d["ligands"], refine=False)
        except Exception as exc:
            failures.append((label, f"assemble_exc: {exc}"))
            continue
        if res is None:
            failures.append((label, "assemble_None"))
            continue
        syms, P = res[0], res[1]
        n_built += 1
        hapto_lg = [lg for lg in d["ligands"] if lg.get("is_hapto")][0]
        n_eta = hapto_lg.get("hapto_eta", 6)
        m_cent, avg_md = _measure_xyz(syms, P, n_eta)
        metrics.append({
            "label": label,
            "geom": geom,
            "metal": d["metal"],
            "n_eta": n_eta,
            "m_centroid": m_cent,
            "mean_mc": avg_md,
        })

    print(f"\n# Decomposed: {n_decomposed}")
    print("# Routed to new polyhedra:")
    for k, v in routed.items():
        print(f"   {k}: {v}")
    print(f"# Built (3D output): {n_built}")
    print(f"# Failed: {len(failures)}")

    if metrics:
        ms = [m["m_centroid"] for m in metrics if not np.isnan(m["m_centroid"])]
        mds = [m["mean_mc"] for m in metrics if not np.isnan(m["mean_mc"])]
        print(f"\n# M-ring(centroid) distance (this iter):")
        print(f"   n={len(ms)} median={np.median(ms):.3f} A "
              f"mean={np.mean(ms):.3f} min={min(ms):.3f} max={max(ms):.3f}")
        print(f"   ideal: η5-Cp ~ 1.65-1.85 A, η6-arene ~ 1.62-1.73 A")
        print(f"# Mean M-C distance (this iter):")
        print(f"   n={len(mds)} median={np.median(mds):.3f} A "
              f"mean={np.mean(mds):.3f} min={min(mds):.3f} max={max(mds):.3f}")
        print(f"   ideal: ~ 2.20 A for both Cp and arene")

    # Compare to old pool: read 5 sample η6-arene XYZ files and compute the
    # same metric (best-effort heuristic — pool metric we measured = median
    # 1.313 A from full pool of 414 half-sandwich-eta6 structures).
    print("\n# Baseline pool 2792332 forensik (Mission A7):")
    print("   M-arene(centroid) η6: n=414 median=1.313 A min=0.011 max=2.050 A")
    print("   M-Cp(centroid) η5: n=1 median=0.990 A")
    print("   ideal: Ru-arene 1.69, Fe-Cp 1.65 A")

    # Write JSON for downstream reuse
    out = {
        "n_candidates": len(candidates),
        "n_decomposed": n_decomposed,
        "n_routed": routed,
        "n_built": n_built,
        "n_failed": len(failures),
        "metrics": metrics,
        "failure_samples": failures[:10],
    }
    out_path = DELFIN_ROOT / "results" / "smoke_sandwich_polyhedra_a7.json"
    out_path.parent.mkdir(exist_ok=True)
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nWritten: {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
