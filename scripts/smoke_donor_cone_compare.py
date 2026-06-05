#!/usr/bin/env python3
"""Donor-cone DoF compare: M-D-X angle RMSE before/after polish with
DELFIN_FFFREE_GRIP_DONOR_CONE=0 vs =1 (both at angle-to-metal=1).

Mirrors scripts/smoke_polish_angle_compare.py which validated the
angle-to-metal layer (commit ebb8cc3) on the 2792332 pool.  Reusing its
mol-from-xyz-block builder keeps the comparison apples-to-apples with
the existing GRIP test harness.

Reports the standard 4 quantities:

    rmse_pre       initial M-D-X RMSE
    rmse_off       polished with DONOR_CONE=0
    rmse_on        polished with DONOR_CONE=1
    delta_off      rmse_off - rmse_pre  (baseline polish drift)
    delta_on       rmse_on  - rmse_pre  (new operator drift)

If the donor-cone DoF helps, ``delta_on`` is more negative than
``delta_off`` (the polish reduces M-D-X deviation more).

USAGE:
    PYTHONHASHSEED=0 DELFIN_GRIP_LIB_PATH=/path/to/grip_lib_v4.npz \\
        python scripts/smoke_donor_cone_compare.py \\
            --archive .../2792332-... --n 300 --parallel 40 \\
            --out /tmp/cone_smoke.json
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import math
import os
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

os.environ.setdefault("PYTHONHASHSEED", "0")

# Suppress noisy rdkit warnings -- the covalent fallback path is expected
# to trigger for TM complexes.
from rdkit import RDLogger  # noqa: E402
RDLogger.DisableLog("rdApp.*")

from rdkit import Chem  # noqa: E402
from rdkit.Chem import RWMol, BondType  # noqa: E402


COV = {
    "H":0.31,"C":0.76,"N":0.71,"O":0.66,"F":0.57,"P":1.07,"S":1.05,
    "Cl":1.02,"Br":1.20,"I":1.39,"Se":1.20,"As":1.19,
    "Sc":1.70,"Ti":1.60,"V":1.53,"Cr":1.39,"Mn":1.50,"Fe":1.42,
    "Co":1.38,"Ni":1.24,"Cu":1.32,"Zn":1.22,"Y":1.90,"Zr":1.75,
    "Nb":1.64,"Mo":1.54,"Tc":1.47,"Ru":1.46,"Rh":1.42,"Pd":1.39,"Ag":1.45,
    "Cd":1.44,"La":2.07,"Hf":1.75,"Ta":1.70,"W":1.62,"Re":1.51,"Os":1.44,
    "Ir":1.41,"Pt":1.36,"Au":1.36,"Hg":1.32,
    "Ce":2.04,"Pr":2.03,"Nd":2.01,"Sm":1.98,"Eu":1.98,"Gd":1.96,
    "Tb":1.94,"Dy":1.92,"Ho":1.92,"Er":1.89,"Tm":1.90,"Yb":1.87,"Lu":1.87,
    "Th":2.06,"U":1.96,
}


TM_SET = frozenset({
    "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
    "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
    "La","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
    "Ce","Pr","Nd","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
    "Th","U",
})


def parse_xyz_first_frame_block(path: Path) -> str:
    lines = path.read_text().splitlines()
    n = int(lines[0].strip())
    return "\n".join(lines[: n + 2])


def parse_block(block: str):
    lines = block.splitlines()
    n = int(lines[0].strip())
    syms = []; xyz = []
    for line in lines[2:2 + n]:
        parts = line.split()
        syms.append(parts[0])
        xyz.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return syms, np.asarray(xyz, dtype=np.float64)


def mol_from_xyz_block(block: str):
    """Same builder as smoke_polish_angle_compare.py (validated by ebb8cc3)."""
    syms, P = parse_block(block)
    n = len(syms)
    if n == 0:
        return None, None
    mol = RWMol()
    for s in syms:
        mol.AddAtom(Chem.Atom(s))
    bonds_added = set()
    for i in range(n):
        ri = P[i]; si = syms[i]; ci = COV.get(si, 0.7)
        for j in range(i + 1, n):
            rj = P[j]; sj = syms[j]; cj = COV.get(sj, 0.7)
            cut = 1.25 * (ci + cj)
            if np.linalg.norm(ri - rj) <= cut:
                if (i, j) not in bonds_added:
                    mol.AddBond(i, j, BondType.SINGLE)
                    bonds_added.add((i, j))
    mol_final = mol.GetMol()
    conf = Chem.Conformer(n)
    for i in range(n):
        conf.SetAtomPosition(i, P[i].tolist())
    mol_final.AddConformer(conf, assignId=True)
    try:
        Chem.SanitizeMol(
            mol_final,
            sanitizeOps=(Chem.SanitizeFlags.SANITIZE_FINDRADICALS
                         | Chem.SanitizeFlags.SANITIZE_KEKULIZE
                         | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
                         | Chem.SanitizeFlags.SANITIZE_SETCONJUGATION
                         | Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION),
            catchErrors=True,
        )
    except Exception:
        pass
    return mol_final, P


def find_metal_donors(mol):
    metal = -1
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in TM_SET:
            metal = int(atom.GetIdx())
            break
    if metal < 0:
        return -1, []
    metal_atom = mol.GetAtomWithIdx(metal)
    donors = sorted(int(n.GetIdx()) for n in metal_atom.GetNeighbors()
                    if n.GetSymbol() not in TM_SET and n.GetSymbol() != "H")
    return metal, donors


def measure_mdx_rmse(mol, P, metal, donors):
    from delfin.fffree.grip_fragment_detect import build_mdx_angle_terms
    if metal < 0 or not donors:
        return None
    terms = build_mdx_angle_terms(mol, metal, donors, include_hydrogens=False)
    if not terms:
        return None
    devs_sq: List[float] = []; devs: List[float] = []
    for t in terms:
        u = P[t.a] - P[t.b]; v = P[t.c] - P[t.b]
        nu = float(np.linalg.norm(u)); nv = float(np.linalg.norm(v))
        if nu < 1e-9 or nv < 1e-9:
            continue
        cos_t = max(-1.0, min(1.0, float(np.dot(u, v) / (nu * nv))))
        th = math.degrees(math.acos(cos_t))
        devs_sq.append((th - t.mu) ** 2)
        devs.append(abs(th - t.mu))
    if not devs_sq:
        return None
    return {
        "n": len(devs),
        "rmse": float(math.sqrt(np.mean(devs_sq))),
        "max_abs": float(np.max(devs)),
    }


def _process_one(args):
    fpath_str, lib_path = args
    fpath = Path(fpath_str)
    os.environ["PYTHONHASHSEED"] = "0"
    if lib_path:
        os.environ["DELFIN_GRIP_LIB_PATH"] = lib_path
    sys.path.insert(0, str(REPO_ROOT))
    try:
        block = parse_xyz_first_frame_block(fpath)
        mol, P0 = mol_from_xyz_block(block)
        if mol is None:
            return {"file": fpath.name, "skip": "mol-build"}
    except Exception as exc:
        return {"file": fpath.name, "error": f"parse: {exc!r}"}
    metal, donors = find_metal_donors(mol)
    if metal < 0 or not donors:
        return {"file": fpath.name, "skip": "no metal/donors"}
    pre = measure_mdx_rmse(mol, P0, metal, donors)
    if pre is None:
        return {"file": fpath.name, "skip": "no M-D-X"}
    # Force angle-to-metal ON for both runs.
    os.environ["DELFIN_FFFREE_GRIP_ANGLE_TO_METAL"] = "1"
    from delfin.fffree.grip_polish import grip_polish

    # OFF run.
    os.environ.pop("DELFIN_FFFREE_GRIP_DONOR_CONE", None)
    try:
        P_off = grip_polish(P0, mol, metal, donors, geom="")
        post_off = measure_mdx_rmse(mol, np.asarray(P_off), metal, donors)
    except Exception as exc:
        return {"file": fpath.name, "error_off": str(exc)}

    # ON run.
    os.environ["DELFIN_FFFREE_GRIP_DONOR_CONE"] = "1"
    try:
        P_on = grip_polish(P0, mol, metal, donors, geom="")
        post_on = measure_mdx_rmse(mol, np.asarray(P_on), metal, donors)
    except Exception as exc:
        return {"file": fpath.name, "error_on": str(exc)}
    finally:
        os.environ.pop("DELFIN_FFFREE_GRIP_DONOR_CONE", None)

    return {
        "file": fpath.name,
        "metal": mol.GetAtomWithIdx(metal).GetSymbol(),
        "n_donors": len(donors),
        "n_mdx": pre["n"],
        "rmse_pre": pre["rmse"],
        "rmse_off": post_off["rmse"] if post_off else None,
        "rmse_on":  post_on["rmse"] if post_on else None,
        "max_pre": pre["max_abs"],
        "max_off": post_off["max_abs"] if post_off else None,
        "max_on":  post_on["max_abs"] if post_on else None,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--archive", required=True)
    ap.add_argument("--n", type=int, default=300)
    ap.add_argument("--out", required=True)
    ap.add_argument("--parallel", type=int, default=1)
    ap.add_argument("--lib", default=os.environ.get(
        "DELFIN_GRIP_LIB_PATH", ""))
    args = ap.parse_args()
    arc = Path(args.archive)
    files = sorted(p for p in arc.iterdir() if p.suffix == ".xyz")[: args.n]
    print(f"[smoke] {len(files)} files, parallel={args.parallel}",
          file=sys.stderr)
    t0 = time.time()
    tasks = [(str(f), args.lib) for f in files]
    results: List[Dict] = []
    if args.parallel <= 1:
        for i, t in enumerate(tasks):
            results.append(_process_one(t))
            if (i + 1) % 20 == 0:
                print(f"  {i+1}/{len(tasks)} "
                      f"t={time.time()-t0:.1f}s", file=sys.stderr)
    else:
        with cf.ProcessPoolExecutor(max_workers=args.parallel) as ex:
            for i, r in enumerate(ex.map(_process_one, tasks, chunksize=1)):
                results.append(r)
                if (i + 1) % 20 == 0:
                    print(f"  {i+1}/{len(tasks)} "
                          f"t={time.time()-t0:.1f}s", file=sys.stderr)
    valid = [r for r in results
             if "rmse_off" in r and "rmse_on" in r
             and r.get("rmse_off") is not None
             and r.get("rmse_on") is not None]
    summary = {
        "n_files": len(files),
        "n_valid": len(valid),
        "mean_rmse_pre": float(np.mean([r["rmse_pre"] for r in valid])) if valid else None,
        "mean_rmse_off": float(np.mean([r["rmse_off"] for r in valid])) if valid else None,
        "mean_rmse_on":  float(np.mean([r["rmse_on"]  for r in valid])) if valid else None,
        "mean_delta_off": float(np.mean([r["rmse_off"] - r["rmse_pre"] for r in valid])) if valid else None,
        "mean_delta_on":  float(np.mean([r["rmse_on"]  - r["rmse_pre"] for r in valid])) if valid else None,
        "mean_max_pre": float(np.mean([r["max_pre"] for r in valid])) if valid else None,
        "mean_max_off": float(np.mean([r["max_off"] for r in valid])) if valid else None,
        "mean_max_on":  float(np.mean([r["max_on"]  for r in valid])) if valid else None,
        "n_files_improved_on_vs_off": sum(1 for r in valid if r["rmse_on"] < r["rmse_off"]),
        "n_files_worsened_on_vs_off": sum(1 for r in valid if r["rmse_on"] > r["rmse_off"]),
        "n_files_equal_on_vs_off":    sum(1 for r in valid if r["rmse_on"] == r["rmse_off"]),
    }
    Path(args.out).write_text(json.dumps(
        {"summary": summary, "per_file": results}, indent=2,
    ))
    print(json.dumps(summary, indent=2))
    print(f"[smoke] wrote {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
