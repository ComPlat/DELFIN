#!/usr/bin/env python3
"""End-to-end polish comparison: M-D-X angle RMSE before/after grip_polish
with the angle-to-metal layer OFF vs ON.

For each XYZ frame, reconstruct an RDKit Mol from the XYZ (DetermineBonds),
identify metal + donors, then run grip_polish twice (OFF / ON).  Measure
M-D-X angle deviation from the VSEPR / library priors after each polish.

The relevant deltas:

    delta_off = rmse(after_off) - rmse(before)   # baseline operator drift
    delta_on  = rmse(after_on)  - rmse(before)   # new operator drift

If the angle-to-metal layer works, ``delta_on`` should be MORE NEGATIVE
than ``delta_off`` (the polish reduces M-D-X deviation more).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path
import numpy as np

# Force env determinism for the polish
os.environ["PYTHONHASHSEED"] = "0"

from rdkit import Chem
from rdkit.Chem import rdDetermineBonds, RWMol, BondType


# Covalent radii (Å)
COV = {
    "H":0.31,"C":0.76,"N":0.71,"O":0.66,"F":0.57,"P":1.07,"S":1.05,
    "Cl":1.02,"Br":1.20,"I":1.39,"Se":1.20,"As":1.19,
    "Sc":1.70,"Ti":1.60,"V":1.53,"Cr":1.39,"Mn":1.50,"Fe":1.42,
    "Co":1.38,"Ni":1.24,"Cu":1.32,"Zn":1.22,"Y":1.90,"Zr":1.75,
    "Nb":1.64,"Mo":1.54,"Ru":1.46,"Rh":1.42,"Pd":1.39,"Ag":1.45,
    "Cd":1.44,"Hf":1.75,"Ta":1.70,"W":1.62,"Re":1.51,"Os":1.44,
    "Ir":1.41,"Pt":1.36,"Au":1.36,"Hg":1.32,
}


TM_SET = frozenset({
    "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
    "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
    "La","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
})


def parse_xyz_first_frame(path: Path):
    lines = path.read_text().splitlines()
    n = int(lines[0].strip())
    return "\n".join(lines[: n + 2])


def parse_syms_coords(block: str):
    lines = block.splitlines()
    n = int(lines[0].strip())
    syms = []; xyz = []
    for line in lines[2:2 + n]:
        parts = line.split()
        syms.append(parts[0])
        xyz.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return syms, np.array(xyz)


def mol_from_xyz_block(block: str):
    """Build a Mol manually via covalent-cutoff bond detection (works for
    TM complexes where DetermineBonds gives up)."""
    syms, P = parse_syms_coords(block)
    n = len(syms)
    if n == 0:
        return None
    mol = RWMol()
    for s in syms:
        mol.AddAtom(Chem.Atom(s))
    # Bond detection: covalent-cutoff 1.25 ×
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
    # Attach a conformer with the coordinates
    conf = Chem.Conformer(n)
    for i in range(n):
        conf.SetAtomPosition(i, P[i].tolist())
    mol_final.AddConformer(conf, assignId=True)
    try:
        Chem.SanitizeMol(mol_final, sanitizeOps=Chem.SanitizeFlags.SANITIZE_FINDRADICALS |
                                              Chem.SanitizeFlags.SANITIZE_KEKULIZE |
                                              Chem.SanitizeFlags.SANITIZE_SETAROMATICITY |
                                              Chem.SanitizeFlags.SANITIZE_SETCONJUGATION |
                                              Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION,
                        catchErrors=True)
    except Exception:
        pass
    return mol_final


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
    """Use the SAME mu values as build_mdx_angle_terms so RMSE measures
    what the L-BFGS is trying to minimise."""
    from delfin.fffree.grip_fragment_detect import build_mdx_angle_terms
    if metal < 0 or not donors:
        return None
    terms = build_mdx_angle_terms(mol, metal, donors, include_hydrogens=False)
    if not terms:
        return None
    devs_sq = []; devs = []
    for t in terms:
        u = P[t.a] - P[t.b]; v = P[t.c] - P[t.b]
        nu = float(np.linalg.norm(u)); nv = float(np.linalg.norm(v))
        if nu < 1e-9 or nv < 1e-9:
            continue
        cos_t = max(-1.0, min(1.0, float(np.dot(u, v) / (nu * nv))))
        import math
        th = math.degrees(math.acos(cos_t))
        devs_sq.append((th - t.mu) ** 2)
        devs.append(abs(th - t.mu))
    if not devs_sq:
        return None
    return {
        "n": len(devs),
        "rmse": float(np.sqrt(np.mean(devs_sq))),
        "max_abs": float(np.max(devs)),
    }


def polish_one(xyz_path: Path):
    block = parse_xyz_first_frame(xyz_path)
    mol = mol_from_xyz_block(block)
    if mol is None:
        return {"file": xyz_path.name, "skip": "mol parse"}
    metal, donors = find_metal_donors(mol)
    if metal < 0 or not donors:
        return {"file": xyz_path.name, "skip": "no metal/donors"}
    P0 = np.array(mol.GetConformer().GetPositions(), dtype=np.float64)
    # Baseline metric
    pre = measure_mdx_rmse(mol, P0, metal, donors)
    if pre is None:
        return {"file": xyz_path.name, "skip": "no M-D-X"}

    from delfin.fffree.grip_polish import grip_polish

    # OFF run
    os.environ.pop("DELFIN_FFFREE_GRIP_ANGLE_TO_METAL", None)
    os.environ.pop("DELFIN_FFFREE_GRIP_DMD_POLYHEDRON", None)
    try:
        P_off = grip_polish(P0, mol, metal, donors, geom="")
        post_off = measure_mdx_rmse(mol, np.asarray(P_off), metal, donors)
    except Exception as exc:
        return {"file": xyz_path.name, "error_off": str(exc)}

    # ON run
    os.environ["DELFIN_FFFREE_GRIP_ANGLE_TO_METAL"] = "1"
    try:
        P_on = grip_polish(P0, mol, metal, donors, geom="")
        post_on = measure_mdx_rmse(mol, np.asarray(P_on), metal, donors)
    except Exception as exc:
        return {"file": xyz_path.name, "error_on": str(exc)}
    finally:
        os.environ.pop("DELFIN_FFFREE_GRIP_ANGLE_TO_METAL", None)

    return {
        "file": xyz_path.name,
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--archive", required=True)
    ap.add_argument("--n", type=int, default=50)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    arc = Path(args.archive)
    files = sorted(p for p in arc.iterdir() if p.suffix == ".xyz")[:args.n]
    print(f"[smoke] {len(files)} files", file=sys.stderr)
    t0 = time.time()
    results = []
    for i, f in enumerate(files):
        try:
            r = polish_one(f)
        except Exception as exc:
            r = {"file": f.name, "fatal": str(exc)}
        results.append(r)
        if (i + 1) % 10 == 0:
            print(f"  {i+1}/{len(files)} t={time.time()-t0:.1f}s", file=sys.stderr)
    valid = [r for r in results if "rmse_off" in r and "rmse_on" in r
             and r["rmse_off"] is not None and r["rmse_on"] is not None]
    summary = {
        "n_files": len(files),
        "n_valid": len(valid),
        "mean_rmse_pre":  float(np.mean([r["rmse_pre"] for r in valid])) if valid else None,
        "mean_rmse_off":  float(np.mean([r["rmse_off"] for r in valid])) if valid else None,
        "mean_rmse_on":   float(np.mean([r["rmse_on"] for r in valid])) if valid else None,
        "mean_delta_off": float(np.mean([r["rmse_off"] - r["rmse_pre"] for r in valid])) if valid else None,
        "mean_delta_on":  float(np.mean([r["rmse_on"]  - r["rmse_pre"] for r in valid])) if valid else None,
        "mean_max_pre":  float(np.mean([r["max_pre"] for r in valid])) if valid else None,
        "mean_max_off":  float(np.mean([r["max_off"] for r in valid])) if valid else None,
        "mean_max_on":   float(np.mean([r["max_on"]  for r in valid])) if valid else None,
        "n_files_improved":     sum(1 for r in valid if r["rmse_on"] < r["rmse_off"]),
        "n_files_worsened":     sum(1 for r in valid if r["rmse_on"] > r["rmse_off"]),
        "n_files_equal":        sum(1 for r in valid if r["rmse_on"] == r["rmse_off"]),
    }
    Path(args.out).write_text(json.dumps(
        {"summary": summary, "per_file": results}, indent=2))
    print(json.dumps(summary, indent=2))
    print(f"wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
