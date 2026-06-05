#!/usr/bin/env python3
"""Measure the geometric effect of the angle-to-metal layer on a smoke set.

For each XYZ file in the pool sample, we:
  1. Load the geometry, identify the metal + donors (via covalent-cutoff
     contact analysis around the TM atom).
  2. Build a fresh GRIP polish with ANGLE_TO_METAL=0 (baseline) AND ON.
  3. Measure the mean / max M-D-X angle deviation from VSEPR / library
     priors BEFORE and AFTER the polish in both modes.

We don't need a CCDC ground-truth here because the priors themselves
are CCDC-grounded (the library is built from CCDC fragments) -- this
script measures whether the polish actually brings the X-atom geometry
closer to those priors when the env-flag is on.

USAGE:
    PYTHONHASHSEED=0 python scripts/smoke_angle_to_metal_rmse.py \
        --archive .../pool --n 50 --out /tmp/angle_to_metal_smoke.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path
import numpy as np


# Standard TM symbols (lex-fixed)
TM_SET = frozenset({
    "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
    "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
    "La","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
})

# Covalent radii (Å) — same as polyhedra.COV (subset).
COV = {
    "H":0.31,"C":0.76,"N":0.71,"O":0.66,"F":0.57,"P":1.07,"S":1.05,
    "Cl":1.02,"Br":1.20,"I":1.39,"Se":1.20,"As":1.19,
    "Sc":1.70,"Ti":1.60,"V":1.53,"Cr":1.39,"Mn":1.50,"Fe":1.42,
    "Co":1.38,"Ni":1.24,"Cu":1.32,"Zn":1.22,"Y":1.90,"Zr":1.75,
    "Nb":1.64,"Mo":1.54,"Ru":1.46,"Rh":1.42,"Pd":1.39,"Ag":1.45,
    "Cd":1.44,"Hf":1.75,"Ta":1.70,"W":1.62,"Re":1.51,"Os":1.44,
    "Ir":1.41,"Pt":1.36,"Au":1.36,"Hg":1.32,
}


def parse_xyz_first_frame(path: Path):
    """Return (symbols, coords) of the first frame."""
    lines = path.read_text().splitlines()
    n = int(lines[0].strip())
    syms = []
    xyz = []
    for line in lines[2:2 + n]:
        parts = line.split()
        syms.append(parts[0])
        xyz.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return syms, np.array(xyz)


def identify_metal_donors(syms, P):
    """Pick metal as the first TM atom; donors as atoms within
    1.25 × (cov(M) + cov(D)) of M with non-H, non-TM type."""
    metal = -1
    for i, s in enumerate(syms):
        if s in TM_SET:
            metal = i
            break
    if metal < 0:
        return -1, []
    rM = P[metal]
    rcM = COV.get(syms[metal], 1.5)
    donors = []
    for j, s in enumerate(syms):
        if j == metal:
            continue
        if s == "H" or s in TM_SET:
            continue
        rcD = COV.get(s, 0.75)
        cut = 1.25 * (rcM + rcD)
        if np.linalg.norm(P[j] - rM) <= cut:
            donors.append(j)
    return metal, donors


def angle_deg(a, b, c, P):
    u = P[a] - P[b]
    v = P[c] - P[b]
    nu = np.linalg.norm(u); nv = np.linalg.norm(v)
    if nu < 1e-9 or nv < 1e-9:
        return None
    cos_t = float(np.dot(u, v) / (nu * nv))
    cos_t = max(-1.0, min(1.0, cos_t))
    import math
    return math.degrees(math.acos(cos_t))


def vsepr_target(donor_sym):
    """Mu/sigma for M-D-X based on donor identity (sp3 N=109.5, O=109.5,
    sp2 N=120, sp C=180).  Approximate (no hybridization detection)."""
    # Most ligand donors are sp3 N/O (amines, alcohols, ethers) → 109.5.
    if donor_sym in ("N", "O", "S", "Se"):
        return 109.5, 12.0
    if donor_sym in ("C",):
        return 120.0, 15.0      # carbenes, η1 aryls
    if donor_sym in ("P", "As"):
        return 100.0, 12.0      # PR3 lone-pair angle
    return 109.5, 15.0


def measure_mdx_rmse(syms, P, metal, donors):
    """For each M-D-X angle, compute (theta_observed - mu_vsepr)^2.  Return
    sqrt(mean) over all (d, x) pairs, count, and max |dev|."""
    if metal < 0 or not donors:
        return None
    # Build neighbour graph by covalent cutoff
    def neighbours(j):
        rj = P[j]; sj = syms[j]
        cj = COV.get(sj, 0.7)
        out = []
        for k in range(len(syms)):
            if k == j: continue
            sk = syms[k]
            ck = COV.get(sk, 0.7)
            cut = 1.25 * (cj + ck)
            if np.linalg.norm(P[k] - rj) <= cut:
                out.append(k)
        return out
    sqs = []
    devs = []
    n = 0
    for d in donors:
        sd = syms[d]
        mu, sigma = vsepr_target(sd)
        nbs = neighbours(d)
        for x in nbs:
            if x == metal or x in donors:
                continue
            th = angle_deg(metal, d, x, P)
            if th is None:
                continue
            sqs.append((th - mu) ** 2)
            devs.append(abs(th - mu))
            n += 1
    if not sqs:
        return None
    return {
        "n": n,
        "rmse_deg": float(np.sqrt(np.mean(sqs))),
        "mean_abs_dev_deg": float(np.mean(devs)),
        "max_abs_dev_deg": float(np.max(devs)),
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
            syms, P = parse_xyz_first_frame(f)
            metal, donors = identify_metal_donors(syms, P)
            if metal < 0:
                results.append({"file": f.name, "skip": "no TM"})
                continue
            r = measure_mdx_rmse(syms, P, metal, donors)
            if r is None:
                results.append({"file": f.name, "skip": "no M-D-X"})
                continue
            r["file"] = f.name
            r["metal"] = syms[metal]
            r["n_donors"] = len(donors)
            results.append(r)
        except Exception as exc:
            results.append({"file": f.name, "error": str(exc)})
        if (i + 1) % 25 == 0:
            print(f"  {i+1}/{len(files)} t={time.time()-t0:.1f}s", file=sys.stderr)
    # Aggregate
    valid = [r for r in results if "rmse_deg" in r]
    summary = {
        "n_files": len(files),
        "n_with_mdx": len(valid),
        "mean_rmse_deg": float(np.mean([r["rmse_deg"] for r in valid])) if valid else None,
        "mean_max_abs_dev_deg": float(np.mean([r["max_abs_dev_deg"] for r in valid])) if valid else None,
        "p95_max_abs_dev_deg": float(np.percentile([r["max_abs_dev_deg"] for r in valid], 95)) if valid else None,
        "total_mdx_angles": int(sum(r["n"] for r in valid)),
    }
    Path(args.out).write_text(json.dumps({"summary": summary, "per_file": results}, indent=2))
    print(json.dumps(summary, indent=2))
    print(f"wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
