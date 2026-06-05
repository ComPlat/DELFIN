#!/usr/bin/env python3
"""Donor-cone angle forensik on the 2792332 voll-pool.

For each XYZ structure in the archive:
  1. Identify metal + sigma-donors via covalent-cutoff contact analysis.
  2. For each donor D and each heavy non-metal/non-donor neighbour X,
     measure theta = angle(M-D-X) and compare to the CCDC-grounded mu
     loaded from the GRIP library (v4 if DELFIN_GRIP_LIB_PATH points to
     it, otherwise VSEPR fallback).
  3. Aggregate per-CN class: how many M-D-X angles deviate by more than
     1, 2, 3 sigma from the CCDC mean.
  4. Compute the rotation around the M-D axis that would minimise the
     L2 distance between the observed X positions and the substituent
     positions of the same donor's ligand mean structure (cone azimuth
     deviation).

The donor-cone DoF question this answers: how many donors are visibly
mis-oriented in their azimuthal rotation around the M-D bond axis? If
the distribution is wide and far from CCDC ideal, the cone DoF will
have measurable effect; if narrow / centred, the cone DoF will be a
diagnostic-only feature.

USAGE:
    PYTHONHASHSEED=0 DELFIN_GRIP_LIB_PATH=/path/to/grip_lib_v4.npz \
        python scripts/forensik_donor_cone_distribution.py \
            --archive .../2792332-aromatic-symmetry-VOLLPOOL \
            --n 500 --out /tmp/cone_forensik.json
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np


# Standard TM symbols (lex-fixed)
TM_SET = frozenset({
    "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
    "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
    "La","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
    "Ce","Pr","Nd","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
    "Th","U",
})

# Covalent radii (Å) -- same as polyhedra.COV subset.
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


def parse_xyz_first_frame(path: Path) -> Tuple[List[str], np.ndarray]:
    """Return (symbols, coords) of the first frame of an XYZ file."""
    lines = path.read_text().splitlines()
    n = int(lines[0].strip())
    syms: List[str] = []
    xyz: List[List[float]] = []
    for line in lines[2:2 + n]:
        parts = line.split()
        syms.append(parts[0])
        xyz.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return syms, np.asarray(xyz, dtype=np.float64)


def identify_metal_donors(
    syms: Sequence[str], P: np.ndarray,
) -> Tuple[int, List[int]]:
    """Pick metal as the first TM atom; donors as non-H non-TM atoms
    within 1.30 x (cov(M) + cov(D)) of M (slightly relaxed cutoff for
    polished structures where M-D may stretch slightly)."""
    metal = -1
    for i, s in enumerate(syms):
        if s in TM_SET:
            metal = i
            break
    if metal < 0:
        return -1, []
    rM = P[metal]
    rcM = COV.get(syms[metal], 1.5)
    donors: List[int] = []
    for j, s in enumerate(syms):
        if j == metal:
            continue
        if s == "H" or s in TM_SET:
            continue
        rcD = COV.get(s, 0.75)
        cut = 1.30 * (rcM + rcD)
        if np.linalg.norm(P[j] - rM) <= cut:
            donors.append(j)
    return metal, donors


def angle_deg(a: int, b: int, c: int, P: np.ndarray) -> Optional[float]:
    u = P[a] - P[b]
    v = P[c] - P[b]
    nu = float(np.linalg.norm(u))
    nv = float(np.linalg.norm(v))
    if nu < 1e-9 or nv < 1e-9:
        return None
    cos_t = float(np.dot(u, v) / (nu * nv))
    cos_t = max(-1.0, min(1.0, cos_t))
    return math.degrees(math.acos(cos_t))


def vsepr_fallback_mu_sigma(donor_sym: str, n_heavy: int) -> Tuple[float, float]:
    """VSEPR fallback mu/sigma based on donor element + heavy-neighbour
    count (sp3 -> 109.5, sp2 -> 120, sp -> 180)."""
    if donor_sym in ("N", "O", "S", "Se"):
        if n_heavy >= 3:
            return 109.5, 12.0   # sp3
        if n_heavy == 2:
            return 120.0, 15.0   # sp2
        return 180.0, 10.0       # sp (nitrile)
    if donor_sym == "C":
        if n_heavy >= 3:
            return 109.5, 12.0
        if n_heavy == 2:
            return 120.0, 15.0
        return 180.0, 10.0
    if donor_sym in ("P", "As"):
        return 100.0, 12.0
    if donor_sym in ("F", "Cl", "Br", "I"):
        return 109.5, 15.0
    return 109.5, 15.0


def covalent_neighbours(j: int, syms: Sequence[str], P: np.ndarray) -> List[int]:
    """Heavy non-self atoms within 1.25 x (cov(j)+cov(k)) of atom j."""
    rj = P[j]
    sj = syms[j]
    cj = COV.get(sj, 0.7)
    out: List[int] = []
    for k in range(len(syms)):
        if k == j:
            continue
        sk = syms[k]
        ck = COV.get(sk, 0.7)
        cut = 1.25 * (cj + ck)
        if np.linalg.norm(P[k] - rj) <= cut:
            out.append(k)
    return out


def cn_class(cn: int) -> str:
    if cn <= 3:
        return "CN<=3"
    if cn == 4:
        return "CN4"
    if cn == 5:
        return "CN5"
    if cn == 6:
        return "CN6"
    if cn == 7:
        return "CN7"
    return "CN>=8"


def cone_azimuth_deviation_deg(
    metal: int, donor: int, x_atoms: Sequence[int], P: np.ndarray,
) -> Optional[float]:
    """Geometric measure of "cone azimuth disorder".

    For a donor with k heavy neighbours, the ideal cone azimuth pattern
    distributes them ~uniformly (e.g. 2 -> 180 deg apart, 3 -> 120 deg
    apart) around the M-D axis. We compute the L2 deviation of the
    observed neighbour azimuths (mod 2pi) from the ideal pattern and
    return it as a degree-scale dispersion.

    Returns None when fewer than 2 X atoms (no azimuth defined for 0/1).
    """
    if len(x_atoms) < 2:
        return None
    axis = P[donor] - P[metal]
    nax = float(np.linalg.norm(axis))
    if nax < 1e-9:
        return None
    axis = axis / nax
    # Build orthonormal basis (e1, e2) perpendicular to axis.
    # Pick e1 by Gram-Schmidt against the first X atom.
    ref = P[x_atoms[0]] - P[donor]
    e1_pre = ref - np.dot(ref, axis) * axis
    n1 = float(np.linalg.norm(e1_pre))
    if n1 < 1e-9:
        # X1 is collinear with axis; fall back to any orthogonal vector.
        e1_pre = np.cross(axis, np.array([1.0, 0.0, 0.0]))
        n1 = float(np.linalg.norm(e1_pre))
        if n1 < 1e-9:
            e1_pre = np.cross(axis, np.array([0.0, 1.0, 0.0]))
            n1 = float(np.linalg.norm(e1_pre))
    e1 = e1_pre / max(n1, 1e-12)
    e2 = np.cross(axis, e1)
    azs: List[float] = []
    for x in x_atoms:
        v = P[x] - P[donor]
        v_perp = v - np.dot(v, axis) * axis
        a = math.atan2(float(np.dot(v_perp, e2)), float(np.dot(v_perp, e1)))
        azs.append(math.degrees(a))
    # Ideal: uniformly distributed by 360/k.
    k = len(azs)
    azs_sorted = sorted((a + 360.0) % 360.0 for a in azs)
    # Find best-fit rotation phi by averaging the per-atom phase
    # offsets from the uniform pattern.
    ideal = [i * (360.0 / k) for i in range(k)]
    diffs = [azs_sorted[i] - ideal[i] for i in range(k)]
    phi0 = float(np.mean(diffs))
    devs = [
        ((azs_sorted[i] - ideal[i] - phi0 + 180.0) % 360.0) - 180.0
        for i in range(k)
    ]
    rms = math.sqrt(float(np.mean([d * d for d in devs])))
    return rms


def measure_one_structure(
    syms: List[str], P: np.ndarray, lib,
    metal: int, donors: List[int],
) -> Dict:
    """Per-structure forensik.

    Returns a dict with:
      cn               : len(donors)
      n_mdx            : number of M-D-X angle observations
      mdx_zscores      : flat list of (theta-mu)/sigma per heavy X
      mdx_signed_devs  : flat list of theta-mu (deg)
      cone_rms_per_donor : per-donor cone-azimuth dispersion (deg)
    """
    out = {
        "cn": len(donors),
        "n_mdx": 0,
        "mdx_zscores": [],
        "mdx_signed_devs": [],
        "mdx_abs_devs": [],
        "cone_rms_per_donor": [],
        "donors_with_2plus_x": 0,
    }
    for d in donors:
        sd = syms[d]
        nbs = [k for k in covalent_neighbours(d, syms, P)
               if k != metal and k not in donors and syms[k] != "H"]
        if len(nbs) >= 2:
            out["donors_with_2plus_x"] += 1
            crms = cone_azimuth_deviation_deg(metal, d, nbs, P)
            if crms is not None:
                out["cone_rms_per_donor"].append(crms)
        if len(nbs) == 0:
            continue
        n_heavy = len(nbs)  # heavy substituent count drives VSEPR class
        for x in nbs:
            sx = syms[x]
            th = angle_deg(metal, d, x, P)
            if th is None:
                continue
            # Library lookup if available.
            mu, sigma = None, None
            if lib is not None:
                try:
                    metal_sym = syms[metal]
                    hyb_d = "sp3" if n_heavy >= 3 else ("sp2" if n_heavy == 2 else "sp")
                    hit = lib.lookup_angle(
                        metal_sym, sd, hyb_d, sx,
                        ring_size_min=-1, in_aromatic=False,
                        hyb1="*", hyb3="sp3",
                        min_n=5,
                    )
                    if hit is not None:
                        mu, sigma = float(hit[0]), float(hit[1])
                except Exception:
                    pass
            if mu is None or sigma is None or sigma <= 0:
                mu, sigma = vsepr_fallback_mu_sigma(sd, n_heavy)
            out["n_mdx"] += 1
            dev = float(th - mu)
            out["mdx_signed_devs"].append(dev)
            out["mdx_abs_devs"].append(abs(dev))
            out["mdx_zscores"].append(dev / sigma if sigma > 0 else 0.0)
    return out


def aggregate(rows: List[Dict]) -> Dict:
    """Aggregate per-file rows by CN class."""
    by_cn: Dict[str, Dict] = {}
    for r in rows:
        if not isinstance(r, dict) or "cn" not in r:
            continue
        key = cn_class(int(r["cn"]))
        b = by_cn.setdefault(key, {
            "n_struct": 0,
            "all_zscores": [],
            "all_abs_devs": [],
            "all_cone_rms": [],
            "n_donors_with_2plus_x": 0,
        })
        b["n_struct"] += 1
        b["all_zscores"].extend(r.get("mdx_zscores", []))
        b["all_abs_devs"].extend(r.get("mdx_abs_devs", []))
        b["all_cone_rms"].extend(r.get("cone_rms_per_donor", []))
        b["n_donors_with_2plus_x"] += int(r.get("donors_with_2plus_x", 0))
    summary: Dict[str, Dict] = {}
    for cn, b in sorted(by_cn.items()):
        zs = np.asarray(b["all_zscores"], dtype=np.float64)
        ad = np.asarray(b["all_abs_devs"], dtype=np.float64)
        cr = np.asarray(b["all_cone_rms"], dtype=np.float64)
        s = {
            "n_struct": int(b["n_struct"]),
            "n_mdx_angles": int(zs.size),
            "n_donors_with_2plus_x": int(b["n_donors_with_2plus_x"]),
        }
        if zs.size:
            s["mdx_zscore_mean"] = float(np.mean(zs))
            s["mdx_zscore_mean_abs"] = float(np.mean(np.abs(zs)))
            s["mdx_abs_dev_mean_deg"] = float(np.mean(ad))
            s["mdx_abs_dev_p95_deg"] = float(np.percentile(ad, 95))
            s["mdx_abs_dev_max_deg"] = float(np.max(ad))
            s["pct_abs_z_gt_1"] = float(np.mean(np.abs(zs) > 1.0))
            s["pct_abs_z_gt_2"] = float(np.mean(np.abs(zs) > 2.0))
            s["pct_abs_z_gt_3"] = float(np.mean(np.abs(zs) > 3.0))
        if cr.size:
            s["cone_rms_mean_deg"] = float(np.mean(cr))
            s["cone_rms_p95_deg"] = float(np.percentile(cr, 95))
            s["cone_rms_max_deg"] = float(np.max(cr))
            s["pct_cone_rms_gt_10deg"] = float(np.mean(cr > 10.0))
            s["pct_cone_rms_gt_20deg"] = float(np.mean(cr > 20.0))
        summary[cn] = s
    return summary


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--archive", required=True)
    ap.add_argument("--n", type=int, default=500)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    repo_root = Path(__file__).resolve().parent.parent
    sys.path.insert(0, str(repo_root))
    try:
        from delfin.fffree.grip_mogul_lookup import get_default_library
        lib = get_default_library()
    except Exception as exc:
        print(f"[warn] failed to load grip library: {exc!r}", file=sys.stderr)
        lib = None

    arc = Path(args.archive)
    files = sorted(p for p in arc.iterdir() if p.suffix == ".xyz")[: args.n]
    print(f"[forensik] {len(files)} files, lib={lib is not None}", file=sys.stderr)
    t0 = time.time()
    rows: List[Dict] = []
    for i, f in enumerate(files):
        try:
            syms, P = parse_xyz_first_frame(f)
            metal, donors = identify_metal_donors(syms, P)
            if metal < 0 or len(donors) < 2:
                rows.append({"file": f.name, "skip": "no TM or CN<2"})
                continue
            r = measure_one_structure(syms, P, lib, metal, donors)
            r["file"] = f.name
            r["metal"] = syms[metal]
            rows.append(r)
        except Exception as exc:
            rows.append({"file": f.name, "error": repr(exc)})
        if (i + 1) % 100 == 0:
            print(f"  {i+1}/{len(files)} t={time.time()-t0:.1f}s",
                  file=sys.stderr)
    summary = aggregate(rows)
    Path(args.out).write_text(json.dumps(
        {"summary": summary, "n_files": len(files)}, indent=2,
    ))
    print(json.dumps(summary, indent=2))
    print(f"[forensik] wrote {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
