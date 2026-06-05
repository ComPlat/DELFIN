#!/usr/bin/env python3
"""Smoke test: measure donor-cone DoF effect on M-D-X RMSE.

For each XYZ file we:
  1. Identify the metal + sigma-donors.
  2. Build a fresh GRIP polish with DONOR_CONE=0 (baseline) and ON.
  3. Measure M-D-X RMSE (vs library/VSEPR priors) and cone-azimuth
     dispersion BEFORE and AFTER the polish in both modes.

The configuration ALWAYS uses the angle-to-metal layer (ON) so the
comparison is "angle-to-metal alone" vs "angle-to-metal + donor-cone".

The script is parallelisable via --parallel; each worker is fully
independent (no shared state) and sets DELFIN_FFFREE_GRIP_DONOR_CONE
internally per call.

USAGE:
    PYTHONHASHSEED=0 DELFIN_GRIP_LIB_PATH=/path/to/grip_lib_v4.npz \\
        python scripts/smoke_donor_cone_rmse.py \\
            --archive .../2792332-aromatic-symmetry-VOLLPOOL \\
            --n 300 --parallel 40 --out /tmp/donor_cone_smoke.json
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
    lines = path.read_text().splitlines()
    n = int(lines[0].strip())
    syms: List[str] = []
    xyz: List[List[float]] = []
    for line in lines[2:2 + n]:
        parts = line.split()
        syms.append(parts[0])
        xyz.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return syms, np.asarray(xyz, dtype=np.float64)


def identify_metal_donors(syms, P) -> Tuple[int, List[int]]:
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


def angle_deg(a, b, c, P):
    u = P[a] - P[b]; v = P[c] - P[b]
    nu = float(np.linalg.norm(u)); nv = float(np.linalg.norm(v))
    if nu < 1e-9 or nv < 1e-9:
        return None
    cos_t = float(np.dot(u, v) / (nu * nv))
    cos_t = max(-1.0, min(1.0, cos_t))
    return math.degrees(math.acos(cos_t))


def vsepr_mu_sigma(donor_sym, n_heavy):
    if donor_sym in ("N", "O", "S", "Se"):
        if n_heavy >= 3:
            return 109.5, 12.0
        if n_heavy == 2:
            return 120.0, 15.0
        return 180.0, 10.0
    if donor_sym == "C":
        if n_heavy >= 3:
            return 109.5, 12.0
        if n_heavy == 2:
            return 120.0, 15.0
        return 180.0, 10.0
    if donor_sym in ("P", "As"):
        return 100.0, 12.0
    return 109.5, 15.0


def covalent_neighbours(j, syms, P):
    rj = P[j]; sj = syms[j]; cj = COV.get(sj, 0.7)
    out: List[int] = []
    for k in range(len(syms)):
        if k == j:
            continue
        sk = syms[k]; ck = COV.get(sk, 0.7)
        cut = 1.25 * (cj + ck)
        if np.linalg.norm(P[k] - rj) <= cut:
            out.append(k)
    return out


def measure_mdx(syms, P, metal, donors, lib):
    """Return (rmse_deg, mean_abs_dev_deg, max_abs_dev_deg, n_angles)
    averaged across all heavy M-D-X angles."""
    if metal < 0 or not donors:
        return None
    sqs: List[float] = []
    devs: List[float] = []
    n = 0
    for d in donors:
        sd = syms[d]
        nbs = [k for k in covalent_neighbours(d, syms, P)
               if k != metal and k not in donors and syms[k] != "H"]
        if not nbs:
            continue
        n_heavy = len(nbs)
        for x in nbs:
            sx = syms[x]
            th = angle_deg(metal, d, x, P)
            if th is None:
                continue
            mu, sigma = None, None
            if lib is not None:
                try:
                    hyb_d = "sp3" if n_heavy >= 3 else ("sp2" if n_heavy == 2 else "sp")
                    hit = lib.lookup_angle(
                        syms[metal], sd, hyb_d, sx,
                        ring_size_min=-1, in_aromatic=False,
                        hyb1="*", hyb3="sp3", min_n=5,
                    )
                    if hit is not None:
                        mu, sigma = float(hit[0]), float(hit[1])
                except Exception:
                    pass
            if mu is None or sigma is None or sigma <= 0:
                mu, sigma = vsepr_mu_sigma(sd, n_heavy)
            sqs.append((th - mu) ** 2)
            devs.append(abs(th - mu))
            n += 1
    if not sqs:
        return None
    return {
        "n": n,
        "rmse_deg": float(math.sqrt(np.mean(sqs))),
        "mean_abs_dev_deg": float(np.mean(devs)),
        "max_abs_dev_deg": float(np.max(devs)),
    }


def _try_polish(syms, P0, metal, donors, *, donor_cone: bool, fpath: Path):
    """Run grip_polish with DONOR_CONE on/off, returns polished P
    (or None on failure)."""
    # Set the env-flag for THIS process; child subprocesses inherit it.
    if donor_cone:
        os.environ["DELFIN_FFFREE_GRIP_DONOR_CONE"] = "1"
    else:
        os.environ.pop("DELFIN_FFFREE_GRIP_DONOR_CONE", None)
    # ALWAYS keep angle-to-metal ON for the comparison.
    os.environ["DELFIN_FFFREE_GRIP_ANGLE_TO_METAL"] = "1"
    try:
        from delfin.fffree.grip_polish import grip_polish
        from delfin.grip_io import perceive_bonds_from_xyz
    except Exception as exc:
        return None, f"import-fail: {exc!r}"
    # Use the production xyz->mol perception pipeline; this populates ring
    # info, hybridisation and aromaticity correctly so the GRIP fragment
    # detector emits real bond/angle priors (not the n_terms=0 short-circuit
    # the synthetic builder produced).
    try:
        mol = perceive_bonds_from_xyz(syms, P0, charge=0)
        if mol is None or mol.GetNumAtoms() != len(syms):
            return None, "mol-perceive-fail"
    except Exception as exc:
        return None, f"mol-perceive-fail: {exc!r}"
    try:
        P_pol = grip_polish(P0, mol, int(metal), [int(d) for d in donors])
        return np.asarray(P_pol, dtype=np.float64), None
    except Exception as exc:
        return None, f"polish-fail: {exc!r}"


def _process_one(args):
    fpath_str, lib_path = args
    fpath = Path(fpath_str)
    if lib_path:
        os.environ["DELFIN_GRIP_LIB_PATH"] = lib_path
    os.environ["PYTHONHASHSEED"] = "0"
    # Lazy import so the lib is fresh per worker.
    try:
        from delfin.fffree.grip_mogul_lookup import get_default_library
        lib = get_default_library()
    except Exception:
        lib = None
    try:
        syms, P = parse_xyz_first_frame(fpath)
    except Exception as exc:
        return {"file": fpath.name, "error": f"parse: {exc!r}"}
    metal, donors = identify_metal_donors(syms, P)
    if metal < 0 or len(donors) < 2:
        return {"file": fpath.name, "skip": "no TM or CN<2"}
    base_mdx = measure_mdx(syms, P, metal, donors, lib)
    if base_mdx is None:
        return {"file": fpath.name, "skip": "no M-D-X"}

    P_off, err_off = _try_polish(syms, P, metal, donors,
                                 donor_cone=False, fpath=fpath)
    if P_off is None:
        return {"file": fpath.name, "error": err_off}
    P_on, err_on = _try_polish(syms, P, metal, donors,
                               donor_cone=True, fpath=fpath)
    if P_on is None:
        return {"file": fpath.name, "error": err_on}
    off_mdx = measure_mdx(syms, P_off, metal, donors, lib)
    on_mdx = measure_mdx(syms, P_on, metal, donors, lib)
    return {
        "file": fpath.name,
        "n_donors": len(donors),
        "metal": syms[metal],
        "init": base_mdx,
        "off": off_mdx,
        "on": on_mdx,
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

    # Ensure delfin is on sys.path for the child processes as well.
    repo_root = Path(__file__).resolve().parent.parent
    sys.path.insert(0, str(repo_root))
    os.environ["PYTHONPATH"] = str(repo_root) + ":" + os.environ.get(
        "PYTHONPATH", "")

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
            for i, r in enumerate(ex.map(_process_one, tasks, chunksize=2)):
                results.append(r)
                if (i + 1) % 20 == 0:
                    print(f"  {i+1}/{len(tasks)} "
                          f"t={time.time()-t0:.1f}s", file=sys.stderr)
    # Aggregate.
    valid = [r for r in results if isinstance(r, dict) and "init" in r]
    def collect(field):
        return [r[field]["rmse_deg"] for r in valid if r[field] is not None]
    def collect_max(field):
        return [r[field]["max_abs_dev_deg"] for r in valid
                if r[field] is not None]
    summary = {
        "n_files": len(files),
        "n_valid": len(valid),
        "init_mean_rmse_deg": float(np.mean(collect("init"))) if valid else None,
        "off_mean_rmse_deg": float(np.mean(collect("off"))) if valid else None,
        "on_mean_rmse_deg": float(np.mean(collect("on"))) if valid else None,
        "init_mean_max_dev_deg": float(np.mean(collect_max("init"))) if valid else None,
        "off_mean_max_dev_deg": float(np.mean(collect_max("off"))) if valid else None,
        "on_mean_max_dev_deg": float(np.mean(collect_max("on"))) if valid else None,
    }
    if valid:
        deltas = [(r["off"]["rmse_deg"] - r["on"]["rmse_deg"])
                  for r in valid if r["off"] is not None and r["on"] is not None]
        if deltas:
            arr = np.asarray(deltas)
            summary["delta_rmse_off_minus_on_mean_deg"] = float(np.mean(arr))
            summary["delta_rmse_off_minus_on_p95_deg"] = float(np.percentile(arr, 95))
            summary["delta_rmse_off_minus_on_p05_deg"] = float(np.percentile(arr, 5))
            summary["n_files_improved"] = int(np.sum(arr > 0))
            summary["n_files_worsened"] = int(np.sum(arr < 0))
            summary["n_files_unchanged"] = int(np.sum(arr == 0))
    Path(args.out).write_text(json.dumps(
        {"summary": summary, "per_file": results}, indent=2,
    ))
    print(json.dumps(summary, indent=2))
    print(f"[smoke] wrote {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
