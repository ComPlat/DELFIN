#!/usr/bin/env python3
"""
G1 — DFT RMSD comparator (multi-frame, differing-frame targeted).

Compares method outputs (UFF baseline, DELFIN-fffree, DELFIN-heal) against xtb
GFN2 relaxed reference geometries on a stratified sample.

Each XYZ file in the pools is a multi-frame container (one frame per
isomer/conformer). To make every xtb run informative, we:

  1. Stratify file selection by CN/hapto class (deterministic).
  2. For each file, identify the set of frame indices where DELFIN-fffree
     produces a geometry that *differs* from the UFF baseline (>1e-3 Å).
  3. Pick the first K differing frames, padded with same-as-UFF frames if
     too few differ. This guarantees the comparator gets at least some
     "constructive" frames to test, while preserving the fallback-fraction
     signal for honesty.
  4. Run `xtb --opt loose --gfn 2 --etemp 300 -P 1` on the UFF frame and
     compute bond-length RMSD, bond-angle RMSD, and heavy-atom Kabsch RMSD
     between each method's frame and the relaxed xtb reference.

Determinism: filename-sorted strata, every-Nth file inside each stratum
(sqrt(class_count) weighting), then frames sorted ascending by index.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np


# -----------------------------------------------------------------------------
# IO helpers.
# -----------------------------------------------------------------------------

def read_xyz_all_frames(path: Path):
    frames = []
    with open(path) as f:
        lines = f.readlines()
    pos = 0
    while pos < len(lines):
        line = lines[pos].strip()
        if not line:
            pos += 1
            continue
        try:
            n = int(line)
        except ValueError:
            pos += 1
            continue
        if pos + 2 + n > len(lines):
            break
        elems = []; coords = []; ok = True
        for i in range(pos + 2, pos + 2 + n):
            parts = lines[i].split()
            if len(parts) < 4:
                ok = False; break
            elems.append(parts[0])
            try:
                coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
            except ValueError:
                ok = False; break
        if ok and len(coords) == n:
            frames.append((elems, np.asarray(coords, dtype=np.float64)))
        pos += n + 2
    return frames


def write_xyz_single(path: Path, elems, coords, comment: str = ""):
    with open(path, "w") as f:
        f.write(f"{len(elems)}\n{comment}\n")
        for e, (x, y, z) in zip(elems, coords):
            f.write(f"{e:2s} {x:15.8f} {y:15.8f} {z:15.8f}\n")


# -----------------------------------------------------------------------------
# Stratification.
# -----------------------------------------------------------------------------

def stratum(name: str) -> str:
    n = name.lower()
    if "hapto" in n or "sandwich" in n:
        return "hapto"
    m = re.search(r"CN(\d+)", name)
    if m:
        cn = int(m.group(1))
        if cn <= 3: return "CN<=3"
        if cn == 4: return "CN4"
        if cn == 5: return "CN5"
        if cn == 6: return "CN6"
        if cn >= 7: return "CN>=7"
    return "unspec"


# -----------------------------------------------------------------------------
# xtb runner.
# -----------------------------------------------------------------------------

def run_xtb(elems, coords, work_root: Path, xtb_bin: str, tag: str,
            gfn: int = 2, opt_level: str = "loose",
            timeout_s: int = 600) -> Optional[tuple[list[str], np.ndarray]]:
    work = Path(tempfile.mkdtemp(prefix=f"xtb_{tag}_", dir=str(work_root)))
    src = work / "in.xyz"
    write_xyz_single(src, elems, coords, comment="G1 xtb input")
    cmd = [
        xtb_bin, "in.xyz", "--opt", opt_level, "--gfn", str(gfn),
        "--etemp", "300", "-P", "1",
    ]
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    env["OMP_STACKSIZE"] = "1G"
    try:
        proc = subprocess.run(cmd, cwd=work, env=env, capture_output=True,
                              text=True, timeout=timeout_s)
    except subprocess.TimeoutExpired:
        shutil.rmtree(work, ignore_errors=True)
        return None
    opt_xyz = work / "xtbopt.xyz"
    if not opt_xyz.exists() or proc.returncode != 0:
        shutil.rmtree(work, ignore_errors=True)
        return None
    frames = read_xyz_all_frames(opt_xyz)
    shutil.rmtree(work, ignore_errors=True)
    if not frames:
        return None
    return frames[0]


# -----------------------------------------------------------------------------
# Geometry analysis.
# -----------------------------------------------------------------------------

COVALENT_R = {
    "H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02,
    "Br": 1.20, "I": 1.39,
    "Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.39, "Fe": 1.32,
    "Co": 1.26, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22,
    "Y": 1.90, "Zr": 1.75, "Nb": 1.64, "Mo": 1.54, "Tc": 1.47, "Ru": 1.46,
    "Rh": 1.42, "Pd": 1.39, "Ag": 1.45, "Cd": 1.44, "Hf": 1.75, "Ta": 1.70,
    "W": 1.62, "Re": 1.51, "Os": 1.44, "Ir": 1.41, "Pt": 1.36, "Au": 1.36,
    "Hg": 1.32,
}

def cov_r(elem: str) -> float:
    return COVALENT_R.get(elem.capitalize(), 1.5)

def bond_list(elems, coords, tol: float = 0.45):
    n = len(elems); bonds = []
    for i in range(n):
        ri = cov_r(elems[i])
        for j in range(i + 1, n):
            rj = cov_r(elems[j])
            d = float(np.linalg.norm(coords[i] - coords[j]))
            if 0.4 < d <= ri + rj + tol:
                bonds.append((i, j, d))
    return bonds

def bond_rmsd(elems_ref, coords_ref, coords_method):
    bonds = bond_list(elems_ref, coords_ref)
    if not bonds:
        return float("nan")
    deltas = []
    for i, j, d_ref in bonds:
        d_m = float(np.linalg.norm(coords_method[i] - coords_method[j]))
        deltas.append((d_m - d_ref) ** 2)
    return float(np.sqrt(np.mean(deltas)))

def angle_rmsd(elems_ref, coords_ref, coords_method, max_angles: int = 2000):
    bonds = bond_list(elems_ref, coords_ref)
    nbr = defaultdict(list)
    for i, j, _ in bonds:
        nbr[i].append(j); nbr[j].append(i)
    triples = []
    for c, nbrs in sorted(nbr.items()):
        nbrs = sorted(set(nbrs))
        for a_i in range(len(nbrs)):
            for b_i in range(a_i + 1, len(nbrs)):
                triples.append((nbrs[a_i], c, nbrs[b_i]))
        if len(triples) >= max_angles:
            break
    triples = triples[:max_angles]
    if not triples:
        return float("nan")

    def ang(coords, a, b, c):
        v1 = coords[a] - coords[b]; v2 = coords[c] - coords[b]
        n1 = np.linalg.norm(v1); n2 = np.linalg.norm(v2)
        if n1 < 1e-6 or n2 < 1e-6:
            return None
        cs = np.clip(np.dot(v1, v2) / (n1 * n2), -1.0, 1.0)
        return float(np.degrees(np.arccos(cs)))

    deltas = []
    for a, b, c in triples:
        ar = ang(coords_ref, a, b, c)
        am = ang(coords_method, a, b, c)
        if ar is None or am is None:
            continue
        deltas.append((am - ar) ** 2)
    if not deltas:
        return float("nan")
    return float(np.sqrt(np.mean(deltas)))

def kabsch_rmsd(coords_a, coords_b, mask=None):
    if mask is None:
        a = coords_a; b = coords_b
    else:
        a = coords_a[mask]; b = coords_b[mask]
    if a.shape != b.shape or a.shape[0] < 3:
        return float("nan")
    ca = a.mean(axis=0); cb = b.mean(axis=0)
    a0 = a - ca; b0 = b - cb
    H = a0.T @ b0
    try:
        U, S, Vt = np.linalg.svd(H)
    except np.linalg.LinAlgError:
        return float("nan")
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.eye(3); D[2, 2] = d
    R = Vt.T @ D @ U.T
    b_rot = b0 @ R.T
    diff = a0 - b_rot
    return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))


# -----------------------------------------------------------------------------
# Sample selection.
# -----------------------------------------------------------------------------

def select_samples(uff_dir: Path, fff_dir: Path, heal_dir: Path,
                   total: int) -> list[tuple[str, str]]:
    uff = set(os.listdir(uff_dir))
    fff = set(os.listdir(fff_dir))
    both = uff & fff
    by_s = defaultdict(list)
    for f in sorted(both):
        if not f.endswith(".xyz"):
            continue
        by_s[stratum(f)].append(f)
    strata = sorted(by_s.keys())
    counts = {s: len(by_s[s]) for s in strata}
    weights = {s: max(1.0, counts[s] ** 0.5) for s in strata}
    wsum = sum(weights.values())
    plan = {s: max(2, int(round(total * weights[s] / wsum))) for s in strata}
    diff = sum(plan.values()) - total
    if diff > 0:
        for s in sorted(plan.keys(), key=lambda x: plan[x], reverse=True):
            if diff <= 0: break
            if plan[s] > 2:
                plan[s] -= 1; diff -= 1
    elif diff < 0:
        for s in sorted(plan.keys(), key=lambda x: -plan[x]):
            if diff >= 0: break
            plan[s] += 1; diff += 1
    out = []
    for s in strata:
        files = by_s[s]
        k = min(plan[s], len(files))
        if k <= 0: continue
        step = max(1, len(files) // k)
        picks = [files[i * step] for i in range(k) if i * step < len(files)]
        for p in picks:
            out.append((p, s))
    return out


# -----------------------------------------------------------------------------
# Per-file processing.
# -----------------------------------------------------------------------------

@dataclass
class FrameResult:
    name: str
    stratum: str
    frame: int
    n_atoms: int
    fff_eq_uff: int = 0
    heal_eq_uff: int = 0
    fff_present: int = 0
    heal_present: int = 0
    uff_bond_rmsd: float = float("nan")
    fff_bond_rmsd: float = float("nan")
    heal_bond_rmsd: float = float("nan")
    uff_angle_rmsd: float = float("nan")
    fff_angle_rmsd: float = float("nan")
    heal_angle_rmsd: float = float("nan")
    uff_kabsch: float = float("nan")
    fff_kabsch: float = float("nan")
    heal_kabsch: float = float("nan")
    xtb_seconds: float = 0.0
    error: str = ""


def coords_equal(a, b, tol: float = 1e-3) -> bool:
    if a.shape != b.shape:
        return False
    return bool(np.max(np.abs(a - b)) < tol)


def pick_frame_indices(uff_frames, fff_frames, heal_frames, target_total: int,
                       prefer_differing: bool):
    """Return list of frame indices to evaluate. If prefer_differing, prioritize
    indices where fff or heal differs from UFF, padded with same-as-UFF."""
    n = len(uff_frames)
    if not prefer_differing:
        step = max(1, n // target_total)
        return [i for i in range(0, n, step)][:target_total]
    differing = []
    same = []
    for i in range(n):
        eu, cu = uff_frames[i]
        is_diff = False
        if i < len(fff_frames):
            ef, cf = fff_frames[i]
            if ef == eu and cf.shape == cu.shape:
                if not coords_equal(cu, cf):
                    is_diff = True
        if i < len(heal_frames):
            eh, ch = heal_frames[i]
            if eh == eu and ch.shape == cu.shape:
                if not coords_equal(cu, ch):
                    is_diff = True
        (differing if is_diff else same).append(i)
    selected = differing[:target_total]
    if len(selected) < target_total:
        # pad with same-as-UFF via even stride
        rem = target_total - len(selected)
        if same:
            step = max(1, len(same) // rem)
            selected.extend(same[::step][:rem])
    return sorted(set(selected))


def process_file(args):
    (name, stratum_lbl, uff_dir, fff_dir, heal_dir, xtb_bin, work_root,
     max_frames_per_file, max_atoms, min_atoms, prefer_differing) = args
    out_results: list[FrameResult] = []
    try:
        uff_frames = read_xyz_all_frames(Path(uff_dir) / name)
    except Exception as e:
        return [FrameResult(name=name, stratum=stratum_lbl, frame=-1, n_atoms=0,
                            error=f"read_uff:{e}")]
    fff_path = Path(fff_dir) / name
    heal_path = Path(heal_dir) / name
    fff_frames = read_xyz_all_frames(fff_path) if fff_path.exists() else []
    heal_frames = read_xyz_all_frames(heal_path) if heal_path.exists() else []
    if not uff_frames:
        return [FrameResult(name=name, stratum=stratum_lbl, frame=-1, n_atoms=0,
                            error="no_uff_frames")]
    selected = pick_frame_indices(uff_frames, fff_frames, heal_frames,
                                  max_frames_per_file, prefer_differing)
    work_root_p = Path(work_root)
    for fi in selected:
        elems_u, c_u = uff_frames[fi]
        if not (max_atoms >= len(elems_u) >= min_atoms):
            continue
        r = FrameResult(name=name, stratum=stratum_lbl, frame=fi, n_atoms=len(elems_u))

        fff_match = None
        heal_match = None
        if fi < len(fff_frames):
            ef, cf = fff_frames[fi]
            if ef == elems_u and cf.shape == c_u.shape:
                fff_match = cf
                r.fff_present = 1
                if coords_equal(c_u, cf):
                    r.fff_eq_uff = 1
        if fi < len(heal_frames):
            eh, ch = heal_frames[fi]
            if eh == elems_u and ch.shape == c_u.shape:
                heal_match = ch
                r.heal_present = 1
                if coords_equal(c_u, ch):
                    r.heal_eq_uff = 1

        tag = f"{name[:6]}_{fi:02d}"
        t0 = time.time()
        relaxed = run_xtb(elems_u, c_u, work_root_p, xtb_bin, tag=tag,
                          timeout_s=600)
        r.xtb_seconds = time.time() - t0
        if relaxed is None:
            r.error = "xtb_failed_or_timeout"
            out_results.append(r)
            continue
        e_r, c_r = relaxed
        if e_r != elems_u:
            r.error = "xtb_element_mismatch"
            out_results.append(r)
            continue

        heavy_mask = np.array([e != "H" for e in elems_u])
        r.uff_bond_rmsd = bond_rmsd(e_r, c_r, c_u)
        r.uff_angle_rmsd = angle_rmsd(e_r, c_r, c_u)
        r.uff_kabsch = kabsch_rmsd(c_r, c_u, heavy_mask)
        if fff_match is not None:
            r.fff_bond_rmsd = bond_rmsd(e_r, c_r, fff_match)
            r.fff_angle_rmsd = angle_rmsd(e_r, c_r, fff_match)
            r.fff_kabsch = kabsch_rmsd(c_r, fff_match, heavy_mask)
        if heal_match is not None:
            r.heal_bond_rmsd = bond_rmsd(e_r, c_r, heal_match)
            r.heal_angle_rmsd = angle_rmsd(e_r, c_r, heal_match)
            r.heal_kabsch = kabsch_rmsd(c_r, heal_match, heavy_mask)
        out_results.append(r)
    return out_results


# -----------------------------------------------------------------------------
# Main.
# -----------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--uff-dir", required=True)
    ap.add_argument("--fff-dir", required=True)
    ap.add_argument("--heal-dir", required=True)
    ap.add_argument("--total", type=int, default=50)
    ap.add_argument("--frames-per-file", type=int, default=4)
    ap.add_argument("--prefer-differing", type=int, default=1,
                    help="1 = prioritize frames where DELFIN differs from UFF")
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--xtb-bin", default="/home/qmchem_max/micromamba/bin/xtb")
    ap.add_argument("--work-root", default="/tmp/G1_xtb_validation/work")
    ap.add_argument("--out-csv", required=True)
    ap.add_argument("--out-json", required=True)
    ap.add_argument("--max-atoms", type=int, default=120)
    ap.add_argument("--min-atoms", type=int, default=6)
    args = ap.parse_args()

    work_root = Path(args.work_root); work_root.mkdir(parents=True, exist_ok=True)

    print(f"[G1] selecting stratified sample (target n={args.total}) ...",
          file=sys.stderr)
    picks = select_samples(Path(args.uff_dir), Path(args.fff_dir),
                           Path(args.heal_dir), args.total)
    print(f"[G1] {len(picks)} files selected; up to {args.frames_per_file} frames each "
          f"(prefer_differing={bool(args.prefer_differing)})",
          file=sys.stderr)

    work_args = [
        (name, s, args.uff_dir, args.fff_dir, args.heal_dir,
         args.xtb_bin, args.work_root, args.frames_per_file,
         args.max_atoms, args.min_atoms, bool(args.prefer_differing))
        for name, s in picks
    ]

    results: list[FrameResult] = []
    t_start = time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futures = {ex.submit(process_file, wa): wa[0] for wa in work_args}
        done = 0
        for fut in as_completed(futures):
            name = futures[fut]
            try:
                rs = fut.result()
            except Exception as e:
                rs = [FrameResult(name=name, stratum="?", frame=-1, n_atoms=0,
                                  error=f"worker:{e}")]
            results.extend(rs)
            done += 1
            if done % 5 == 0 or done == len(work_args):
                elapsed = time.time() - t_start
                ok = sum(1 for r in results if not r.error)
                print(f"[G1] {done}/{len(work_args)} files  frames={len(results)} ok={ok}  ({elapsed:.0f}s)",
                      file=sys.stderr)

    results.sort(key=lambda r: (r.name, r.frame))

    Path(args.out_csv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_csv, "w") as f:
        f.write("name,stratum,frame,n_atoms,fff_present,heal_present,"
                "fff_eq_uff,heal_eq_uff,xtb_seconds,error,"
                "uff_bond,fff_bond,heal_bond,"
                "uff_angle,fff_angle,heal_angle,"
                "uff_kabsch,fff_kabsch,heal_kabsch\n")
        for r in results:
            f.write(f"{r.name},{r.stratum},{r.frame},{r.n_atoms},"
                    f"{r.fff_present},{r.heal_present},"
                    f"{r.fff_eq_uff},{r.heal_eq_uff},"
                    f"{r.xtb_seconds:.2f},{r.error},"
                    f"{r.uff_bond_rmsd:.5f},{r.fff_bond_rmsd:.5f},{r.heal_bond_rmsd:.5f},"
                    f"{r.uff_angle_rmsd:.4f},{r.fff_angle_rmsd:.4f},{r.heal_angle_rmsd:.4f},"
                    f"{r.uff_kabsch:.5f},{r.fff_kabsch:.5f},{r.heal_kabsch:.5f}\n")

    by_strat = defaultdict(list)
    for r in results:
        if r.error: continue
        by_strat[r.stratum].append(r)

    def mean_finite(vals):
        v = [x for x in vals if x == x]
        return float(np.mean(v)) if v else float("nan")
    def median_finite(vals):
        v = [x for x in vals if x == x]
        return float(np.median(v)) if v else float("nan")

    summary = {
        "config": {
            "uff_dir": args.uff_dir, "fff_dir": args.fff_dir,
            "heal_dir": args.heal_dir,
            "total_files_requested": args.total,
            "total_files_processed": len({(r.name) for r in results}),
            "total_frames": len(results),
            "total_frames_ok": sum(1 for r in results if not r.error),
            "max_atoms": args.max_atoms, "min_atoms": args.min_atoms,
            "frames_per_file": args.frames_per_file,
            "prefer_differing": bool(args.prefer_differing),
            "xtb_bin": args.xtb_bin,
            "xtb_method": "GFN2-xTB --opt loose --etemp 300",
        },
        "per_stratum": {}, "overall": {}, "honesty_signals": {},
    }

    for s in sorted(by_strat.keys()):
        rs = by_strat[s]
        rs_fff = [r for r in rs if r.fff_present]
        rs_heal = [r for r in rs if r.heal_present]
        rs_fff_diff = [r for r in rs_fff if not r.fff_eq_uff]
        rs_heal_diff = [r for r in rs_heal if not r.heal_eq_uff]
        summary["per_stratum"][s] = {
            "n_frames": len(rs),
            "n_fff_present": len(rs_fff),
            "n_fff_differing": len(rs_fff_diff),
            "n_heal_present": len(rs_heal),
            "n_heal_differing": len(rs_heal_diff),
            "uff_bond_mean": mean_finite([r.uff_bond_rmsd for r in rs]),
            "fff_bond_mean_all": mean_finite([r.fff_bond_rmsd for r in rs_fff]),
            "fff_bond_mean_diff": mean_finite([r.fff_bond_rmsd for r in rs_fff_diff]),
            "heal_bond_mean_all": mean_finite([r.heal_bond_rmsd for r in rs_heal]),
            "heal_bond_mean_diff": mean_finite([r.heal_bond_rmsd for r in rs_heal_diff]),
            "uff_angle_mean": mean_finite([r.uff_angle_rmsd for r in rs]),
            "fff_angle_mean_all": mean_finite([r.fff_angle_rmsd for r in rs_fff]),
            "fff_angle_mean_diff": mean_finite([r.fff_angle_rmsd for r in rs_fff_diff]),
            "heal_angle_mean_all": mean_finite([r.heal_angle_rmsd for r in rs_heal]),
            "heal_angle_mean_diff": mean_finite([r.heal_angle_rmsd for r in rs_heal_diff]),
            "uff_kabsch_mean": mean_finite([r.uff_kabsch for r in rs]),
            "fff_kabsch_mean_all": mean_finite([r.fff_kabsch for r in rs_fff]),
            "fff_kabsch_mean_diff": mean_finite([r.fff_kabsch for r in rs_fff_diff]),
            "heal_kabsch_mean_all": mean_finite([r.heal_kabsch for r in rs_heal]),
            "heal_kabsch_mean_diff": mean_finite([r.heal_kabsch for r in rs_heal_diff]),
        }
    all_ok = [r for rs in by_strat.values() for r in rs]
    all_fff = [r for r in all_ok if r.fff_present]
    all_fff_diff = [r for r in all_fff if not r.fff_eq_uff]
    all_heal = [r for r in all_ok if r.heal_present]
    all_heal_diff = [r for r in all_heal if not r.heal_eq_uff]
    summary["overall"] = {
        "n_frames": len(all_ok),
        "n_fff_present": len(all_fff),
        "n_fff_differing": len(all_fff_diff),
        "n_heal_present": len(all_heal),
        "n_heal_differing": len(all_heal_diff),
        "uff_bond_mean": mean_finite([r.uff_bond_rmsd for r in all_ok]),
        "fff_bond_mean_all": mean_finite([r.fff_bond_rmsd for r in all_fff]),
        "fff_bond_mean_diff": mean_finite([r.fff_bond_rmsd for r in all_fff_diff]),
        "heal_bond_mean_all": mean_finite([r.heal_bond_rmsd for r in all_heal]),
        "heal_bond_mean_diff": mean_finite([r.heal_bond_rmsd for r in all_heal_diff]),
        "uff_bond_median": median_finite([r.uff_bond_rmsd for r in all_ok]),
        "fff_bond_median_all": median_finite([r.fff_bond_rmsd for r in all_fff]),
        "fff_bond_median_diff": median_finite([r.fff_bond_rmsd for r in all_fff_diff]),
        "uff_angle_mean": mean_finite([r.uff_angle_rmsd for r in all_ok]),
        "fff_angle_mean_all": mean_finite([r.fff_angle_rmsd for r in all_fff]),
        "fff_angle_mean_diff": mean_finite([r.fff_angle_rmsd for r in all_fff_diff]),
        "heal_angle_mean_all": mean_finite([r.heal_angle_rmsd for r in all_heal]),
        "heal_angle_mean_diff": mean_finite([r.heal_angle_rmsd for r in all_heal_diff]),
        "uff_kabsch_mean": mean_finite([r.uff_kabsch for r in all_ok]),
        "fff_kabsch_mean_all": mean_finite([r.fff_kabsch for r in all_fff]),
        "fff_kabsch_mean_diff": mean_finite([r.fff_kabsch for r in all_fff_diff]),
        "heal_kabsch_mean_all": mean_finite([r.heal_kabsch for r in all_heal]),
        "heal_kabsch_mean_diff": mean_finite([r.heal_kabsch for r in all_heal_diff]),
    }
    n_fff = len(all_fff); n_heal = len(all_heal)
    summary["honesty_signals"] = {
        "frac_fff_identical_to_uff": (sum(r.fff_eq_uff for r in all_fff) / n_fff) if n_fff else 0.0,
        "frac_heal_identical_to_uff": (sum(r.heal_eq_uff for r in all_heal) / n_heal) if n_heal else 0.0,
        "note": (
            "Frame selection prefers indices where DELFIN differs from UFF "
            "(prefer_differing=True), so frac_*_identical is BIASED LOW vs the "
            "full pool. Use these summary stats to compare methods on the "
            "frames where DELFIN actually exercised its constructive logic. "
            "The whole-pool fallback fraction is a separate measurement; see "
            "scripts/dft_rmsd_pool_fallback.py."
        ),
    }
    Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_json, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[G1] wrote {args.out_csv}", file=sys.stderr)
    print(f"[G1] wrote {args.out_json}", file=sys.stderr)


if __name__ == "__main__":
    main()
