#!/usr/bin/env python3
"""
G1 — fast pool-wide UFF-fallback fraction estimator.

For every file present in BOTH the UFF and the fffree (or heal) pool, count
frames where the DELFIN output is byte-identical to UFF (within 1e-3 Å). The
output is a single JSON with per-stratum and overall fallback fractions.

This complements `dft_rmsd_comparator.py`: that tool runs xtb on a stratified,
differing-biased sample (~200 frames), while this tool answers "in the entire
pool, how often did DELFIN actually deviate?".
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np


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


def read_xyz_all_frames(path):
    frames = []
    with open(path) as f:
        lines = f.readlines()
    pos = 0
    while pos < len(lines):
        line = lines[pos].strip()
        if not line:
            pos += 1; continue
        try:
            n = int(line)
        except ValueError:
            pos += 1; continue
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


def process_file(args):
    name, stratum_lbl, uff_dir, delfin_dir = args
    try:
        u = read_xyz_all_frames(os.path.join(uff_dir, name))
        d = read_xyz_all_frames(os.path.join(delfin_dir, name))
    except Exception:
        return (name, stratum_lbl, 0, 0, 0, 0)
    total = 0; same = 0; mismatched = 0
    for i in range(min(len(u), len(d))):
        eu, cu = u[i]; ed, cd = d[i]
        total += 1
        if eu != ed or cu.shape != cd.shape:
            mismatched += 1
            continue
        if np.max(np.abs(cu - cd)) < 1e-3:
            same += 1
    return (name, stratum_lbl, total, same, mismatched, len(u))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--uff-dir", required=True)
    ap.add_argument("--delfin-dir", required=True,
                    help="fffree or heal pool")
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--out-json", required=True)
    ap.add_argument("--label", default="fffree")
    args = ap.parse_args()

    uff = set(os.listdir(args.uff_dir))
    dlf = set(os.listdir(args.delfin_dir))
    common = sorted([f for f in (uff & dlf) if f.endswith(".xyz")])
    print(f"[G1-fallback] {len(common)} common files in {args.label} pool",
          file=sys.stderr)

    work = [(f, stratum(f), args.uff_dir, args.delfin_dir) for f in common]
    by_strat = defaultdict(lambda: {"total": 0, "same": 0, "mismatched": 0,
                                    "files": 0})
    overall = {"total": 0, "same": 0, "mismatched": 0, "files": 0}

    done = 0
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for fut in as_completed([ex.submit(process_file, w) for w in work]):
            name, s, total, same, mismatched, _ = fut.result()
            by_strat[s]["total"] += total
            by_strat[s]["same"] += same
            by_strat[s]["mismatched"] += mismatched
            by_strat[s]["files"] += 1
            overall["total"] += total
            overall["same"] += same
            overall["mismatched"] += mismatched
            overall["files"] += 1
            done += 1
            if done % 500 == 0:
                print(f"[G1-fallback] {done}/{len(common)}", file=sys.stderr)

    out = {
        "label": args.label,
        "uff_dir": args.uff_dir,
        "delfin_dir": args.delfin_dir,
        "overall": {
            **overall,
            "frac_identical": overall["same"] / overall["total"] if overall["total"] else 0.0,
            "frac_mismatched": overall["mismatched"] / overall["total"] if overall["total"] else 0.0,
        },
        "per_stratum": {
            s: {
                **d,
                "frac_identical": d["same"] / d["total"] if d["total"] else 0.0,
                "frac_mismatched": d["mismatched"] / d["total"] if d["total"] else 0.0,
            } for s, d in by_strat.items()
        },
    }
    Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_json, "w") as f:
        json.dump(out, f, indent=2)
    print(f"[G1-fallback] wrote {args.out_json}", file=sys.stderr)


if __name__ == "__main__":
    main()
