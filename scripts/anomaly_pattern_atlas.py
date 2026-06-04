#!/usr/bin/env python3
"""Mission 4 — CCDC anomaly pattern atlas.

For a target voll-pool archive, run CCDC's GeometryAnalyser on a stratified
random sample (default 500 structures) and catalog the top-N fragment
patterns with the highest mean |z|.

For each top pattern, also dump a few example refcodes (the pool filename
typically encodes a CCDC refcode) and a suggested fix-direction:

    "build"  : feature class is wrong by construction (e.g. metal–donor
               geometry, hapto-mode, polyhedron) — fix at build time
    "polish" : feature class is sensitive to internals (bond / angle /
               torsion) — addressable by GRIP / Mogul-polish stage
    "library": fragment underrepresented in our COD library —
               expand fragment lookup table

USAGE
    PYTHONHASHSEED=0 python scripts/anomaly_pattern_atlas.py \\
        --archive /path/to/eaee0c2-new-stack-VOLLPOOL \\
        --n 500 --workers 6 \\
        --out paper_data/anomaly_pattern_atlas.csv

Heavy I/O; ~30 s/structure (CCDC analyse_molecule).
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import random
import re
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

# --------------------- CCDC + obabel paths ---------------------
CCDC_PY = Path(
    "/home/qmchem_max/CCDC/CSD_2026.1/ccdc-software/csd-python-api/run_csd_python_api"
)
OBABEL = "/home/qmchem_max/micromamba/bin/obabel"

# --------------------- Heuristic fix-direction lookup ---------------------
# Pattern: substring in CCDC fragment_label → fix direction
_FIX_DIRECTION_RULES = [
    # build-time issues — these are coordination-sphere / hapto
    (re.compile(r"\bM-(N|C|O|P|S|Cl|Br|I|F)\b"), "build",
     "Metal–donor bond — build-time hapto/coordination fix"),
    (re.compile(r"\beta[2-9]"), "build", "Hapto mode — build-time fix"),
    (re.compile(r"\bcarbonyl\b"), "polish", "Carbonyl backbone — polish"),
    # polish — internals
    (re.compile(r"\bring\b"), "polish", "Ring puckering / planarity — polish"),
    (re.compile(r"\baromatic\b"), "polish", "Aromatic flatness — polish"),
    (re.compile(r"^C-C$|^C=C$|^C#C$|^Ar-C$"), "polish",
     "C-C internal bond — polish"),
    (re.compile(r"^C-(N|O|S|P|F|Cl)$|^C=(N|O|S)$"), "polish",
     "C–X internal bond — polish"),
    (re.compile(r"^N-H$|^O-H$|^C-H$"), "polish", "X–H internal bond — polish"),
    # library fallback
    (re.compile(r"\bL0\b|\bgeneralised\b", re.IGNORECASE), "library",
     "Generalised — COD library under-represented"),
]


def classify_fix(label: str) -> tuple[str, str]:
    for pat, direction, reason in _FIX_DIRECTION_RULES:
        if pat.search(label or ""):
            return direction, reason
    return "polish", "default — internals polish"


# --------------------- XYZ → MOL2 ---------------------
def extract_first_frame(xyz_path: Path, dst: Path) -> bool:
    try:
        lines = xyz_path.read_text().splitlines()
        if len(lines) < 3:
            return False
        try:
            na = int(lines[0].strip())
        except Exception:
            return False
        if len(lines) < na + 2:
            return False
        dst.write_text("\n".join(lines[: na + 2]) + "\n")
        return True
    except Exception:
        return False


def xyz_to_mol2(xyz_path: Path, mol2_path: Path) -> bool:
    try:
        r = subprocess.run(
            [OBABEL, str(xyz_path), "-O", str(mol2_path)],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=30,
        )
        return r.returncode == 0 and mol2_path.exists() and mol2_path.stat().st_size > 0
    except Exception:
        return False


# --------------------- CCDC sub-script ---------------------
CCDC_SCRIPT = r"""
import json, sys, time
import ccdc.io as io
import ccdc.conformer as conf
mol2_path = sys.argv[1]; out_path = sys.argv[2]
try:
    mr = io.MoleculeReader(mol2_path)
    mol = next(iter(mr))
    ga = conf.GeometryAnalyser()
    t0 = time.time()
    an = ga.analyse_molecule(mol)
    elapsed = time.time() - t0
    out_features = []
    for axis, attr in (("bonds", "analysed_bonds"),
                       ("angles", "analysed_angles"),
                       ("torsions", "analysed_torsions"),
                       ("rings", "analysed_rings")):
        for f in (getattr(an, attr, []) or []):
            try:
                if not f.enough_hits:
                    continue
                z = f.z_score
                if z is None:
                    continue
                out_features.append({
                    "axis": axis,
                    "frag": str(f.fragment_label) if f.fragment_label else "",
                    "value": float(f.value) if f.value is not None else None,
                    "median": float(f.median) if f.median is not None else None,
                    "z": float(z),
                    "nhits": int(f.nhits),
                    "unusual": bool(f.unusual),
                    "generalised": bool(f.generalised) if hasattr(f,'generalised') else False,
                })
            except Exception:
                continue
    out = {"ok": True, "elapsed": elapsed, "features": out_features}
except Exception as e:
    out = {"ok": False, "error": f"{type(e).__name__}: {e}"}
open(out_path, "w").write(json.dumps(out))
"""


def run_ccdc(mol2: Path, out_json: Path, timeout: float = 180.0) -> dict:
    script = out_json.parent / "_ccdc_runner.py"
    script.write_text(CCDC_SCRIPT)
    try:
        r = subprocess.run(
            [str(CCDC_PY), str(script), str(mol2), str(out_json)],
            stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, timeout=timeout,
        )
        if r.returncode != 0:
            return {"ok": False, "error": f"exit {r.returncode}"}
        return json.loads(out_json.read_text())
    except subprocess.TimeoutExpired:
        return {"ok": False, "error": "timeout"}
    except Exception as e:
        return {"ok": False, "error": f"{type(e).__name__}: {e}"}


def _refcode_from_filename(name: str):
    stem = Path(name).stem
    m = re.findall(r"[A-Z]{6}", stem)
    return m[-1] if m else None


def process_one(args):
    xyz_path, work_dir, idx = args
    tmp = Path(work_dir) / f"w{idx}"
    tmp.mkdir(parents=True, exist_ok=True)
    frame = tmp / "frame.xyz"
    mol2 = tmp / "frame.mol2"
    ccdc_out = tmp / "ccdc.json"

    if not extract_first_frame(Path(xyz_path), frame):
        return {"file": xyz_path, "ok": False, "stage": "extract"}
    if not xyz_to_mol2(frame, mol2):
        return {"file": xyz_path, "ok": False, "stage": "obabel"}
    res = run_ccdc(mol2, ccdc_out, timeout=180.0)
    res["file"] = xyz_path
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--archive", required=True)
    ap.add_argument("--n", type=int, default=500)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--work-dir", default="/tmp/anomaly_atlas")
    ap.add_argument("--out", required=True, help="CSV path for top-N atlas")
    ap.add_argument("--top-n", type=int, default=20)
    ap.add_argument("--jsonl", default=None)
    args = ap.parse_args()

    arc = Path(args.archive)
    files = sorted(arc.rglob("*.xyz"))
    if not files:
        print(f"no XYZ in {arc}", file=sys.stderr)
        return 1
    random.seed(args.seed)
    sample = random.sample(files, min(args.n, len(files)))
    print(f"sampled {len(sample)} of {len(files)} from {arc.name}",
          file=sys.stderr)

    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    from multiprocessing import Pool
    tasks = [(str(p), str(work_dir), i) for i, p in enumerate(sample)]
    by_frag = defaultdict(lambda: {
        "zs": [], "n_unusual": 0, "n_total": 0,
        "example_files": [],
        "example_refcodes": [],
        "axis_counter": Counter(),
    })
    n_ok = 0
    t0 = time.time()
    with Pool(args.workers) as pool:
        for i, rec in enumerate(pool.imap_unordered(process_one, tasks)):
            if not rec.get("ok"):
                continue
            n_ok += 1
            fname = rec["file"]
            rc = _refcode_from_filename(fname)
            for f in rec.get("features", []):
                frag = (f.get("frag") or "<no-frag>").strip()
                axis = f.get("axis", "?")
                z = float(f.get("z", 0.0) or 0.0)
                d = by_frag[(axis, frag)]
                d["zs"].append(abs(z))
                d["n_total"] += 1
                if f.get("unusual"):
                    d["n_unusual"] += 1
                d["axis_counter"][axis] += 1
                if len(d["example_files"]) < 5:
                    d["example_files"].append(Path(fname).name)
                    if rc and rc not in d["example_refcodes"]:
                        d["example_refcodes"].append(rc)
            if (i + 1) % 20 == 0:
                print(f"  {i+1}/{len(sample)} ok={n_ok} "
                      f"elapsed={time.time()-t0:.1f}s",
                      file=sys.stderr, flush=True)

    rows = []
    for (axis, frag), d in by_frag.items():
        if d["n_total"] < 5:
            continue
        mean_abs_z = sum(d["zs"]) / len(d["zs"]) if d["zs"] else 0.0
        max_abs_z = max(d["zs"]) if d["zs"] else 0.0
        fix_dir, fix_reason = classify_fix(frag)
        rows.append({
            "axis": axis,
            "fragment": frag,
            "n_observations": d["n_total"],
            "n_unusual": d["n_unusual"],
            "frac_unusual": d["n_unusual"] / max(d["n_total"], 1),
            "mean_abs_z": round(mean_abs_z, 3),
            "max_abs_z": round(max_abs_z, 3),
            "example_refcodes": ",".join(d["example_refcodes"][:5]),
            "example_files": "|".join(d["example_files"][:5]),
            "fix_direction": fix_dir,
            "fix_reason": fix_reason,
        })
    rows.sort(key=lambda r: -r["mean_abs_z"])
    rows_top = rows[:args.top_n]

    out_csv = Path(args.out)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())
                           if rows else ["axis", "fragment"])
        w.writeheader()
        w.writerows(rows_top)
    print(f"WROTE {out_csv} ({len(rows_top)} rows)", file=sys.stderr)

    if args.jsonl:
        out_jsonl = Path(args.jsonl)
        out_jsonl.parent.mkdir(parents=True, exist_ok=True)
        out_jsonl.write_text("\n".join(json.dumps(r) for r in rows) + "\n")
        print(f"WROTE {out_jsonl} (full {len(rows)} rows)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
