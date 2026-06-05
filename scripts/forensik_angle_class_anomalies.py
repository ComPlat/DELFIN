#!/usr/bin/env python3
"""Forensik: classify CCDC angle anomalies into three classes per structure
(D-M-D, M-D-X, X-Y-Z organic), measure CCDC z-score distribution for each.

This addresses the question "where do we lose accuracy on angles?" — and
specifically the user's hypothesis that M-D-X angles (and donor-position
freedom around M) are under-served by GRIP.

The script:
  * samples N XYZ files from a pool
  * for each: convert to mol2, run CCDC GeometryAnalyser, identify the metal
    atom and its bonded donors
  * for each "analysed_angle" record classify by which atoms are M/D/X
  * report mean |z|, max |z|, fraction unusual (|z|>=2) per class

USAGE (delfin env):
    PYTHONHASHSEED=0 python3 scripts/forensik_angle_class_anomalies.py \
        --archive /home/qmchem_max/agent_workspace/quality_framework/xyz_archive/2792332-aromatic-symmetry-VOLLPOOL \
        --n 100 --workers 8 \
        --out /tmp/angle_class_forensik.json
"""
from __future__ import annotations

import argparse
import json
import os
import random
import subprocess
import sys
import time
from pathlib import Path
from collections import Counter, defaultdict

CCDC_PY = "/home/qmchem_max/CCDC/CSD_2026.1/ccdc-software/csd-python-api/run_csd_python_api"
OBABEL = "/home/qmchem_max/micromamba/bin/obabel"


# Standard TM set (lex-fixed)
TM_SET = frozenset({
    "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
    "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
    "La","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
    "Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
    "Th","U","Np","Pu",
})


def extract_first_frame(xyz_path: Path, dst: Path) -> bool:
    try:
        lines = xyz_path.read_text().splitlines()
        if len(lines) < 3:
            return False
        try:
            n = int(lines[0].strip())
        except Exception:
            return False
        if len(lines) < n + 2:
            return False
        dst.write_text("\n".join(lines[: n + 2]) + "\n")
        return True
    except Exception:
        return False


def xyz_to_mol2(xyz: Path, mol2: Path) -> bool:
    try:
        r = subprocess.run(
            [OBABEL, "-ixyz", str(xyz), "-omol2", "-O", str(mol2)],
            capture_output=True, timeout=30,
        )
        return mol2.exists() and mol2.stat().st_size > 0
    except Exception:
        return False


def ccdc_analyse(mol2: Path) -> dict:
    """Run CCDC GeometryAnalyser on mol2; return dict with bonds/angles + atom syms."""
    script = f"""
import json, sys
from ccdc.io import MoleculeReader
from ccdc.conformer import GeometryAnalyser
try:
    mol = list(MoleculeReader('{mol2}'))[0]
except Exception as e:
    print(json.dumps({{"error": str(e)}})); sys.exit(0)
ga = GeometryAnalyser()
ana = ga.analyse_molecule(mol)
out = {{"atoms": [], "angles": []}}
# atom symbols (index order matches the mol2 / our XYZ)
for atom in mol.atoms:
    out['atoms'].append(atom.atomic_symbol)
for a in (ana.analysed_angles or []):
    if a.unusual is None or a.nhits < 15:
        continue
    out['angles'].append({{
        "atoms": list(a.atom_indices),
        "z": float(a.z_score) if a.z_score is not None else None,
        "value": float(a.value),
        "mean": float(a.mean),
        "std": float(a.standard_deviation),
        "nhits": int(a.nhits),
        "unusual": bool(a.unusual),
    }})
print(json.dumps(out))
"""
    try:
        r = subprocess.run(
            [CCDC_PY, "-c", script],
            capture_output=True, timeout=180, text=True,
        )
        if r.returncode != 0:
            return {"error": r.stderr.splitlines()[-1] if r.stderr else "rc!=0"}
        line = r.stdout.strip().splitlines()[-1]
        return json.loads(line)
    except Exception as e:
        return {"error": str(e)}


def classify_angle(ang_atoms: list, atom_syms: list, metal_idx: int, donor_idxs: set) -> str:
    """Classify angle by atoms a-b-c.  Returns class label."""
    a, b, c = ang_atoms[0], ang_atoms[1], ang_atoms[2]
    # b is the central atom by Mogul convention
    is_metal = [i == metal_idx for i in (a, b, c)]
    is_donor = [i in donor_idxs for i in (a, b, c)]
    if is_metal[1]:
        # D-M-D' or D-M-X
        if is_donor[0] and is_donor[2]:
            return "D-M-D"
        return "D-M-X"   # very rare, but mathematically possible
    if is_metal[0] or is_metal[2]:
        # M is at endpoint, central is some atom
        if is_donor[1]:
            return "M-D-X"        # central = donor
        else:
            return "M-X-Y"        # central = non-donor (probably won't appear)
    if is_donor[1]:
        return "X-D-X"            # central = donor, but no metal — donor-internal
    return "X-Y-Z"                # pure organic


def detect_metal_donors(atom_syms: list, ang_records: list) -> tuple:
    """Pick metal as first atom whose symbol is in TM_SET; donors as atoms
    that appear adjacent to the metal in at least one analysed angle (i.e.
    bonded to it as judged by CCDC)."""
    metal = -1
    for i, s in enumerate(atom_syms):
        if s in TM_SET:
            metal = i
            break
    donors: set = set()
    if metal < 0:
        return metal, donors
    for rec in ang_records:
        a = rec["atoms"]
        # if metal is at position 1 (centre) both ends are donors
        if a[1] == metal:
            donors.add(a[0]); donors.add(a[2])
        # if metal is at endpoint, the centre is a donor
        elif a[0] == metal:
            donors.add(a[1])
        elif a[2] == metal:
            donors.add(a[1])
    return metal, donors


def _job(arg):
    return process_one(*arg)


def process_one(xyz_path: Path, tmpdir: Path) -> dict:
    name = xyz_path.stem
    frame = tmpdir / f"{name}.xyz"
    mol2 = tmpdir / f"{name}.mol2"
    if not extract_first_frame(xyz_path, frame):
        return {"file": name, "error": "frame extract"}
    if not xyz_to_mol2(frame, mol2):
        return {"file": name, "error": "obabel"}
    res = ccdc_analyse(mol2)
    if "error" in res:
        return {"file": name, "error": res["error"]}
    metal, donors = detect_metal_donors(res["atoms"], res["angles"])
    if metal < 0:
        return {"file": name, "skip": "no TM"}

    classes = Counter()
    by_class = defaultdict(list)   # z values per class
    for rec in res["angles"]:
        cls = classify_angle(rec["atoms"], res["atoms"], metal, donors)
        z = rec.get("z")
        if z is None:
            continue
        classes[cls] += 1
        by_class[cls].append(abs(z))
    summary = {"file": name, "metal_sym": res["atoms"][metal], "n_donors": len(donors),
               "classes": dict(classes), "by_class_meanabs_z": {},
               "by_class_maxabs_z": {}, "by_class_unusual_count": {}}
    for cls, zs in by_class.items():
        summary["by_class_meanabs_z"][cls] = sum(zs) / len(zs)
        summary["by_class_maxabs_z"][cls] = max(zs)
        summary["by_class_unusual_count"][cls] = sum(1 for z in zs if z >= 2.0)
    return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--archive", required=True)
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    arc = Path(args.archive)
    files = sorted(p for p in arc.iterdir() if p.suffix == ".xyz")
    random.seed(args.seed)
    if len(files) > args.n:
        files = random.sample(files, args.n)
    files.sort()
    print(f"[forensik] processing {len(files)} files", file=sys.stderr)

    tmproot = Path("/tmp/angle_class_forensik_work")
    tmproot.mkdir(exist_ok=True)
    results = []
    t0 = time.time()
    # serial first (CCDC subprocess is heavyweight) — workers via Pool
    from multiprocessing import Pool
    work = [(f, tmproot / f"w{i}") for i, f in enumerate(files)]
    for _, d in work:
        d.mkdir(exist_ok=True)

    if args.workers > 1:
        with Pool(args.workers) as pool:
            for i, r in enumerate(pool.imap_unordered(_job, work)):
                results.append(r)
                if (i + 1) % 10 == 0:
                    print(f"  {i+1}/{len(files)} t={time.time()-t0:.1f}s", file=sys.stderr)
    else:
        for i, a in enumerate(work):
            results.append(_job(a))
            if (i + 1) % 10 == 0:
                print(f"  {i+1}/{len(files)} t={time.time()-t0:.1f}s", file=sys.stderr)

    # aggregate per class
    agg = defaultdict(lambda: {"n_angles": 0, "n_unusual": 0, "sum_z": 0.0,
                                "max_z": 0.0, "n_files_with": 0})
    ok = 0
    for r in results:
        if "error" in r or "skip" in r:
            continue
        ok += 1
        for cls, n in r["classes"].items():
            agg[cls]["n_angles"] += n
            agg[cls]["n_unusual"] += r["by_class_unusual_count"].get(cls, 0)
            agg[cls]["sum_z"] += r["by_class_meanabs_z"].get(cls, 0.0) * n
            agg[cls]["max_z"] = max(agg[cls]["max_z"], r["by_class_maxabs_z"].get(cls, 0.0))
            agg[cls]["n_files_with"] += 1
    summary = {"n_files_total": len(files), "n_files_ok": ok, "by_class": {}}
    for cls, d in agg.items():
        n = d["n_angles"]
        summary["by_class"][cls] = {
            "n_angles": n,
            "n_unusual": d["n_unusual"],
            "frac_unusual": d["n_unusual"] / n if n else 0.0,
            "mean_abs_z": d["sum_z"] / n if n else 0.0,
            "max_abs_z": d["max_z"],
            "n_files_with": d["n_files_with"],
        }

    Path(args.out).write_text(json.dumps(
        {"summary": summary, "per_file": results}, indent=2))
    print("WROTE", args.out, file=sys.stderr)
    print(json.dumps(summary["by_class"], indent=2))


if __name__ == "__main__":
    main()
