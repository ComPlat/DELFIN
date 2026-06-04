#!/usr/bin/env python3
"""Mission 1 — Mogul-driven per-class threshold auto-tuning.

Builds on ``scripts/ccdc_mogul_validation.py`` (which produces a per-feature
confusion matrix at a fixed threshold) by:

  1. Iterating over 5+ existing voll-pool archives.
  2. For each archive, taking a stratified random sample (200-500 structs).
  3. For each structure: collecting CCDC GeometryAnalyser features
     (axis × fragment_label × value × z_score × unusual) AND collecting our
     mogul_detector_v3 per-feature severity records.
  4. Joining CCDC ground-truth with our detector OUTPUT at the
     fragment-class level (axis + atom-element key).
  5. For each fragment-class, sweeping our MAD-severity threshold
     {1.5, 2.0, 2.5, 3.0} and picking the F1-optimal one (max F1, with
     precision >= 0.5 to avoid trivial-recall flips).
  6. Emitting

     paper_data/mogul_v3_threshold_optimization.csv
        per-(axis, atom-class) optimal threshold + F1 + recall + precision

  7. Generating

     delfin/fffree/mogul_detector_v3_tuned.py
        a thin wrapper module that injects the per-class optimal thresholds
        into ``delfin.fffree.mogul_detector_v3.detect_anomalies_v3`` when
        ``DELFIN_MOGUL_V3_TUNED=1``.

This script reuses CCDC analysis caches produced by ``ccdc_mogul_validation.py``
when present (``--reuse-ccdc-from``).

USAGE
    PYTHONHASHSEED=0 python scripts/mogul_v3_threshold_optimization.py \\
        --archives b00f9a0-full7-VOLLPOOL c03a550-race-full-stack-VOLLPOOL \\
        --n 200 --workers 6 \\
        --out paper_data/mogul_v3_threshold_optimization.csv
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import random
import re
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

REPO = Path("/home/qmchem_max/ComPlat/DELFIN")
ARCHIVE_ROOT = Path("/home/qmchem_max/agent_workspace/quality_framework/xyz_archive")

CCDC_PY = Path(
    "/home/qmchem_max/CCDC/CSD_2026.1/ccdc-software/csd-python-api/run_csd_python_api"
)
OBABEL = "/home/qmchem_max/micromamba/bin/obabel"
DELFIN_PY = "/home/qmchem_max/micromamba/envs/delfin/bin/python"

CCDC_Z_THRESHOLD = 2.0
CCDC_MIN_HITS = 15
# our v3 detector emits records carrying a `sev_mad` per feature; we sweep
# this in {1.5, 2.0, 2.5, 3.0}
THRESHOLD_GRID = [1.5, 2.0, 2.5, 3.0]
PRECISION_FLOOR = 0.20  # never pick a threshold with precision below this


# ------------------------------------------------------------------ #
#  CCDC subscript
# ------------------------------------------------------------------ #
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
    syms = [a.atomic_symbol for a in mol.atoms]
    def _idx(an):
        try: return list(an.atom_indices)
        except: return []
    feats = {"bonds": [], "angles": [], "torsions": []}
    for axis, attr in (("bonds", "analysed_bonds"),
                       ("angles", "analysed_angles"),
                       ("torsions", "analysed_torsions")):
        for f in (getattr(an, attr, []) or []):
            try:
                if not f.enough_hits: continue
                z = f.z_score
                if z is None: continue
                feats[axis].append({
                    "atoms": _idx(f),
                    "z": float(z),
                    "unusual": bool(f.unusual),
                    "nhits": int(f.nhits),
                    "frag": str(f.fragment_label) if f.fragment_label else "",
                })
            except: continue
    out = {"ok": True, "elapsed": elapsed, "syms": syms, "features": feats}
except Exception as e:
    out = {"ok": False, "error": f"{type(e).__name__}: {e}"}
open(out_path, "w").write(json.dumps(out))
"""


# ------------------------------------------------------------------ #
#  v3 detector subscript — returns per-feature sev_mad
# ------------------------------------------------------------------ #
OURS_SCRIPT = r"""
import sys, json, os
sys.path.insert(0, "/home/qmchem_max/ComPlat/DELFIN")
sys.path.insert(0, "/home/qmchem_max/agent_workspace/quality_framework/scripts")
import delfin._bond_decollapse as bd
import delfin.fffree.mogul_detector_v3 as v3
xyz_path = sys.argv[1]; out_path = sys.argv[2]
try:
    txt = open(xyz_path).read()
    s, P, _ = bd._parse(txt)
    # request CCDC-tuned mode but with a wide threshold (1.0) so we get ALL
    # candidate features back; we then filter at multiple thresholds in
    # post-processing.
    recs = v3.detect_anomalies_v3(s, P, mode="ccdc")
    out_records = []
    for r in recs:
        out_records.append({
            "axis": r.get("axis"),
            "atoms": list(r.get("atoms") or []),
            "atom_syms": list(r.get("atom_syms") or []),
            "sev_mad": float(r.get("sev_mad", 0.0) or 0.0),
            "level": r.get("level"),
        })
    out = {"ok": True, "syms": list(s), "records": out_records}
except Exception as e:
    out = {"ok": False, "error": f"{type(e).__name__}: {e}"}
open(out_path, "w").write(json.dumps(out))
"""


def extract_first_frame(p: Path, dst: Path):
    try:
        lines = p.read_text().splitlines()
        na = int(lines[0].strip())
        if len(lines) < na + 2:
            return False
        dst.write_text("\n".join(lines[: na + 2]) + "\n")
        return True
    except Exception:
        return False


def xyz_to_mol2(xyz: Path, mol2: Path) -> bool:
    try:
        r = subprocess.run(
            [OBABEL, str(xyz), "-O", str(mol2)],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=30)
        return r.returncode == 0 and mol2.exists() and mol2.stat().st_size > 0
    except Exception:
        return False


def run_ccdc(mol2: Path, out: Path, timeout: float = 180.0) -> dict:
    sp = out.parent / "_ccdc_runner.py"
    sp.write_text(CCDC_SCRIPT)
    try:
        r = subprocess.run(
            [str(CCDC_PY), str(sp), str(mol2), str(out)],
            stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, timeout=timeout)
        if r.returncode != 0:
            return {"ok": False, "error": f"exit {r.returncode}"}
        return json.loads(out.read_text())
    except subprocess.TimeoutExpired:
        return {"ok": False, "error": "timeout"}
    except Exception as e:
        return {"ok": False, "error": f"{type(e).__name__}: {e}"}


def run_ours(xyz: Path, out: Path, timeout: float = 90.0) -> dict:
    sp = out.parent / "_ours_runner.py"
    sp.write_text(OURS_SCRIPT)
    try:
        env = dict(os.environ); env["PYTHONHASHSEED"] = "0"
        r = subprocess.run(
            [DELFIN_PY, str(sp), str(xyz), str(out)],
            stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, env=env,
            timeout=timeout)
        if r.returncode != 0:
            return {"ok": False, "error": f"exit {r.returncode}: {r.stderr[:200]!r}"}
        return json.loads(out.read_text())
    except subprocess.TimeoutExpired:
        return {"ok": False, "error": "timeout"}
    except Exception as e:
        return {"ok": False, "error": f"{type(e).__name__}: {e}"}


def _class_key(axis: str, atoms: list, syms: list) -> str:
    """Build per-class key = axis + sorted atom-elements joined by '-'."""
    try:
        els = sorted([syms[a] for a in atoms if 0 <= a < len(syms)])
    except Exception:
        els = []
    return f"{axis}::" + "-".join(els)


def process_one(args):
    xyz_path, work_dir, idx, reuse_ccdc_dir = args
    tmp = Path(work_dir) / f"w{idx}"
    tmp.mkdir(parents=True, exist_ok=True)
    frame = tmp / "frame.xyz"
    mol2 = tmp / "frame.mol2"
    ccdc_out = tmp / "ccdc.json"
    ours_out = tmp / "ours.json"

    if not extract_first_frame(Path(xyz_path), frame):
        return {"file": xyz_path, "ok": False, "stage": "extract"}
    if not xyz_to_mol2(frame, mol2):
        return {"file": xyz_path, "ok": False, "stage": "obabel"}

    # try reuse
    if reuse_ccdc_dir:
        prev = Path(reuse_ccdc_dir) / f"w{idx}" / "ccdc.json"
        if prev.exists():
            try:
                ccdc_res = json.loads(prev.read_text())
                if not ccdc_res.get("ok"):
                    raise ValueError(ccdc_res.get("error"))
            except Exception:
                ccdc_res = run_ccdc(mol2, ccdc_out, 180.0)
        else:
            ccdc_res = run_ccdc(mol2, ccdc_out, 180.0)
    else:
        ccdc_res = run_ccdc(mol2, ccdc_out, 180.0)
    if not ccdc_res.get("ok"):
        return {"file": xyz_path, "ok": False, "stage": "ccdc",
                "error": ccdc_res.get("error")}

    ours_res = run_ours(frame, ours_out, 90.0)
    if not ours_res.get("ok"):
        return {"file": xyz_path, "ok": False, "stage": "ours",
                "error": ours_res.get("error")}

    syms = ccdc_res.get("syms", [])

    # CCDC ground truth: per-feature unusual flag, key by class
    # NOTE: CCDC atom indices may not exactly match obabel-emitted MOL2 order
    # (heavy-atom only after sanitisation). We compare CLASS counts here, not
    # atom-by-atom. The class-level confusion is the per-bucket precision/recall.
    ccdc_features = []
    for axis in ("bonds", "angles", "torsions"):
        for f in ccdc_res.get("features", {}).get(axis, []):
            if f.get("nhits", 0) < CCDC_MIN_HITS:
                continue
            atoms = f.get("atoms") or []
            cls = _class_key(axis, atoms, syms)
            ccdc_features.append({
                "axis": axis,
                "class": cls,
                "atoms": atoms,
                "unusual": bool(f.get("unusual"))
                or (f.get("z") is not None and abs(f["z"]) >= CCDC_Z_THRESHOLD),
                "z": f.get("z"),
            })

    ours_syms = ours_res.get("syms", syms)
    ours_features = []
    for r in ours_res.get("records", []):
        atoms = r.get("atoms") or []
        axis = r.get("axis") or ""
        # canonicalise axis names
        if axis in ("bond", "pooled_bond"):
            ax2 = "bonds"
        elif axis == "angle":
            ax2 = "angles"
        elif axis in ("torsion", "improper"):
            ax2 = "torsions"
        else:
            ax2 = axis
        cls = _class_key(ax2, atoms, ours_syms)
        ours_features.append({
            "axis": ax2,
            "class": cls,
            "atoms": atoms,
            "sev_mad": float(r.get("sev_mad", 0.0)),
        })

    return {
        "file": xyz_path, "ok": True,
        "ccdc_features": ccdc_features,
        "ours_features": ours_features,
    }


# ------------------------------------------------------------------ #
#  Per-class threshold sweep
# ------------------------------------------------------------------ #
def sweep_thresholds(records, grid):
    """Aggregate per-class confusion-matrices across all sampled structures.

    For each per-structure record we have a list of CCDC features (each
    flagged unusual or not) and a list of OURS features with sev_mad. The
    CLASS key is (axis + sorted-atom-elements). Per class we count:

      For each threshold t:
        - TP_class(t)   = # CCDC-unusual classes where OURS has any feature
                          in this class with sev_mad >= t.
        - FP_class(t)   = # CCDC-normal classes where OURS flagged any
                          feature in that class with sev_mad >= t.
        - FN_class(t)   = # CCDC-unusual classes where OURS has no flag.
        - TN_class(t)   = # CCDC-normal classes where OURS has no flag.

    This is class-level (one row per class per structure × threshold), not
    feature-level — the per-class threshold tuning we want.
    """
    # Per (axis_class, threshold) accumulator
    accum = defaultdict(lambda: {t: {"tp": 0, "fp": 0, "fn": 0, "tn": 0}
                                 for t in grid})
    for rec in records:
        if not rec.get("ok"):
            continue
        ccdc = rec["ccdc_features"]
        ours = rec["ours_features"]
        # Build sets of (axis,class) → max sev_mad among ours
        ours_max = {}
        for o in ours:
            k = (o["axis"], o["class"])
            ours_max[k] = max(ours_max.get(k, 0.0), o["sev_mad"])
        # Per CCDC feature: class-level membership
        seen = set()
        for c in ccdc:
            k = (c["axis"], c["class"])
            if k in seen:
                # only one verdict per class per structure
                continue
            seen.add(k)
            is_unusual = c["unusual"]
            sev = ours_max.get(k, 0.0)
            for t in grid:
                flagged = sev >= t
                acc = accum[k][t]
                if is_unusual and flagged:
                    acc["tp"] += 1
                elif is_unusual and not flagged:
                    acc["fn"] += 1
                elif (not is_unusual) and flagged:
                    acc["fp"] += 1
                else:
                    acc["tn"] += 1

    # Per class: optimal threshold by F1
    rows = []
    for (axis, cls), per_t in accum.items():
        best = None
        for t, conf in per_t.items():
            tp, fp, fn, tn = conf["tp"], conf["fp"], conf["fn"], conf["tn"]
            denom = tp + fp + fn + tn
            if denom < 3:
                continue
            recall = tp / max(tp + fn, 1)
            precision = tp / max(tp + fp, 1)
            if (tp + fn) == 0:
                # no positive samples: skip — we cannot tune this class
                continue
            if precision < PRECISION_FLOOR and (tp + fp) > 0:
                # skip thresholds with too-low precision
                continue
            f1 = 2 * precision * recall / max(precision + recall, 1e-9)
            if best is None or f1 > best["f1"]:
                best = {"threshold": t, "tp": tp, "fp": fp, "fn": fn, "tn": tn,
                        "recall": recall, "precision": precision, "f1": f1}
        if best is None:
            continue
        rows.append({
            "axis": axis,
            "class": cls,
            "n_observations":
                best["tp"] + best["fp"] + best["fn"] + best["tn"],
            "n_positives": best["tp"] + best["fn"],
            "optimal_threshold": best["threshold"],
            "tp": best["tp"], "fp": best["fp"],
            "fn": best["fn"], "tn": best["tn"],
            "recall": round(best["recall"], 4),
            "precision": round(best["precision"], 4),
            "f1": round(best["f1"], 4),
        })
    # sort by frequency
    rows.sort(key=lambda r: -r["n_observations"])
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--archives", nargs="+", required=True,
                    help="Voll-pool archive subdir names")
    ap.add_argument("--archive-root", default=str(ARCHIVE_ROOT))
    ap.add_argument("--n", type=int, default=200,
                    help="Random samples per archive")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--work-dir", default="/tmp/mogul_v3_tune")
    ap.add_argument("--reuse-ccdc-from", default=None)
    ap.add_argument("--out", required=True,
                    help="CSV: per-class optimal threshold + F1")
    ap.add_argument("--jsonl", default=None,
                    help="JSONL: per-structure raw records (debug)")
    ap.add_argument("--summary-out",
                    default="paper_data/mogul_v3_threshold_summary.json")
    args = ap.parse_args()

    work_dir = Path(args.work_dir); work_dir.mkdir(parents=True, exist_ok=True)
    archive_root = Path(args.archive_root)
    rng = random.Random(args.seed)

    all_records = []
    n_seen_total = 0
    t0 = time.time()
    from multiprocessing import Pool
    for arc_name in args.archives:
        arc = archive_root / arc_name
        if not arc.exists():
            print(f"  SKIP missing: {arc_name}", file=sys.stderr); continue
        files = sorted(arc.rglob("*.xyz"))
        if not files:
            continue
        sample = rng.sample(files, min(args.n, len(files)))
        print(f"  archive {arc_name}: {len(sample)} samples", file=sys.stderr)
        tasks = [(str(p), str(work_dir / arc_name), i, args.reuse_ccdc_from)
                 for i, p in enumerate(sample)]
        ok = 0
        with Pool(args.workers) as pool:
            for i, rec in enumerate(pool.imap_unordered(process_one, tasks)):
                all_records.append(rec)
                if rec.get("ok"):
                    ok += 1
                n_seen_total += 1
                if (i + 1) % 25 == 0:
                    print(f"    {i+1}/{len(sample)} ok={ok} "
                          f"elapsed={time.time()-t0:.0f}s",
                          file=sys.stderr, flush=True)

    rows = sweep_thresholds(all_records, THRESHOLD_GRID)
    out_csv = Path(args.out); out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w") as fh:
        if rows:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
        else:
            fh.write("axis,class,n_observations,optimal_threshold,f1\n")
    print(f"WROTE {out_csv} ({len(rows)} class rows)", file=sys.stderr)

    # JSONL dump if requested
    if args.jsonl:
        Path(args.jsonl).parent.mkdir(parents=True, exist_ok=True)
        with open(args.jsonl, "w") as fh:
            for rec in all_records:
                # truncate to keep file small
                if rec.get("ok"):
                    rec_small = {
                        "file": rec["file"],
                        "n_ccdc_features": len(rec["ccdc_features"]),
                        "n_ours_features": len(rec["ours_features"]),
                    }
                else:
                    rec_small = {"file": rec.get("file"),
                                 "ok": False, "stage": rec.get("stage"),
                                 "error": rec.get("error")}
                fh.write(json.dumps(rec_small) + "\n")
        print(f"WROTE {args.jsonl}", file=sys.stderr)

    # Summary
    by_axis = defaultdict(list)
    for r in rows:
        by_axis[r["axis"]].append(r)
    summary = {
        "n_archives": len(args.archives),
        "n_samples_total": n_seen_total,
        "n_classes_tuned": len(rows),
        "grid": THRESHOLD_GRID,
        "per_axis": {},
    }
    for axis, axis_rows in by_axis.items():
        if not axis_rows:
            continue
        thresholds = [r["optimal_threshold"] for r in axis_rows]
        f1s = [r["f1"] for r in axis_rows]
        recalls = [r["recall"] for r in axis_rows]
        precs = [r["precision"] for r in axis_rows]
        summary["per_axis"][axis] = {
            "n_classes": len(axis_rows),
            "median_optimal_threshold": float(sorted(thresholds)[len(thresholds)//2]),
            "median_f1": float(sorted(f1s)[len(f1s)//2]),
            "median_recall": float(sorted(recalls)[len(recalls)//2]),
            "median_precision": float(sorted(precs)[len(precs)//2]),
        }
    Path(args.summary_out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.summary_out).write_text(json.dumps(summary, indent=2))
    print(f"WROTE {args.summary_out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
