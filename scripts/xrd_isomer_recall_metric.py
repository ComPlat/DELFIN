#!/usr/bin/env python3
"""Mission 3 — XRD-isomer + conformer recall on voll-pool archives.

For each archive listed (or auto-discovered), this script:
  1. Loads the CCDC ground-truth dump
     (default ``ccdc_ground_truth_48.json``).
  2. Walks each archive, groups XYZ files by CCDC refcode (inferred
     from filename), and applies :func:`delfin.fffree.xrd_recall_metric.score_pool`.
  3. Aggregates per-archive isomer / conformer recall and per-refcode
     records.
  4. Writes a CSV trajectory ``xrd_isomer_conformer_recall_trajectory.csv``
     for SI Table.

USAGE
    PYTHONHASHSEED=0 python scripts/xrd_isomer_recall_metric.py \\
        --archives b00f9a0-full7-VOLLPOOL c03a550-race-full-stack-VOLLPOOL \\
        --out paper_data/xrd_isomer_conformer_recall_trajectory.csv

Default behaviour (no --archives): score the canonical voll-pool list.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

DEFAULT_ARCHIVE_ROOT = Path(
    "/home/qmchem_max/agent_workspace/quality_framework/xyz_archive"
)
DEFAULT_REF = Path(
    "/home/qmchem_max/agent_workspace/quality_framework/reference/"
    "ccdc_ground_truth_48.json"
)
DEFAULT_ARCHIVES = [
    "b00f9a0-full7-VOLLPOOL",
    "c03a550-race-full-stack-VOLLPOOL",
    "ec7fb0d-full9-fblock-VOLLPOOL",
    "f8c9905-v2lib-construction-VOLLPOOL",
    "0df554d-fffree-BEST3-fullpool",
    "fb1ae9a-fffree-PURE_TRACK3-fullpool",
    "065f6f4-fffree-FULLPOOL",
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--archives", nargs="*", default=None,
                    help="Pool archive subdirs under DEFAULT_ARCHIVE_ROOT "
                    "(default: canonical 7)")
    ap.add_argument("--archive-root", default=str(DEFAULT_ARCHIVE_ROOT))
    ap.add_argument("--ref", default=str(DEFAULT_REF))
    ap.add_argument("--out-csv",
                    default="paper_data/xrd_isomer_conformer_recall_trajectory.csv")
    ap.add_argument("--out-jsonl",
                    default="paper_data/xrd_recall_per_refcode.jsonl")
    ap.add_argument("--rmsd-threshold-A", type=float, default=0.5)
    args = ap.parse_args()

    import delfin.fffree.xrd_recall_metric as xrd

    ref = xrd.load_reference(args.ref)
    if ref is None:
        print(f"ERROR: reference dump not found / unreadable: {args.ref}",
              file=sys.stderr)
        return 2
    n_records = sum(1 for r in ref.get("records", []) if r.get("ok"))
    print(f"Loaded {n_records} CCDC ground-truth records from {args.ref}",
          file=sys.stderr)

    archive_root = Path(args.archive_root)
    archives = args.archives or DEFAULT_ARCHIVES
    rows = []
    per_ref_lines = []
    for a in archives:
        p = archive_root / a
        if not p.exists():
            print(f"  SKIP (missing): {a}", file=sys.stderr)
            rows.append({
                "archive": a,
                "exists": 0,
                "isomer_recall": 0.0,
                "conformer_recall": 0.0,
                "n_refcodes_scored": 0,
                "n_refcodes_missing": 0,
                "rmsd_threshold_A": args.rmsd_threshold_A,
            })
            continue
        t0 = time.time()
        res = xrd.score_pool(p, ref=ref, rmsd_threshold=args.rmsd_threshold_A)
        if "error" in res:
            print(f"  ERROR {a}: {res['error']}", file=sys.stderr)
            continue
        dt = time.time() - t0
        rows.append({
            "archive": a,
            "exists": 1,
            "isomer_recall": res["isomer_recall_mean"],
            "conformer_recall": res["conformer_recall_mean"],
            "n_refcodes_scored": res["n_refcodes_scored"],
            "n_refcodes_missing": res["n_refcodes_missing"],
            "rmsd_threshold_A": res["rmsd_threshold_A"],
            "elapsed_s": round(dt, 1),
        })
        print(f"  {a}: iso={res['isomer_recall_mean']:.3f} "
              f"conf={res['conformer_recall_mean']:.3f} "
              f"scored={res['n_refcodes_scored']} "
              f"missing={res['n_refcodes_missing']} ({dt:.1f}s)",
              file=sys.stderr)
        for pr in res["per_refcode"]:
            pr["archive"] = a
            per_ref_lines.append(json.dumps(pr))

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else
                           ["archive", "isomer_recall", "conformer_recall"])
        w.writeheader()
        w.writerows(rows)
    print(f"WROTE {out_csv}", file=sys.stderr)

    out_jsonl = Path(args.out_jsonl)
    out_jsonl.parent.mkdir(parents=True, exist_ok=True)
    out_jsonl.write_text("\n".join(per_ref_lines) + "\n")
    print(f"WROTE {out_jsonl}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
