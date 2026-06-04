#!/usr/bin/env python
"""Mission 4: Measure full-CCDC XRD-recall on existing pool archives.

Outputs:
  paper_data/xrd_recall_full_ccdc_trajectory.csv
  paper_data/xrd_recall_per_smiles_<archive>.jsonl
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("PYTHONHASHSEED", "0")

# Force worktree to take precedence over any pip-editable installed delfin
WORKTREE_REPO = "/home/qmchem_max/ComPlat/DELFIN/.claude/worktrees/agent-aca2076a88bd461f8"
if Path(WORKTREE_REPO).exists():
    # Remove any pre-existing delfin module imports, prepend worktree
    for k in list(sys.modules.keys()):
        if k == "delfin" or k.startswith("delfin."):
            del sys.modules[k]
    sys.path.insert(0, WORKTREE_REPO)

# Load the module directly via file-path to avoid editable-install shadowing
import importlib.util
_M = importlib.util.spec_from_file_location(
    "xrd_recall_metric_full",
    Path(WORKTREE_REPO) / "delfin" / "fffree" / "xrd_recall_metric_full.py",
)
XRDM = importlib.util.module_from_spec(_M)
_M.loader.exec_module(XRDM)


DEFAULT_ARCHIVE_ROOT = "/home/qmchem_max/agent_workspace/quality_framework/xyz_archive"
DEFAULT_TABLE = "/home/qmchem_max/agent_workspace/quality_framework/reference/per_smiles_ccdc_families.jsonl"
DEFAULT_INDEX = "/home/qmchem_max/agent_workspace/quality_framework/reference/ccdc_tmc_index.jsonl"
DEFAULT_MASTER = "/home/qmchem_max/agent_workspace/quality_framework/pools/smiles_master_v3_plus.txt"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--archives", nargs="+", default=[
        "b00f9a0-full7-VOLLPOOL",
        "c03a550-race-full-stack-VOLLPOOL",
        "ec7fb0d-full9-fblock-VOLLPOOL",
    ])
    ap.add_argument("--archive-root", default=DEFAULT_ARCHIVE_ROOT)
    ap.add_argument("--table-path", default=DEFAULT_TABLE)
    ap.add_argument("--index-path", default=DEFAULT_INDEX)
    ap.add_argument("--master-pool", default=DEFAULT_MASTER)
    ap.add_argument("--rmsd-A", type=float, default=0.5)
    ap.add_argument("--no-per-smiles", action="store_true",
                    help="Skip per-SMILES detailed report (3-5min/archive)")
    ap.add_argument("--paper-data-dir", default=str(
        Path(__file__).resolve().parent.parent.parent.parent /
        "ComPlat" / "DELFIN" / ".claude" / "worktrees" /
        "agent-aca2076a88bd461f8" / "paper_data"
    ))
    args = ap.parse_args()

    Path(args.paper_data_dir).mkdir(parents=True, exist_ok=True)

    print(f"[mission-4] table:  {args.table_path}", flush=True)
    print(f"[mission-4] index:  {args.index_path}", flush=True)
    print(f"[mission-4] master: {args.master_pool}", flush=True)

    table = XRDM.load_ccdc_family_table(args.table_path)
    index = XRDM.load_ccdc_tmc_index(args.index_path)
    master_map = XRDM.load_master_label_to_smiles(args.master_pool)

    print(f"[mission-4] loaded {len(table)} SMILES family records, "
          f"{len(index)} TMC index records, {len(master_map)} master-label entries",
          flush=True)

    csv_path = Path(args.paper_data_dir) / "xrd_recall_full_ccdc_trajectory.csv"
    with csv_path.open("w") as fcsv:
        w = csv.writer(fcsv)
        w.writerow(["archive", "n_xyz", "n_smiles_in_archive",
                    "n_smiles_with_family", "isomer_recall",
                    "conformer_recall", "rmsd_threshold_A", "elapsed_s"])

        for archive_label in args.archives:
            archive = Path(args.archive_root) / archive_label
            if not archive.exists():
                print(f"[mission-4] SKIP missing: {archive}", flush=True)
                continue
            t0 = time.time()
            by_smi = XRDM.group_archive_by_smiles(str(archive),
                                                  master_label_map=master_map)
            n_xyz = sum(len(v) for v in by_smi.values())
            n_with_fam = sum(1 for s in by_smi if table.get(s, {}).get("n_matches", 0) > 0)
            print(f"\n[mission-4] === {archive_label}  n_xyz={n_xyz} "
                  f"n_smiles={len(by_smi)} with-family={n_with_fam}", flush=True)

            if n_with_fam == 0:
                print(f"  -> no SMILES resolve to a CCDC family (master_map may be incomplete)",
                      flush=True)
                w.writerow([archive_label, n_xyz, len(by_smi), 0,
                            "nan", "nan", args.rmsd_A, round(time.time() - t0, 1)])
                continue

            print(f"  [stage] isomer-recall ({n_with_fam} SMILES with family)...", flush=True)
            t_iso = time.time()
            iso = XRDM.compute_isomer_recall(by_smi, table)
            print(f"    iso = {iso:.4f}  ({time.time()-t_iso:.1f}s)", flush=True)
            print(f"  [stage] conformer-recall (RMSD vs CCDC heavy-atom positions)...", flush=True)
            t_conf = time.time()
            conf = XRDM.compute_conformer_recall(by_smi, table, index,
                                                  rmsd_threshold=args.rmsd_A,
                                                  max_emitted_per_smiles=10)
            print(f"    conf = {conf:.4f}  ({time.time()-t_conf:.1f}s)", flush=True)
            print(f"  iso_recall = {iso:.3f}  conf_recall = {conf:.3f}", flush=True)

            # Per-SMILES detailed report (slow ~5 min; skip if requested)
            if args.no_per_smiles:
                print(f"  per-SMILES report skipped (--no-per-smiles)", flush=True)
            else:
                print(f"  [stage] per-SMILES detailed report...", flush=True)
                t_per = time.time()
                per = XRDM.per_smiles_recall_report(by_smi, table, index,
                                                    rmsd_threshold=args.rmsd_A)
                per_path = Path(args.paper_data_dir) / f"xrd_recall_per_smiles_{archive_label}.jsonl"
                with per_path.open("w") as fjs:
                    for r in per:
                        fjs.write(json.dumps(r) + "\n")
                print(f"  per-SMILES → {per_path.name}  ({time.time()-t_per:.1f}s)",
                      flush=True)

            w.writerow([archive_label, n_xyz, len(by_smi), n_with_fam,
                        round(iso, 4) if iso == iso else "nan",
                        round(conf, 4) if conf == conf else "nan",
                        args.rmsd_A, round(time.time() - t0, 1)])
            fcsv.flush()

    print(f"\n[mission-4] DONE → {csv_path}", flush=True)


if __name__ == "__main__":
    main()
