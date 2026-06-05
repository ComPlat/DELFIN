#!/usr/bin/env python
"""Measure XRD-recall BEFORE and AFTER CCDC cleaning, on identical refcode set.

For honest comparison:
  - Take the same 1000 refcodes selected for cleaning.
  - For each refcode, look up its SMILES association via family table.
  - Group DELFIN emissions by SMILES (production archive read-only).
  - Compute conformer-recall against:
      A) RAW positions from ccdc_tmc_index.jsonl (UNCLEANED, has counter-ions etc.)
      B) CLEAN positions from ccdc_tmc_index_cleaned.jsonl (counter-ions stripped)
  - Report:
      isomer_recall (raw vs clean, but isomer label is CN-based — usually same)
      conformer_recall @ 0.5 / 1.0 / 2.0 Å (raw vs clean)
      per-metal-class breakdown

This script imports nothing from CCDC. It reads only the pre-built JSONLs.

Output:
  - xrd_recall_before_after.json (full report)
  - xrd_recall_before_after_per_class.csv (per-metal breakdown)

Usage:
  python xrd_recall_before_after.py \\
      --archive /home/qmchem_max/agent_workspace/quality_framework/xyz_archive/2792332-aromatic-symmetry-VOLLPOOL \\
      --refcodes /tmp/refcodes_to_clean_1000.txt \\
      --out /home/qmchem_max/agent_workspace/quality_framework/reports/xrd_recall_before_after.json
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("PYTHONHASHSEED", "0")

# Make delfin importable
sys.path.insert(0, "/home/qmchem_max/ComPlat/DELFIN")
sys.path.insert(0, "/home/qmchem_max/ComPlat/DELFIN/.claude/worktrees/agent-ac52cd3a778f239b9")

try:
    import numpy as np
except Exception:
    print("ERROR: numpy required", file=sys.stderr)
    raise

from delfin.fffree.xrd_recall_metric_full import (
    load_ccdc_family_table,
    load_master_label_to_smiles,
    group_archive_by_smiles,
)
from delfin.fffree.xrd_recall_metric import (
    classify_isomer,
    kabsch_rmsd_heavy,
    parse_xyz,
)


def load_index_subset(path: str, refcode_subset) -> dict:
    """Load TMC-index JSONL, keep only records whose refcode is in subset."""
    idx = {}
    with Path(path).open() as fh:
        for line in fh:
            try:
                r = json.loads(line)
            except Exception:
                continue
            rc = r.get("refcode")
            if rc in refcode_subset:
                idx[rc] = r
    return idx


def _classify(syms, P):
    """Best-effort isomer classification."""
    try:
        return classify_isomer(syms, P)
    except Exception:
        return "unknown"


def _rmsd(syms_a, Pa, syms_b, Pb):
    try:
        v = kabsch_rmsd_heavy(syms_a, Pa, syms_b, Pb)
        return v
    except Exception:
        return float("nan")


def compute_recall_against_index(by_smi, family_table, index, *,
                                 rmsd_thresholds=(0.5, 1.0, 2.0),
                                 max_emitted=50):
    """For each SMILES with family matches, compute recall against `index`.

    Returns: per_smiles list + global aggregates.
    """
    per_smiles = []
    # accumulators
    iso_recall_sum = 0.0
    iso_recall_n = 0
    conf_hits = {t: 0 for t in rmsd_thresholds}
    conf_total = 0  # total ccdc family members evaluated

    # Per-metal: by metal_element of family match
    per_metal_acc = defaultdict(lambda: {
        "n_family": 0, "n_smiles": 0,
        **{f"hit_{t}": 0 for t in rmsd_thresholds},
        "best_rmsd_sum": 0.0, "best_rmsd_n": 0,
    })

    for smi, xyz_paths in by_smi.items():
        fam = family_table.get(smi)
        if not fam or fam.get("n_matches", 0) == 0:
            continue

        # Pre-parse DELFIN emissions for this SMILES (cap)
        emitted_parsed = []
        delfin_isos = set()
        for xp in list(xyz_paths)[:max_emitted]:
            try:
                syms, P = parse_xyz(Path(xp).read_text())
                if syms and P is not None and len(P) > 0:
                    emitted_parsed.append((syms, P))
                    delfin_isos.add(_classify(syms, P))
            except Exception:
                continue

        if not emitted_parsed:
            continue

        # Restrict family matches to those whose refcode is in our subset index
        in_scope_matches = [m for m in fam["matches"] if m["refcode"] in index]
        if not in_scope_matches:
            continue

        # Isomer recall (CCDC-side unique isomers vs delfin-emitted, by CN match)
        ccdc_iso_tuples = set((m["isomer_label"], m.get("cn", -1))
                              for m in in_scope_matches)
        iso_matched = 0
        for ciso, ccn in ccdc_iso_tuples:
            # match if delfin emitted CN<n>- with same n, OR direct substring
            ok = False
            for d in delfin_isos:
                if d == ciso:
                    ok = True; break
                if d.startswith(f"CN{ccn}-") or d.startswith(f"CN{ccn}_"):
                    ok = True; break
                # geometry word match
                for w in ("OH", "TBP", "SPY", "TET", "SP", "PBP", "TTP",
                          "linear", "tetrahedral", "octahedral"):
                    if w in d and w in ciso:
                        ok = True; break
                if ok: break
            if ok:
                iso_matched += 1
        iso_r = iso_matched / max(len(ccdc_iso_tuples), 1)
        iso_recall_sum += iso_r
        iso_recall_n += 1

        # Conformer recall: for each family ref, best RMSD across emitted
        smi_per_metal = defaultdict(lambda: {
            "n_family": 0,
            **{f"hit_{t}": 0 for t in rmsd_thresholds},
            "best_rmsd_sum": 0.0, "best_rmsd_n": 0,
        })
        smi_hits = {t: 0 for t in rmsd_thresholds}
        smi_total = 0
        for fam_rec in in_scope_matches:
            rc = fam_rec["refcode"]
            ccdc_rec = index.get(rc)
            if not ccdc_rec:
                continue
            ref_syms = ccdc_rec["symbols"]
            ref_P = np.asarray(ccdc_rec["positions"], dtype=float)
            best = float("inf")
            for syms, P in emitted_parsed:
                r = _rmsd(syms, P, ref_syms, ref_P)
                if math.isfinite(r) and r < best:
                    best = r
                    # early-exit if already below smallest threshold? keep for stats
            smi_total += 1
            conf_total += 1
            metal = fam_rec.get("metal", "?")
            for t in rmsd_thresholds:
                if math.isfinite(best) and best < t:
                    conf_hits[t] += 1
                    smi_hits[t] += 1
                    smi_per_metal[metal][f"hit_{t}"] += 1
                    per_metal_acc[metal][f"hit_{t}"] += 1
            per_metal_acc[metal]["n_family"] += 1
            smi_per_metal[metal]["n_family"] += 1
            if math.isfinite(best):
                per_metal_acc[metal]["best_rmsd_sum"] += best
                per_metal_acc[metal]["best_rmsd_n"] += 1
                smi_per_metal[metal]["best_rmsd_sum"] += best
                smi_per_metal[metal]["best_rmsd_n"] += 1

        per_smiles.append({
            "smiles_label": fam.get("smiles_label"),
            "smiles": smi[:120],
            "n_emitted": len(xyz_paths),
            "n_family_in_scope": smi_total,
            "isomer_recall": iso_r,
            **{f"conformer_recall_{t}A": (smi_hits[t] / smi_total) if smi_total else 0.0
               for t in rmsd_thresholds},
        })

    global_metrics = {
        "n_smiles_with_family": iso_recall_n,
        "n_family_total": conf_total,
        "isomer_recall_mean": (iso_recall_sum / iso_recall_n)
        if iso_recall_n else float("nan"),
        **{f"conformer_recall_{t}A": (conf_hits[t] / conf_total)
           if conf_total else float("nan")
           for t in rmsd_thresholds},
    }
    # Finalize per-metal averages
    per_metal_out = {}
    for metal, acc in per_metal_acc.items():
        n_f = acc["n_family"]
        per_metal_out[metal] = {
            "n_family": n_f,
            **{f"conformer_recall_{t}A": (acc[f"hit_{t}"] / n_f) if n_f else 0.0
               for t in rmsd_thresholds},
            "mean_best_rmsd": (acc["best_rmsd_sum"] / acc["best_rmsd_n"])
            if acc["best_rmsd_n"] else float("nan"),
        }
    return per_smiles, global_metrics, per_metal_out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--archive", required=True,
                    help="Production DELFIN pool directory (read-only)")
    ap.add_argument("--families",
                    default="/home/qmchem_max/agent_workspace/quality_framework/"
                            "reference/per_smiles_ccdc_families.jsonl")
    ap.add_argument("--index-raw",
                    default="/home/qmchem_max/agent_workspace/quality_framework/"
                            "reference/ccdc_tmc_index.jsonl")
    ap.add_argument("--index-clean",
                    default="/home/qmchem_max/agent_workspace/quality_framework/"
                            "reference/ccdc_tmc_index_cleaned.jsonl")
    ap.add_argument("--refcodes", required=True,
                    help="newline-delimited refcode list (the 1000 we cleaned)")
    ap.add_argument("--master-pool",
                    default="/home/qmchem_max/agent_workspace/quality_framework/"
                            "pools/smiles_master_v3_plus.txt")
    ap.add_argument("--out", required=True)
    ap.add_argument("--rmsd-thresholds", default="0.5,1.0,2.0")
    args = ap.parse_args()

    rmsd_thresholds = tuple(float(x) for x in args.rmsd_thresholds.split(","))
    print(f"[before-after] RMSD thresholds: {rmsd_thresholds}")

    # 1. Refcode subset (1000)
    refcode_subset = set()
    for line in Path(args.refcodes).read_text().splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            refcode_subset.add(line)
    print(f"[before-after] refcode subset: {len(refcode_subset)}")

    # 2. Load family table and restrict to families with >=1 in-scope match
    family_table = load_ccdc_family_table(args.families)
    fam_in_scope = {}
    for smi, fam in family_table.items():
        in_scope = [m for m in fam.get("matches", [])
                    if m["refcode"] in refcode_subset]
        if in_scope:
            # shallow copy with filtered matches
            f = dict(fam)
            f["matches"] = in_scope
            f["n_matches"] = len(in_scope)
            fam_in_scope[smi] = f
    print(f"[before-after] SMILES with matches in subset: {len(fam_in_scope)}")

    # 3. Load BEFORE (raw) and AFTER (clean) indices, restricted to subset
    print("[before-after] loading raw index ...")
    index_raw = load_index_subset(args.index_raw, refcode_subset)
    print(f"[before-after] raw index entries: {len(index_raw)}")
    print("[before-after] loading clean index ...")
    index_clean = load_index_subset(args.index_clean, refcode_subset)
    print(f"[before-after] clean index entries: {len(index_clean)}")

    # 4. Group archive by SMILES (use master pool label map for resolution)
    master_map = load_master_label_to_smiles(args.master_pool)
    print(f"[before-after] master label-to-smiles entries: {len(master_map)}")
    by_smi = group_archive_by_smiles(args.archive, master_label_map=master_map)
    # Restrict to SMILES that have family matches in subset
    by_smi_in_scope = {smi: paths for smi, paths in by_smi.items()
                       if smi in fam_in_scope}
    print(f"[before-after] archive SMILES with in-scope family: {len(by_smi_in_scope)}")
    print(f"[before-after] total emitted xyz files in scope: "
          f"{sum(len(v) for v in by_smi_in_scope.values())}")

    if not by_smi_in_scope:
        print("ERROR: no in-scope SMILES — check archive path / master pool",
              file=sys.stderr)
        # Continue with empty so report still emits
        before_per_smiles, before_global, before_per_metal = [], {}, {}
        after_per_smiles, after_global, after_per_metal = [], {}, {}
    else:
        # 5. Compute BEFORE (raw)
        print("[before-after] computing BEFORE (raw uncleaned) ...")
        before_per_smiles, before_global, before_per_metal = compute_recall_against_index(
            by_smi_in_scope, fam_in_scope, index_raw,
            rmsd_thresholds=rmsd_thresholds,
        )
        print(f"[before-after] BEFORE: {before_global}")

        # 6. Compute AFTER (clean)
        print("[before-after] computing AFTER (cleaned) ...")
        after_per_smiles, after_global, after_per_metal = compute_recall_against_index(
            by_smi_in_scope, fam_in_scope, index_clean,
            rmsd_thresholds=rmsd_thresholds,
        )
        print(f"[before-after] AFTER: {after_global}")

    # 7. Compose report
    report = {
        "archive": args.archive,
        "n_refcodes_in_subset": len(refcode_subset),
        "n_raw_index_loaded": len(index_raw),
        "n_clean_index_loaded": len(index_clean),
        "n_smiles_with_family_in_scope": len(fam_in_scope),
        "n_smiles_with_emissions_in_scope": len(by_smi_in_scope) if by_smi_in_scope else 0,
        "rmsd_thresholds_A": list(rmsd_thresholds),
        "BEFORE": {
            "global": before_global,
            "per_metal": before_per_metal,
        },
        "AFTER": {
            "global": after_global,
            "per_metal": after_per_metal,
        },
        "per_smiles_top20": sorted(
            after_per_smiles,
            key=lambda x: -x.get(f"conformer_recall_{rmsd_thresholds[0]}A", 0),
        )[:20] if after_per_smiles else [],
        "delta_global": {
            f"d_conformer_recall_{t}A": (
                after_global.get(f"conformer_recall_{t}A", float("nan"))
                - before_global.get(f"conformer_recall_{t}A", float("nan"))
            ) for t in rmsd_thresholds
        } if after_global and before_global else {},
        "delta_isomer_recall": (
            after_global.get("isomer_recall_mean", float("nan"))
            - before_global.get("isomer_recall_mean", float("nan"))
        ) if after_global and before_global else float("nan"),
    }
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(report, indent=2, default=str))
    print(f"[before-after] wrote {args.out}")

    # 8. Per-metal CSV
    csv_path = Path(args.out).with_suffix("").as_posix() + "_per_metal.csv"
    with open(csv_path, "w", newline="") as fh:
        cols = ["metal", "n_family"]
        for t in rmsd_thresholds:
            cols.extend([f"BEFORE_conf_recall_{t}A", f"AFTER_conf_recall_{t}A",
                         f"delta_conf_recall_{t}A"])
        cols.extend(["BEFORE_mean_best_rmsd", "AFTER_mean_best_rmsd",
                     "delta_mean_best_rmsd"])
        w = csv.writer(fh)
        w.writerow(cols)
        metals = sorted(set(list(before_per_metal.keys())
                            + list(after_per_metal.keys())))
        for m in metals:
            b = before_per_metal.get(m, {})
            a = after_per_metal.get(m, {})
            row = [m, b.get("n_family", a.get("n_family", 0))]
            for t in rmsd_thresholds:
                br = b.get(f"conformer_recall_{t}A", float("nan"))
                ar = a.get(f"conformer_recall_{t}A", float("nan"))
                row.extend([f"{br:.4f}", f"{ar:.4f}",
                            f"{ar - br:+.4f}" if math.isfinite(ar - br) else "nan"])
            br_rmsd = b.get("mean_best_rmsd", float("nan"))
            ar_rmsd = a.get("mean_best_rmsd", float("nan"))
            row.extend([f"{br_rmsd:.3f}", f"{ar_rmsd:.3f}",
                        f"{ar_rmsd - br_rmsd:+.3f}"
                        if math.isfinite(ar_rmsd - br_rmsd) else "nan"])
            w.writerow(row)
    print(f"[before-after] wrote {csv_path}")

    # 9. Headlines
    print()
    print("=" * 64)
    print("XRD-RECALL BEFORE vs AFTER CLEANING (sample of "
          f"{len(refcode_subset)} refcodes)")
    print("=" * 64)
    if before_global and after_global:
        print(f"  isomer_recall:        "
              f"BEFORE={before_global.get('isomer_recall_mean', float('nan')):.4f}   "
              f"AFTER={after_global.get('isomer_recall_mean', float('nan')):.4f}   "
              f"delta="
              f"{after_global.get('isomer_recall_mean', float('nan'))-before_global.get('isomer_recall_mean', float('nan')):+.4f}")
        for t in rmsd_thresholds:
            b = before_global.get(f"conformer_recall_{t}A", float("nan"))
            a = after_global.get(f"conformer_recall_{t}A", float("nan"))
            print(f"  conformer_recall@{t}A:  "
                  f"BEFORE={b:.4f}   AFTER={a:.4f}   delta={a - b:+.4f}")
    print("=" * 64)


if __name__ == "__main__":
    main()
