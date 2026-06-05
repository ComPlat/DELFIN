#!/usr/bin/env python
"""Pick a stratified-per-metal sample of CCDC refcodes for cleaning.

Source: the existing ccdc_tmc_index.jsonl. Constraint: only refcodes that
appear in per_smiles_ccdc_families.jsonl as matches (so the before/after
xrd-recall comparison is honest — same matches, same SMILES).

Strategy:
  - Build the set of "in-scope" refcodes (those referenced by family table)
  - Bucket them by metal_element via the TMC index
  - Round-robin sample to fill --n entries, weighted lightly by metal frequency
    (sqrt-scaled so rare metals also get representation)
  - Output: one refcode per line, deterministic

Usage:
  python pick_cleaning_sample.py --n 1000 --out /tmp/refcodes_to_clean.txt
"""
from __future__ import annotations

import argparse
import json
import math
import os
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("PYTHONHASHSEED", "0")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--families", default=
                    "/home/qmchem_max/agent_workspace/quality_framework/reference/"
                    "per_smiles_ccdc_families.jsonl")
    ap.add_argument("--tmc-index", default=
                    "/home/qmchem_max/agent_workspace/quality_framework/reference/"
                    "ccdc_tmc_index.jsonl")
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--out", required=True)
    ap.add_argument("--prefer-disorder", action="store_true",
                    help="bias selection towards entries flagged has_disorder=True "
                         "(stresses the cleaning pipeline more aggressively)")
    args = ap.parse_args()

    # 1. Collect in-scope refcodes from families
    in_scope = set()
    fam_path = Path(args.families)
    with fam_path.open() as fh:
        for line in fh:
            try:
                rec = json.loads(line)
            except Exception:
                continue
            for m in rec.get("matches", []):
                rc = m.get("refcode")
                if rc:
                    in_scope.add(rc)
    print(f"[pick] in-scope refcodes from families: {len(in_scope)}")

    # 2. Look up metal/disorder for each
    by_metal = defaultdict(list)
    n_seen = 0
    with Path(args.tmc_index).open() as fh:
        for line in fh:
            try:
                rec = json.loads(line)
            except Exception:
                continue
            rc = rec.get("refcode")
            if rc not in in_scope:
                continue
            n_seen += 1
            metal = rec.get("metal_element", "?")
            # tag for prefer-disorder boost
            key = metal
            if args.prefer_disorder and rec.get("has_disorder"):
                by_metal[key].insert(0, rc)  # disorder first
            else:
                by_metal[key].append(rc)
    print(f"[pick] resolved {n_seen} refcodes; metals={len(by_metal)}")

    # 3. Sqrt-weighted round-robin
    weights = {m: math.sqrt(len(lst)) for m, lst in by_metal.items()}
    total_w = sum(weights.values())
    quotas = {m: max(1, int(round(args.n * w / total_w))) for m, w in weights.items()}
    # Clamp to available pool size
    for m in list(quotas.keys()):
        quotas[m] = min(quotas[m], len(by_metal[m]))
    # Trim to args.n
    total = sum(quotas.values())
    if total > args.n:
        # remove from largest buckets first
        sorted_metals = sorted(quotas, key=lambda m: -quotas[m])
        excess = total - args.n
        for m in sorted_metals:
            if excess <= 0:
                break
            cut = min(excess, max(0, quotas[m] - 1))
            quotas[m] -= cut
            excess -= cut
    elif total < args.n:
        # pad from largest pools
        sorted_metals = sorted(weights, key=lambda m: -weights[m])
        deficit = args.n - total
        for m in sorted_metals:
            if deficit <= 0:
                break
            extra = min(deficit, len(by_metal[m]) - quotas[m])
            quotas[m] += extra
            deficit -= extra

    # 4. Emit refcodes
    selected = []
    for m in sorted(by_metal.keys()):
        q = quotas.get(m, 0)
        if q <= 0:
            continue
        selected.extend(by_metal[m][:q])
    print(f"[pick] selected {len(selected)} refcodes across {len(quotas)} metals")
    # Print quota breakdown
    sorted_quotas = sorted(quotas.items(), key=lambda kv: -kv[1])[:15]
    print(f"[pick] quotas top15: {sorted_quotas}")

    Path(args.out).write_text("\n".join(selected) + "\n")
    print(f"[pick] wrote {args.out}")


if __name__ == "__main__":
    main()
