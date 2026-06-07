#!/usr/bin/env python3
"""extract_worst_files.py — top-N WORST files per metric axis.

Existing detector battery emits AGGREGATE percentages.  This module computes
per-file values for the diagnostic axes we care about, then ranks them.

Per-file axes (deterministic, locally computed):
  - min_md            minimum metal-donor distance (severe if < 1.5)
  - max_md            maximum metal-donor distance
  - cn                coordination number
  - min_donor_donor   shortest donor-donor pair distance through metal
  - n_h_clash         number of H-H pairs < 1.4 Å (non-bonded)
  - n_ch_short        number of C-H bonds < 0.9 Å
  - n_atom_overlap    number of heavy-atom pairs < 0.8 Å
  - bug_class_count   total bug classes that fire for this file
  - donor_drift_count number of donors in the 2.4-3.0 Å bracket
  - polyhedron_mismatch_score  1 if polyhedron classifier fires, 0 otherwise

For each top-N file the record carries:
  file_path, refcode (heuristic), metal, cn, polyhedron_class, smiles
  (from manifest pool), dispatch path (from forensik TSV), axis_value,
  contributing_atoms (small list for the top axis).

Author: hmaximilian <hmaximilian496@gmail.com>
"""
from __future__ import annotations
import argparse
import itertools
import json
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from diagnostic_common import (  # noqa: E402
    list_archive_xyz, summarise_structure, load_manifest_pool,
    load_dispatch_forensik, refcode_from_filename, label_from_filename,
    HALOGENS, dist,
)
from classify_bug_pattern import classify_atoms  # noqa: E402


# ---------------------------------------------------------------------------
# Per-file axis computation
# ---------------------------------------------------------------------------
def _per_file_axes(rec: Dict) -> Dict:
    """Compute every per-file axis used by extract_worst_files.

    rec is the output of summarise_structure().  Returns a dict mapping
    axis-name → numeric value (NaN/None for unavailable).
    """
    axes: Dict = {
        "min_md": None,
        "max_md": None,
        "mean_md": None,
        "cn": rec.get("cn", 0),
        "min_donor_donor": None,
        "n_h_clash": 0,
        "n_ch_short": 0,
        "n_atom_overlap": 0,
        "donor_drift_count": 0,
        "polyhedron_mismatch_score": 0,
        "bug_class_count": 0,
    }
    if rec.get("error"):
        return axes
    md = rec.get("donor_distances") or []
    if md:
        axes["min_md"] = min(md)
        axes["max_md"] = max(md)
        axes["mean_md"] = sum(md) / len(md)
    # donor-donor min
    dp = rec.get("donor_positions") or []
    if len(dp) >= 2:
        axes["min_donor_donor"] = min(
            dist(a, b) for a, b in itertools.combinations(dp, 2)
        )
    # H-H clashes + C-H short + heavy overlap on demand (re-walk atoms)
    atoms = rec.get("atoms") or []
    syms = [a[0] for a in atoms]
    pos = [a[1] for a in atoms]
    # bonds: re-use the bond rule from common (lazy — only if we need ch_short)
    from diagnostic_common import covalent_bonds  # local import for speed
    bonds = covalent_bonds(atoms)
    bondset = {(i, j) for i, j in bonds}
    # ch_short
    ch_short = 0
    for i, j in bonds:
        a, b = syms[i], syms[j]
        if (a == "C" and b == "H") or (a == "H" and b == "C"):
            if dist(pos[i], pos[j]) < 0.9:
                ch_short += 1
    axes["n_ch_short"] = ch_short
    # h_h_clash
    h_idx = [i for i, s in enumerate(syms) if s == "H"]
    n_hh = 0
    for a, b in itertools.combinations(h_idx, 2):
        if (a, b) in bondset or (b, a) in bondset:
            continue
        if dist(pos[a], pos[b]) < 1.4:
            n_hh += 1
    axes["n_h_clash"] = n_hh
    # atom overlap
    heavy = [i for i, s in enumerate(syms) if s != "H"]
    n_ov = 0
    for a, b in itertools.combinations(heavy, 2):
        if dist(pos[a], pos[b]) < 0.8:
            n_ov += 1
    axes["n_atom_overlap"] = n_ov
    # donor drift bracket
    axes["donor_drift_count"] = sum(1 for d in md if 2.4 < d < 3.0)
    # bug-class fire counts
    hits = classify_atoms(atoms)
    axes["bug_class_count"] = len(hits)
    if "polyhedron_mismatch_d8_t4" in hits or "polyhedron_mismatch_d10_sp4" in hits:
        axes["polyhedron_mismatch_score"] = 1
    return axes


# Direction "worse" = which extreme is bad (per-axis sort direction).
# +1 = higher is worse (sort DESCENDING).  -1 = lower is worse (sort ASCENDING).
AXIS_DIRECTION = {
    "min_md": -1,                  # lower is worse
    "max_md": +1,                  # higher is worse
    "mean_md": +1,                 # only meaningful if we say "drift up bad"
    "cn": 0,                       # neither direction meaningful as 'worst'
    "min_donor_donor": -1,         # lower = collapse
    "n_h_clash": +1,
    "n_ch_short": +1,
    "n_atom_overlap": +1,
    "donor_drift_count": +1,
    "polyhedron_mismatch_score": +1,
    "bug_class_count": +1,
}


# ---------------------------------------------------------------------------
# Top-N selection
# ---------------------------------------------------------------------------
def extract_worst_files(
    archive: str,
    axis: str,
    top_n: int = 50,
    *,
    master_pool: Optional[str] = None,
    dispatch_tsv: Optional[str] = None,
    sample: int = 0,
    seed: int = 0,
) -> List[Dict]:
    """Top-N worst files in 'archive' along 'axis'.

    axis must be a key in AXIS_DIRECTION.  Returns a list of dicts with
    file_path, smiles, axis_value, refcode, metal, cn, polyhedron_class,
    dispatch_path, contributing_atoms.
    """
    if axis not in AXIS_DIRECTION:
        raise SystemExit(
            f"unknown axis '{axis}'.  known: {sorted(AXIS_DIRECTION)}"
        )
    direction = AXIS_DIRECTION[axis]
    if direction == 0:
        raise SystemExit(
            f"axis '{axis}' is not directional; cannot rank 'worst'."
        )

    paths = list_archive_xyz(archive)
    if sample and sample < len(paths):
        import random
        rng = random.Random(seed)
        paths = sorted(rng.sample(paths, sample))

    pool = load_manifest_pool(master_pool) if master_pool else {}
    dispatch = load_dispatch_forensik(dispatch_tsv) if dispatch_tsv else {}

    rows: List[Dict] = []
    for p in paths:
        rec = summarise_structure(p)
        axes = _per_file_axes(rec)
        v = axes.get(axis)
        if v is None:
            continue
        label = label_from_filename(p)
        smi = pool.get(label, "")
        row = {
            "file_path": p,
            "label": label,
            "refcode": refcode_from_filename(p),
            "smiles": smi,
            "metal": rec.get("metal_elem"),
            "cn": axes.get("cn"),
            "polyhedron_class": rec.get("polyhedron_class"),
            "axis": axis,
            "axis_value": v,
            "dispatch_path": dispatch.get(smi, "") if smi else "",
            "all_axes": axes,
        }
        # Per-axis "contributing atom indices" hint for downstream forensik
        row["contributing_atoms"] = _contrib_atoms(rec, axis, axes)
        rows.append(row)
    # sort: descending for +1, ascending for -1
    rows.sort(key=lambda r: r["axis_value"], reverse=(direction > 0))
    return rows[:top_n]


def _contrib_atoms(rec: Dict, axis: str, axes: Dict) -> List[int]:
    """Per-axis contributing atom indices (small heuristic)."""
    if rec.get("error"):
        return []
    if axis == "min_md":
        # the donor with the smallest distance
        md = rec.get("donor_distances") or []
        if not md:
            return []
        j_min = min(range(len(md)), key=lambda i: md[i])
        return [rec.get("metal_idx"), rec["donor_indices"][j_min]]
    if axis == "max_md":
        md = rec.get("donor_distances") or []
        if not md:
            return []
        j_max = max(range(len(md)), key=lambda i: md[i])
        return [rec.get("metal_idx"), rec["donor_indices"][j_max]]
    if axis == "min_donor_donor":
        dp = rec.get("donor_positions") or []
        di = rec.get("donor_indices") or []
        best = (1e9, None, None)
        for a, b in itertools.combinations(range(len(dp)), 2):
            d = dist(dp[a], dp[b])
            if d < best[0]:
                best = (d, di[a], di[b])
        return [x for x in (best[1], best[2]) if x is not None]
    # default: no-op
    return []


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("archive", help="path to xyz_archive/<name>/")
    ap.add_argument("--axis", default="min_md",
                    choices=sorted(k for k, v in AXIS_DIRECTION.items() if v != 0))
    ap.add_argument("--top-n", type=int, default=20)
    ap.add_argument("--master-pool", default="")
    ap.add_argument("--dispatch-tsv", default="")
    ap.add_argument("--sample", type=int, default=0,
                    help="if >0, sample this many files deterministically")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--json", action="store_true",
                    help="emit JSON to stdout instead of a table")
    args = ap.parse_args()

    rows = extract_worst_files(
        args.archive, args.axis, args.top_n,
        master_pool=args.master_pool or None,
        dispatch_tsv=args.dispatch_tsv or None,
        sample=args.sample,
        seed=args.seed,
    )
    if args.json:
        # Strip non-JSONable items
        out = []
        for r in rows:
            rr = {k: v for k, v in r.items() if k not in {"all_axes"}}
            rr["all_axes"] = {k: v for k, v in r["all_axes"].items()
                              if isinstance(v, (int, float, type(None)))}
            out.append(rr)
        json.dump(out, sys.stdout, indent=2, default=str)
        sys.stdout.write("\n")
        return
    print(f"WORST {len(rows)} files in '{args.archive}' along axis '{args.axis}'")
    print(f"{'rank':>4} {'value':>10} {'metal':>5} {'CN':>3} {'poly':>7} "
          f"{'refcode':>8} {'dispatch':>14} file")
    for i, r in enumerate(rows, 1):
        v = r["axis_value"]
        vs = f"{v:.3f}" if isinstance(v, float) else str(v)
        print(f"{i:4d} {vs:>10} {str(r['metal'] or ''):>5} "
              f"{str(r['cn'] or ''):>3} {str(r['polyhedron_class'] or ''):>7} "
              f"{str(r['refcode'] or ''):>8} "
              f"{str(r['dispatch_path'] or ''):>14} "
              f"{os.path.basename(r['file_path'])}")


if __name__ == "__main__":
    main()
