#!/usr/bin/env python3
"""bug_census.py — full archive bug census + per-class breakdown.

Walks every .xyz file in an archive, runs classify_bug_pattern on each, and
aggregates the results by bug class.  Produces a human-readable report AND a
machine-readable JSON for the iter-gate diagnostics.

Public API:
  census(archive, *, master_pool=None, dispatch_tsv=None, sample=0, seed=0)
    -> Dict[bug_class_name -> {
         count, frac_files, example_files, affected_classes, dispatch_paths
       }]

Author: hmaximilian <hmaximilian496@gmail.com>
"""
from __future__ import annotations
import argparse
import json
import os
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from diagnostic_common import (  # noqa: E402
    list_archive_xyz, summarise_structure, load_manifest_pool,
    load_dispatch_forensik, label_from_filename, refcode_from_filename,
)
from classify_bug_pattern import classify_atoms, KNOWN_PATTERNS  # noqa: E402


# ---------------------------------------------------------------------------
def census(
    archive: str,
    *,
    master_pool: Optional[str] = None,
    dispatch_tsv: Optional[str] = None,
    sample: int = 0,
    seed: int = 0,
    example_n: int = 3,
) -> Dict:
    """Run classify_bug_pattern over every file in 'archive' + aggregate.

    Returns a dict with keys 'archive', 'n_files_scanned', 'n_files_clean',
    'n_files_with_bugs', and 'bug_classes' = {name → {count, ...}}.
    """
    paths = list_archive_xyz(archive)
    if sample and sample < len(paths):
        rng = random.Random(seed)
        paths = sorted(rng.sample(paths, sample))

    pool = load_manifest_pool(master_pool) if master_pool else {}
    dispatch = load_dispatch_forensik(dispatch_tsv) if dispatch_tsv else {}

    class_count: Counter = Counter()
    class_files: Dict[str, List[str]] = defaultdict(list)
    class_classes: Dict[str, Counter] = defaultdict(Counter)  # (metal,CN) → count
    class_dispatch: Dict[str, Counter] = defaultdict(Counter)
    n_clean = 0
    n_with_bugs = 0
    n_no_metal = 0
    n_parse_failed = 0

    for p in paths:
        rec = summarise_structure(p)
        atoms = rec.get("atoms") or []
        hits = classify_atoms(atoms)
        if not hits:
            n_clean += 1
            continue
        if hits == ["parse_failed"]:
            n_parse_failed += 1
            continue
        if hits == ["no_metal"]:
            n_no_metal += 1
            continue
        n_with_bugs += 1
        label = label_from_filename(p)
        smi = pool.get(label, "")
        disp = dispatch.get(smi, "unknown") if smi else "unknown"
        metal = rec.get("metal_elem") or "?"
        cn = rec.get("cn") if rec.get("cn") is not None else "?"
        cls_key = f"{metal},CN{cn}"
        for h in hits:
            class_count[h] += 1
            if len(class_files[h]) < example_n:
                class_files[h].append(
                    refcode_from_filename(p) or os.path.basename(p)
                )
            class_classes[h][cls_key] += 1
            class_dispatch[h][disp] += 1

    bug_classes: Dict[str, Dict] = {}
    total_scanned = len(paths)
    for name in sorted(class_count, key=lambda k: -class_count[k]):
        bug_classes[name] = {
            "count": class_count[name],
            "frac_files": (class_count[name] / total_scanned
                           if total_scanned else 0.0),
            "severity": KNOWN_PATTERNS.get(name, {}).get("severity", "?"),
            "description": KNOWN_PATTERNS.get(name, {}).get("description", ""),
            "example_files": list(class_files[name]),
            "affected_classes": class_classes[name].most_common(5),
            "dispatch_paths": dict(class_dispatch[name]),
        }
    return {
        "archive": archive,
        "n_files_scanned": total_scanned,
        "n_files_clean": n_clean,
        "n_files_with_bugs": n_with_bugs,
        "n_files_no_metal": n_no_metal,
        "n_files_parse_failed": n_parse_failed,
        "bug_classes": bug_classes,
    }


def _fmt_pct(num: int, denom: int) -> str:
    if not denom:
        return "n/a"
    return f"{100 * num / denom:.1f}%"


def render_human(report: Dict, *, max_classes: int = 30) -> str:
    out: List[str] = []
    arch = report["archive"]
    n_total = report["n_files_scanned"]
    out.append(f"=== Bug Census for {os.path.basename(arch.rstrip('/'))} ===")
    out.append(f"  scanned         : {n_total}")
    out.append(f"  clean           : {report['n_files_clean']} "
               f"({_fmt_pct(report['n_files_clean'], n_total)})")
    out.append(f"  with-bugs       : {report['n_files_with_bugs']} "
               f"({_fmt_pct(report['n_files_with_bugs'], n_total)})")
    out.append(f"  no-metal        : {report['n_files_no_metal']}")
    out.append(f"  parse-failed    : {report['n_files_parse_failed']}")
    out.append("")
    bugs = report["bug_classes"]
    if not bugs:
        out.append("  (no bug classes fired)")
        return "\n".join(out)
    out.append(f"{'Bug class':<35} {'Sev':<8} {'Count':>7} {'Frac':>7} "
               f"{'Examples':<26} {'Top class':<14} {'Top dispatch':<14}")
    out.append("-" * 120)
    for i, (name, info) in enumerate(bugs.items()):
        if i >= max_classes:
            break
        ex = ",".join(info["example_files"][:3])
        if len(ex) > 24:
            ex = ex[:24] + ".."
        # affected_classes[0] is (cls_str, count) or absent
        ac = info["affected_classes"]
        ac_top = f"{ac[0][0]} ({ac[0][1]})" if ac else ""
        if len(ac_top) > 12:
            ac_top = ac_top[:12] + ".."
        # dispatch
        dp = info["dispatch_paths"]
        if dp:
            top_dp = max(dp.items(), key=lambda kv: kv[1])
            tot = sum(dp.values())
            top_dp_s = f"{top_dp[0]} {100*top_dp[1]//tot}%"
        else:
            top_dp_s = ""
        out.append(f"{name:<35} {info['severity']:<8} {info['count']:>7d} "
                   f"{_fmt_pct(info['count'], n_total):>7} "
                   f"{ex:<26} {ac_top:<14} {top_dp_s:<14}")
    if len(bugs) > max_classes:
        out.append(f"... and {len(bugs) - max_classes} more classes")
    return "\n".join(out)


def diff_reports(prev: Dict, new: Dict) -> List[Dict]:
    """Per-class diff between two census reports."""
    prev_b = prev.get("bug_classes", {})
    new_b = new.get("bug_classes", {})
    all_classes = sorted(set(prev_b) | set(new_b))
    rows = []
    for k in all_classes:
        pc = prev_b.get(k, {}).get("count", 0)
        nc = new_b.get(k, {}).get("count", 0)
        if pc == 0 and nc == 0:
            continue
        delta = nc - pc
        rel = (delta / pc) if pc else None
        rows.append({
            "bug_class": k,
            "prev": pc,
            "new": nc,
            "delta": delta,
            "rel": rel,
        })
    rows.sort(key=lambda r: -abs(r["delta"]))
    return rows


def render_diff(rows: List[Dict]) -> str:
    out = [f"{'Bug class':<35} {'PREV':>8} {'NEW':>8} {'DELTA':>8} {'REL':>8} {'mark':<8}"]
    out.append("-" * 76)
    for r in rows:
        rel = r["rel"]
        rel_s = f"{rel*100:+.0f}%" if rel is not None else "n/a"
        if r["delta"] < 0:
            mark = "FIXED"
        elif r["delta"] == 0:
            mark = "(same)"
        else:
            mark = "WORSE"
        out.append(f"{r['bug_class']:<35} {r['prev']:>8d} {r['new']:>8d} "
                   f"{r['delta']:>+8d} {rel_s:>8} {mark:<8}")
    return "\n".join(out)


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("archive", help="xyz_archive/<name>/")
    ap.add_argument("--master-pool", default="")
    ap.add_argument("--dispatch-tsv", default="")
    ap.add_argument("--sample", type=int, default=0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out-json", default=None,
                    help="write the full report dict to this path")
    ap.add_argument("--diff", default=None,
                    help="path to a prior census JSON for diff report")
    args = ap.parse_args()

    rep = census(
        args.archive,
        master_pool=args.master_pool or None,
        dispatch_tsv=args.dispatch_tsv or None,
        sample=args.sample, seed=args.seed,
    )
    print(render_human(rep))
    if args.out_json:
        os.makedirs(os.path.dirname(os.path.abspath(args.out_json)) or ".",
                    exist_ok=True)
        with open(args.out_json, "w") as fh:
            json.dump(rep, fh, indent=2)
        print(f"\nJSON written: {args.out_json}")
    if args.diff:
        with open(args.diff) as fh:
            prev = json.load(fh)
        rows = diff_reports(prev, rep)
        print("\n=== Bug Census Diff vs PREV ===")
        print(render_diff(rows))


if __name__ == "__main__":
    main()
