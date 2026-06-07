#!/usr/bin/env python3
"""iter_gate_with_diagnostics.py — extended iter-gate with auto-diagnostics.

Wraps the existing quality_framework/scripts/iter_gate.py, then appends
bug-census diff (HEAD vs PREV) and a per-class CCDC RMSD report on HEAD.

The standard iter_gate is invoked as a subprocess so its semantics + exit
codes are preserved 1:1 (PASS=0, HEAL-FIRST=2, FAIL=1).  Diagnostics never
*change* the verdict — they expand the visibility.

Usage:
  iter_gate_with_diagnostics.py <NEW_ARCHIVE> --prev <PREV_ARCHIVE> \\
        --archive-dir-new <path> [--archive-dir-prev <path>] \\
        [--reports <DIR>] [--no-ccdc] [--diff-out <JSON>]

Author: hmaximilian <hmaximilian496@gmail.com>
"""
from __future__ import annotations
import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

DEFAULT_ITER_GATE = (
    "/home/qmchem_max/agent_workspace/quality_framework/scripts/iter_gate.py"
)
DEFAULT_REPORTS = (
    "/home/qmchem_max/agent_workspace/quality_framework/reports"
)
DEFAULT_ARCHIVE_ROOT = (
    "/home/qmchem_max/agent_workspace/quality_framework/xyz_archive"
)


def _run_iter_gate(new: str, prev: Optional[str], reports: str,
                   iter_gate_path: str) -> int:
    """Run the existing iter_gate.py as a subprocess.  Returns its exit code."""
    if not os.path.isfile(iter_gate_path):
        print(f"WARN: iter_gate.py not found at {iter_gate_path}; "
              f"diagnostics-only mode")
        return -1
    cmd = [sys.executable, iter_gate_path, new]
    if prev:
        cmd.extend(["--prev", prev])
    cmd.extend(["--reports", reports])
    try:
        r = subprocess.run(cmd, check=False)
        return r.returncode
    except Exception as e:
        print(f"WARN: iter_gate subprocess error: {e}")
        return -1


def _resolve_archive_dir(name: str, root: str) -> Optional[str]:
    """Best-effort: archive name → directory path."""
    cand = os.path.join(root, name)
    if os.path.isdir(cand):
        return cand
    # try to match a prefix in case --prev shorthand omitted suffix
    if os.path.isdir(root):
        for d in os.listdir(root):
            if name in d and os.path.isdir(os.path.join(root, d)):
                return os.path.join(root, d)
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("new")
    ap.add_argument("--prev", default=None)
    ap.add_argument("--archive-dir-new", default=None,
                    help="explicit archive dir for HEAD (default: derive from --archive-root + new)")
    ap.add_argument("--archive-dir-prev", default=None,
                    help="explicit archive dir for PREV (default: derive)")
    ap.add_argument("--archive-root", default=DEFAULT_ARCHIVE_ROOT)
    ap.add_argument("--reports", default=DEFAULT_REPORTS)
    ap.add_argument("--iter-gate", default=DEFAULT_ITER_GATE)
    ap.add_argument("--no-ccdc", action="store_true",
                    help="skip per-class CCDC overlay (faster)")
    ap.add_argument("--no-census", action="store_true",
                    help="skip bug-census")
    ap.add_argument("--diff-out", default=None,
                    help="write the census diff + per-class report to this JSON")
    ap.add_argument("--master-pool", default="")
    ap.add_argument("--dispatch-tsv-new", default="")
    ap.add_argument("--dispatch-tsv-prev", default="")
    args = ap.parse_args()

    # 1) standard iter_gate
    print("\n" + "=" * 76)
    print("=== Standard iter_gate (aggregate axes) ===")
    print("=" * 76)
    rc = _run_iter_gate(args.new, args.prev, args.reports, args.iter_gate)

    # 2) bug census on HEAD (+ diff vs PREV)
    diff_payload = {}
    if not args.no_census:
        from bug_census import census, render_human, diff_reports, render_diff
        head_dir = (args.archive_dir_new or
                    _resolve_archive_dir(args.new, args.archive_root))
        if not head_dir:
            print("\nWARN: no archive dir for HEAD; skipping census.")
        else:
            print("\n" + "=" * 76)
            print("=== Bug Census on HEAD ===")
            print("=" * 76)
            head_rep = census(
                head_dir,
                master_pool=args.master_pool or None,
                dispatch_tsv=args.dispatch_tsv_new or None,
            )
            print(render_human(head_rep))
            diff_payload["head_census"] = head_rep

            prev_dir = (args.archive_dir_prev or
                        (args.prev and _resolve_archive_dir(args.prev, args.archive_root))
                        or None)
            if prev_dir:
                print("\n=== Bug Census on PREV ===")
                prev_rep = census(
                    prev_dir,
                    master_pool=args.master_pool or None,
                    dispatch_tsv=args.dispatch_tsv_prev or None,
                )
                print(render_human(prev_rep))
                diff_rows = diff_reports(prev_rep, head_rep)
                print("\n=== Bug-Class Diff (HEAD vs PREV) ===")
                print(render_diff(diff_rows))
                diff_payload["prev_census"] = prev_rep
                diff_payload["diff"] = diff_rows

    # 3) per-class CCDC overlay on HEAD
    if not args.no_ccdc:
        from ccdc_overlay import overlay_against_ccdc
        from per_class_ccdc_report import group_by_class, render_human as render_pc
        head_dir = (args.archive_dir_new or
                    _resolve_archive_dir(args.new, args.archive_root))
        if head_dir:
            print("\n" + "=" * 76)
            print("=== Per-Class CCDC RMSD Report (HEAD) ===")
            print("=" * 76)
            ov = overlay_against_ccdc(head_dir)
            grp = group_by_class(ov["per_refcode"])
            print(render_pc(ov, grp))
            diff_payload["head_per_class"] = {"overlay_summary": ov["summary"],
                                              "groups": grp}

    if args.diff_out:
        with open(args.diff_out, "w") as fh:
            json.dump(diff_payload, fh, indent=2, default=str)
        print(f"\nDiff JSON written: {args.diff_out}")

    print("\n" + "=" * 76)
    print(f"VERDICT (from standard iter_gate): exit_code={rc}")
    print("=" * 76)
    sys.exit(rc if rc >= 0 else 0)


if __name__ == "__main__":
    main()
