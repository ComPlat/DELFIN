#!/usr/bin/env python3
"""Build anomaly atlas CSV from existing CCDC analysis cache.

Fallback for ``anomaly_pattern_atlas.py`` when the live CCDC run
times out before the parent finishes — reads the per-worker
``ccdc.json`` files left in the work-dir and aggregates them.

USAGE
    python scripts/anomaly_atlas_from_cache.py \\
        --work-dir /tmp/anomaly_atlas \\
        --out paper_data/anomaly_pattern_atlas.csv
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path


_FIX_DIRECTION_RULES = [
    (re.compile(r"\bM-(N|C|O|P|S|Cl|Br|I|F)\b"), "build",
     "Metal-donor bond — build-time hapto/coordination fix"),
    (re.compile(r"\beta[2-9]"), "build", "Hapto mode — build-time fix"),
    (re.compile(r"\bcarbonyl\b"), "polish", "Carbonyl backbone — polish"),
    (re.compile(r"\bring\b"), "polish", "Ring puckering / planarity — polish"),
    (re.compile(r"\baromatic\b"), "polish", "Aromatic flatness — polish"),
    (re.compile(r"^C-C$|^C=C$|^C#C$|^Ar-C$"), "polish",
     "C-C internal bond — polish"),
    (re.compile(r"^C-(N|O|S|P|F|Cl)$|^C=(N|O|S)$"), "polish",
     "C-X internal bond — polish"),
    (re.compile(r"^N-H$|^O-H$|^C-H$"), "polish",
     "X-H internal bond — polish"),
    (re.compile(r"\bL0\b|\bgeneralised\b", re.IGNORECASE), "library",
     "Generalised — COD library under-represented"),
]


def classify_fix(label: str):
    for pat, direction, reason in _FIX_DIRECTION_RULES:
        if pat.search(label or ""):
            return direction, reason
    return "polish", "default — internals polish"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--top-n", type=int, default=20)
    args = ap.parse_args()

    wdir = Path(args.work_dir)
    by_frag = defaultdict(lambda: {
        "zs": [], "n_unusual": 0, "n_total": 0,
        "example_files": [],
        "example_refcodes": [],
        "axis_counter": Counter(),
    })
    n_ok = 0
    for d in sorted(wdir.glob("w*")):
        cf = d / "ccdc.json"
        if not cf.exists() or cf.stat().st_size < 50:
            continue
        try:
            rec = json.loads(cf.read_text())
        except Exception:
            continue
        if not rec.get("ok"):
            continue
        n_ok += 1
        # try to infer original filename
        frame = d / "frame.xyz"
        fname = "unknown"
        if frame.exists():
            try:
                tx = frame.read_text(errors="ignore")
                m = re.search(r"smi=([^ \n]+)", tx)
                if m:
                    fname = m.group(1)
            except Exception:
                pass
        rcm = re.findall(r"[A-Z]{6}", fname)
        rc = rcm[-1] if rcm else None
        for f in rec.get("features", []):
            frag = (f.get("frag") or "<no-frag>").strip()
            axis = f.get("axis", "?")
            z = float(f.get("z", 0.0) or 0.0)
            d_acc = by_frag[(axis, frag)]
            d_acc["zs"].append(abs(z))
            d_acc["n_total"] += 1
            if f.get("unusual"):
                d_acc["n_unusual"] += 1
            d_acc["axis_counter"][axis] += 1
            if len(d_acc["example_files"]) < 5:
                d_acc["example_files"].append(fname)
                if rc and rc not in d_acc["example_refcodes"]:
                    d_acc["example_refcodes"].append(rc)
    print(f"aggregated from {n_ok} CCDC results", file=sys.stderr)

    rows = []
    for (axis, frag), d in by_frag.items():
        if d["n_total"] < 3:
            continue
        mean_abs_z = sum(d["zs"]) / len(d["zs"]) if d["zs"] else 0.0
        max_abs_z = max(d["zs"]) if d["zs"] else 0.0
        fix_dir, fix_reason = classify_fix(frag)
        rows.append({
            "axis": axis,
            "fragment": frag,
            "n_observations": d["n_total"],
            "n_unusual": d["n_unusual"],
            "frac_unusual": round(d["n_unusual"] / max(d["n_total"], 1), 3),
            "mean_abs_z": round(mean_abs_z, 3),
            "max_abs_z": round(max_abs_z, 3),
            "example_refcodes": ",".join(d["example_refcodes"][:5]),
            "example_files": "|".join(d["example_files"][:5]),
            "fix_direction": fix_dir,
            "fix_reason": fix_reason,
        })
    rows.sort(key=lambda r: -r["mean_abs_z"])
    rows_top = rows[:args.top_n]
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        if rows:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows_top)
    print(f"WROTE {out} ({len(rows_top)} rows from {len(rows)} classes)",
          file=sys.stderr)
    # also write full
    full = out.with_name(out.stem + "_full.csv")
    with full.open("w") as fh:
        if rows:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
    print(f"WROTE {full} (full {len(rows)} classes)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
