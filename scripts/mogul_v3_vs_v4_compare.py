"""Compare v3 vs v4 anomaly detector on a sample of XYZ files.

Reports per-axis severity stats and counts of:
 - records present in v3 but eliminated by v4 (torsion-multimodal)
 - records present in v4 only (X-H, M-D angles — new channels)
 - records shared between v3 and v4
"""
from __future__ import annotations

import argparse
import json
import os
import random
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, "/home/qmchem_max/ComPlat/DELFIN")
sys.path.insert(0, "/home/qmchem_max/ComPlat/DELFIN/scripts")

import delfin._bond_decollapse as _bd  # noqa: E402
import delfin.fffree.mogul_detector_v3 as v3  # noqa: E402
import delfin.fffree.mogul_detector_v4 as v4  # noqa: E402
from mogul_v4_forensik import parse_xyz  # noqa: E402


def fingerprint(r):
    """Stable record key — (axis, sorted atoms, center)."""
    return (r["axis"], tuple(sorted(r["atoms"])), r.get("center"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--pool",
        default="/home/qmchem_max/agent_workspace/quality_framework/xyz_archive/2792332-aromatic-symmetry-VOLLPOOL",
    )
    ap.add_argument("--sample", type=int, default=200)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument(
        "--out",
        default="/home/qmchem_max/ComPlat/DELFIN/paper_data/mogul_v3_vs_v4_compare.json",
    )
    ap.add_argument(
        "--summary",
        default="/home/qmchem_max/ComPlat/DELFIN/paper_data/mogul_v3_vs_v4_compare.md",
    )
    args = ap.parse_args()

    os.environ.setdefault("PYTHONHASHSEED", "0")
    random.seed(args.seed)
    np.random.seed(args.seed)

    pool = Path(args.pool)
    files = sorted([p for p in pool.iterdir() if p.suffix == ".xyz"])
    if len(files) > args.sample:
        files = random.sample(files, args.sample)
    files.sort()

    # warm v5
    print(f"Loading v5 lib ...", flush=True)
    v4.load_v5_lookup()
    idx = v3._legacy_load_index()
    print(f"Fragment index size: {len(idx)}", flush=True)

    print(f"Processing {len(files)} files", flush=True)
    t0 = time.time()

    sum_v3 = 0
    sum_v4 = 0
    sum_only_v3 = 0
    sum_only_v4 = 0
    sum_shared = 0
    axis_v3 = Counter()
    axis_v4 = Counter()
    source_v4 = Counter()
    only_v3_axis = Counter()
    only_v4_source = Counter()
    sev_only_v3 = []
    sev_only_v4 = []
    sev_shared_old = []
    sev_shared_new = []
    n_ok = 0
    n_err = 0

    for i, f in enumerate(files):
        try:
            syms, P = parse_xyz(str(f))
            if len(syms) < 3:
                n_err += 1
                continue
            r3 = v3.detect_anomalies_v3(syms, P, mode="ccdc", index=idx)
            r4 = v4.detect_anomalies_v4(syms, P, use_v4=True, index=idx)
        except Exception as e:
            n_err += 1
            continue
        n_ok += 1

        sum_v3 += len(r3)
        sum_v4 += len(r4)
        for r in r3:
            axis_v3[r["axis"]] += 1
        for r in r4:
            axis_v4[r["axis"]] += 1
            source_v4[r.get("source", "v3")] += 1

        fps3 = {fingerprint(r): r for r in r3}
        fps4 = {fingerprint(r): r for r in r4}
        only3_keys = set(fps3) - set(fps4)
        only4_keys = set(fps4) - set(fps3)
        shared_keys = set(fps3) & set(fps4)

        sum_only_v3 += len(only3_keys)
        sum_only_v4 += len(only4_keys)
        sum_shared += len(shared_keys)

        for k in only3_keys:
            r = fps3[k]
            only_v3_axis[r["axis"]] += 1
            sev_only_v3.append(float(r.get("sev_mad", 0)))
        for k in only4_keys:
            r = fps4[k]
            only_v4_source[r.get("source", "?")] += 1
            sev_only_v4.append(float(r.get("sev_mad", 0)))
        for k in shared_keys:
            sev_shared_old.append(float(fps3[k].get("sev_mad", 0)))
            sev_shared_new.append(float(fps4[k].get("sev_mad", 0)))

        if (i + 1) % 50 == 0:
            print(f"  ...{i+1}/{len(files)} elapsed={time.time()-t0:.1f}s",
                  flush=True)

    def stats(arr):
        if not arr:
            return {"n": 0}
        a = np.asarray(arr, dtype=float)
        return {
            "n": int(len(a)),
            "mean": float(a.mean()),
            "median": float(np.median(a)),
            "p95": float(np.percentile(a, 95)),
            "max": float(a.max()),
        }

    summary = {
        "pool": str(pool),
        "sample_size": len(files),
        "n_ok": n_ok,
        "n_err": n_err,
        "elapsed_s": round(time.time() - t0, 1),
        "total_v3": sum_v3,
        "total_v4": sum_v4,
        "only_v3": sum_only_v3,
        "only_v4": sum_only_v4,
        "shared": sum_shared,
        "axis_v3": dict(axis_v3),
        "axis_v4": dict(axis_v4),
        "source_v4": dict(source_v4),
        "only_v3_axis": dict(only_v3_axis),
        "only_v4_source": dict(only_v4_source),
        "sev_only_v3": stats(sev_only_v3),
        "sev_only_v4": stats(sev_only_v4),
        "sev_shared_old": stats(sev_shared_old),
        "sev_shared_new": stats(sev_shared_new),
    }

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(summary, fh, indent=2)

    md = ["# Mogul v3 vs v4 Comparison", ""]
    md.append(f"- pool: `{pool.name}`")
    md.append(f"- sample size: **{summary['sample_size']}** "
              f"({n_ok} ok, {n_err} err)")
    md.append(f"- elapsed: {summary['elapsed_s']} s")
    md.append("")
    md.append("## Total anomaly counts")
    md.append(f"- v3: **{sum_v3}**")
    md.append(f"- v4: **{sum_v4}**")
    md.append(f"- shared (v3 ∩ v4): {sum_shared}")
    md.append(f"- only in v3 (v4 eliminates): **{sum_only_v3}**")
    md.append(f"- only in v4 (new detections): **{sum_only_v4}**")
    md.append("")
    md.append("## Axis breakdown")
    md.append(f"- v3: {dict(axis_v3)}")
    md.append(f"- v4: {dict(axis_v4)}")
    md.append("")
    md.append("## v4 record sources")
    md.append(f"  {dict(source_v4)}")
    md.append("")
    md.append("## Records dropped by v4 (per axis)")
    md.append(f"  {dict(only_v3_axis)}")
    md.append("")
    md.append("## Records added by v4 (per source)")
    md.append(f"  {dict(only_v4_source)}")
    md.append("")
    md.append("## Severity statistics")
    md.append(f"- v3-only (dropped by v4): {summary['sev_only_v3']}")
    md.append(f"- v4-only (new by v4): {summary['sev_only_v4']}")
    md.append(f"- shared, old sev: {summary['sev_shared_old']}")
    md.append(f"- shared, new sev: {summary['sev_shared_new']}")

    Path(args.summary).write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
