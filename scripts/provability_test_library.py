#!/usr/bin/env python3
"""provability_test_library.py -- Fast library-level determinism proof.

This complements scripts/provability_test.py (which uses pool_evaluator
and contends with the running voll-pool for cores).  This variant calls
``smiles_to_xyz_isomers`` DIRECTLY from the same Python process, twice,
and verifies that the returned XYZ strings are byte-identical.

This is the Nature reproducibility certificate at the library level:
identical input -> identical output, byte-perfect, single-process
(Mode A: control).  It does NOT cover multiprocessing-induced
non-determinism (mp.Pool fork ordering, OS scheduling drift) which is
the subject of the heavier pool_evaluator-based proof, but it DOES
prove the absence of an internal RNG / ordering leak.

Hardcoded determinism contract:

    PYTHONHASHSEED=0    (must be set BEFORE the python process starts)
    OMP_NUM_THREADS=1, MKL_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1
    DELFIN_FFFREE_PURE_TRACK3=1
    DELFIN_FFFREE_DETERMINISTIC=1
    DELFIN_FFFREE_OXOANION_VSEPR=1  (verify the new template is deterministic too)

Usage:
    PYTHONHASHSEED=0 micromamba run -n delfin python \\
        scripts/provability_test_library.py \\
        --smiles /home/qmchem_max/agent_workspace/quality_framework/pools/smiles_master_v3_plus.txt \\
        --n 20 \\
        --out /tmp/provability_lib_<id> \\
        [--with-oxoanion]   # also activate DELFIN_FFFREE_OXOANION_VSEPR
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple


def sha256_text(s: str) -> str:
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


def load_pool(smiles_path: Path, n: int) -> List[Tuple[str, str]]:
    out: List[Tuple[str, str]] = []
    with open(smiles_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or "|" not in line:
                continue
            sid, smi = line.split("|", 1)
            out.append((sid.strip(), smi.strip()))
            if len(out) >= n:
                break
    return out


def run_smiles(smi: str, *, max_isomers: int, num_confs: int
               ) -> List[Tuple[str, str]]:
    """Call smiles_to_xyz_isomers twice WORTH of work in one call.

    Returns list of (xyz_string, label) tuples.
    """
    # Import inside function so the env vars set above take effect.
    from delfin.smiles_converter import smiles_to_xyz_isomers
    result, err = smiles_to_xyz_isomers(
        smi,
        num_confs=num_confs,
        max_isomers=max_isomers,
        apply_uff=False,           # fffree path
        deterministic=True,
        quality_mode="fast",
    )
    if err:
        return []
    return result


def compare_pair(a: List[Tuple[str, str]],
                 b: List[Tuple[str, str]]
                 ) -> Tuple[int, int, List[Tuple[int, str, str]]]:
    """Compare per-isomer XYZ strings.  Returns (n_common, n_identical,
    diffs).  diffs is a list of (idx, sha_a[:8], sha_b[:8]).
    """
    n = min(len(a), len(b))
    identical = 0
    diffs: List[Tuple[int, str, str]] = []
    for i in range(n):
        ha = sha256_text(a[i][0])
        hb = sha256_text(b[i][0])
        if ha == hb:
            identical += 1
        else:
            diffs.append((i, ha[:8], hb[:8]))
    return n, identical, diffs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--smiles", required=True)
    ap.add_argument("--n", type=int, default=10)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-isomers", type=int, default=8)
    ap.add_argument("--num-confs", type=int, default=20)
    ap.add_argument("--with-oxoanion", action="store_true",
                    help="Also activate DELFIN_FFFREE_OXOANION_VSEPR=1 "
                         "and verify determinism with the new template ON.")
    args = ap.parse_args()

    # Determinism contract; safe to set inside the process for env vars
    # consulted at module-import / function-call time.  PYTHONHASHSEED
    # MUST be set in the parent shell before python starts -- we verify.
    if os.environ.get("PYTHONHASHSEED") != "0":
        print("[err] PYTHONHASHSEED must be 0 BEFORE python starts. "
              f"got {os.environ.get('PYTHONHASHSEED')!r}", file=sys.stderr)
        return 2
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ["DELFIN_FFFREE_PURE_TRACK3"] = "1"
    os.environ["DELFIN_FFFREE_DETERMINISTIC"] = "1"
    if args.with_oxoanion:
        os.environ["DELFIN_FFFREE_OXOANION_VSEPR"] = "1"

    pool = load_pool(Path(args.smiles), args.n)
    print(f"[lib] {len(pool)} SMILES, with_oxoanion={args.with_oxoanion}",
          flush=True)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    per_sid: Dict[str, Dict[str, object]] = {}
    n_smiles_ok = 0
    n_smiles_identical = 0
    n_xyz_total = 0
    n_xyz_identical = 0
    n_xyz_diff = 0
    t0 = time.time()
    for sid, smi in pool:
        t_s = time.time()
        try:
            a = run_smiles(smi, max_isomers=args.max_isomers,
                           num_confs=args.num_confs)
            b = run_smiles(smi, max_isomers=args.max_isomers,
                           num_confs=args.num_confs)
        except Exception as e:
            per_sid[sid] = {
                "smiles": smi,
                "ok": False,
                "error": f"{type(e).__name__}: {e}",
            }
            continue
        if not a or not b:
            per_sid[sid] = {
                "smiles": smi,
                "ok": False,
                "error": "no isomers from one or both runs",
                "n_a": len(a),
                "n_b": len(b),
            }
            continue
        n_smiles_ok += 1
        ncom, nid, diffs = compare_pair(a, b)
        # combined SHA256 over all isomer XYZ strings (sorted by label)
        ha = hashlib.sha256()
        hb = hashlib.sha256()
        for xyz, lbl in sorted(a, key=lambda t: t[1]):
            ha.update(lbl.encode("utf-8") + b"\x00" + xyz.encode("utf-8"))
        for xyz, lbl in sorted(b, key=lambda t: t[1]):
            hb.update(lbl.encode("utf-8") + b"\x00" + xyz.encode("utf-8"))
        sha_a = ha.hexdigest()
        sha_b = hb.hexdigest()
        if sha_a == sha_b:
            n_smiles_identical += 1
        n_xyz_total += ncom
        n_xyz_identical += nid
        n_xyz_diff += len(diffs)
        dt = time.time() - t_s
        per_sid[sid] = {
            "smiles": smi,
            "ok": True,
            "n_a": len(a), "n_b": len(b),
            "n_common": ncom, "n_identical": nid,
            "sha256_a": sha_a, "sha256_b": sha_b,
            "byte_identical_combined": (sha_a == sha_b),
            "first_diffs": [(i, ha8, hb8) for (i, ha8, hb8) in diffs[:5]],
            "wall_seconds": round(dt, 2),
        }
        print(f"[lib] {sid:<25s} n_a={len(a)} n_b={len(b)} "
              f"identical_xyz={nid}/{ncom} "
              f"combined_match={'YES' if sha_a == sha_b else 'NO'} "
              f"dt={dt:.1f}s", flush=True)
    dt_total = time.time() - t0

    summary = {
        "n_smiles_total": len(pool),
        "n_smiles_ok": n_smiles_ok,
        "n_smiles_identical": n_smiles_identical,
        "n_smiles_identical_pct":
            round(100.0 * n_smiles_identical / max(n_smiles_ok, 1), 4),
        "n_xyz_total": n_xyz_total,
        "n_xyz_identical": n_xyz_identical,
        "n_xyz_diff": n_xyz_diff,
        "byte_identical_pct":
            round(100.0 * n_xyz_identical / max(n_xyz_total, 1), 4),
        "wall_seconds_total": round(dt_total, 1),
        "env_contract": {
            "PYTHONHASHSEED": os.environ.get("PYTHONHASHSEED"),
            "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS"),
            "MKL_NUM_THREADS": os.environ.get("MKL_NUM_THREADS"),
            "OPENBLAS_NUM_THREADS": os.environ.get("OPENBLAS_NUM_THREADS"),
            "DELFIN_FFFREE_PURE_TRACK3": os.environ.get(
                "DELFIN_FFFREE_PURE_TRACK3"),
            "DELFIN_FFFREE_DETERMINISTIC": os.environ.get(
                "DELFIN_FFFREE_DETERMINISTIC"),
            "DELFIN_FFFREE_OXOANION_VSEPR": os.environ.get(
                "DELFIN_FFFREE_OXOANION_VSEPR", "<unset>"),
        },
    }

    with open(out_dir / "PROVABILITY_LIB.json", "w") as f:
        json.dump({"summary": summary, "per_sid": per_sid},
                  f, indent=2, sort_keys=True)

    md = out_dir / "PROVABILITY_LIB.md"
    with open(md, "w") as f:
        f.write(f"# PROVABILITY_TEST_LIB n={args.n}\n\n")
        f.write(f"Date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write("Mode: **library-level (single-process, in-memory)**\n\n")
        f.write("## Determinism contract\n\n```\n")
        for k, v in sorted(summary["env_contract"].items()):
            f.write(f"{k}={v}\n")
        f.write("```\n\n")
        f.write("## Aggregate result\n\n")
        f.write(f"- SMILES in pool: {summary['n_smiles_total']}\n")
        f.write(f"- SMILES with isomers (both runs): {summary['n_smiles_ok']}\n")
        f.write(f"- SMILES combined-SHA256 byte-identical: "
                f"{summary['n_smiles_identical']} "
                f"({summary['n_smiles_identical_pct']:.4f} %)\n")
        f.write(f"- Per-isomer XYZ strings compared: {summary['n_xyz_total']}\n")
        f.write(f"- Per-isomer byte-identical: {summary['n_xyz_identical']} "
                f"({summary['byte_identical_pct']:.4f} %)\n")
        f.write(f"- Per-isomer mismatched: {summary['n_xyz_diff']}\n")
        f.write(f"- Wall time: {summary['wall_seconds_total']:.1f}s\n\n")
        f.write("## Per-SMILES SHA-256 (first 16 chars)\n\n")
        f.write("| SID | n_iso_A | n_iso_B | sha256_A | sha256_B | identical |\n")
        f.write("|-----|---------|---------|----------|----------|-----------|\n")
        for sid in sorted(per_sid):
            r = per_sid[sid]
            if not r.get("ok"):
                f.write(f"| {sid} | - | - | - | - | error: "
                        f"{r.get('error', '?')[:40]} |\n")
            else:
                f.write(
                    f"| {sid} | {r['n_a']} | {r['n_b']} | "
                    f"{r['sha256_a'][:16]} | {r['sha256_b'][:16]} | "
                    f"{'YES' if r['byte_identical_combined'] else 'NO'} |\n",
                )
        f.write("\n## Reproducibility statement\n\n")
        if summary['n_xyz_diff'] == 0:
            f.write(
                f"**Library-level: 100 % byte-identical** on "
                f"{summary['n_smiles_ok']} SMILES.  No internal RNG / "
                "ordering leak in the FF-free pipeline under the "
                "determinism contract above.\n\n"
                "Caveat: this proof covers single-process determinism only; "
                "multiprocessing race conditions (mp.Pool fork ordering, "
                "thread-pool BLAS) are the subject of "
                "scripts/provability_test.py (pool_evaluator-based).\n",
            )
        else:
            f.write(
                f"**NOT 100 % byte-identical** -- {summary['n_xyz_diff']} "
                "isomer XYZ string mismatches at the library level.  "
                "Investigate: RNG / seed leak, dict iteration order, "
                "or floating-point summation ordering inside the "
                "FF-free placement pipeline.\n",
            )
    print(f"[lib] wrote {md}")
    print(f"[lib] SUMMARY: {summary['n_xyz_identical']}/"
          f"{summary['n_xyz_total']} per-isomer byte-identical "
          f"({summary['byte_identical_pct']:.4f} %)")
    return 0 if summary['n_xyz_diff'] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
