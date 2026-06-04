"""scripts/grace_determinism_proof.py — Task D: determinism + reproducibility.

Two-run byte-identity check for GRACE on 50 SMILES drawn deterministically
from the 200-SMILES set (every 4th entry).  Compares
``_result_to_tuple(r1) == _result_to_tuple(r2)`` for each SMILES.

Also performs a cache-clear simulation: between the two runs, optionally
clear the Mogul library cache (the only cache GRACE touches) and verify
the output is still byte-identical.

Output: ``paper_data/grace_determinism_proof.txt`` with PASS/FAIL counts.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sys
import time
from typing import Any, Dict, List, Tuple

os.environ.setdefault("PYTHONHASHSEED", "0")
os.environ.setdefault("DELFIN_FFFREE_GRACE_ENABLE", "1")

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, os.pardir))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from delfin.fffree.grace_ensemble import grace_enumerate  # noqa: E402
from scripts.grace_benchmark_200_curation import all_entries  # noqa: E402

try:
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
except Exception:
    pass


def _result_to_tuple(r) -> Tuple:
    """Bytes-comparable representation of a GraceResult, timing excluded."""
    flat: List[Tuple] = []
    for iso_id in sorted(r.per_isomer.keys()):
        for c in r.per_isomer[iso_id]:
            flat.append((
                int(iso_id),
                int(c.ring_id),
                int(c.rotamer_id),
                tuple(round(float(x), 6) for x in c.rotamer_offsets),
                str(c.label),
                str(c.method),
                int(c.clash_count),
                round(float(c.severity), 5),
                round(float(c.cshm), 5),
                round(float(c.score), 5),
                tuple(c.syms),
                # Coords rounded to 5 decimal places so floating-point noise
                # does not invalidate the equality test.
                tuple(tuple(round(float(v), 5) for v in row) for row in c.P),
            ))
    return tuple(flat)


def _hash_tuple(t: Tuple) -> str:
    """SHA-256 of the JSON of the tuple; used as a compact digest in the report."""
    return hashlib.sha256(
        json.dumps(t, sort_keys=True, default=str).encode("utf-8")
    ).hexdigest()[:16]


def _clear_cache_dirs():
    """Remove any user-cache directories GRACE / Mogul touch."""
    candidates = [
        os.path.expanduser("~/.cache/mogul"),
        os.path.expanduser("~/.cache/delfin"),
        os.path.expanduser("~/.cache/grace"),
    ]
    for p in candidates:
        if os.path.isdir(p):
            try:
                shutil.rmtree(p)
            except Exception:
                pass


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=os.path.join(_ROOT, "paper_data",
                                                  "grace_determinism_proof.txt"))
    ap.add_argument("--n", type=int, default=50,
                     help="Number of SMILES to sample (every k-th entry).")
    ap.add_argument("--budget-s", type=float, default=15.0)
    ap.add_argument("--max-per-isomer", type=int, default=5)
    ap.add_argument("--max-isomers", type=int, default=6)
    ap.add_argument("--max-rotamer-tuples", type=int, default=9)
    ap.add_argument("--clear-cache-between", action="store_true",
                     help="Clear ~/.cache/{mogul,delfin,grace} between the "
                          "two runs.")
    args = ap.parse_args()

    entries = all_entries()
    # Deterministic sample: every k-th entry so we get a diverse subset.
    step = max(1, len(entries) // args.n)
    sample = entries[::step][: args.n]

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fh = open(args.out, "w")
    fh.write(f"# GRACE determinism proof — n={len(sample)}  "
             f"method=lm  budget={args.budget_s}s  clear_cache={args.clear_cache_between}\n")
    fh.write("# " + "=" * 76 + "\n")
    fh.write(f"# Determinism contract: same SMILES + same env -> byte-identical "
             f"per-candidate (iso_id, ring_id, rot_id, label, method, "
             f"clash, severity[5dp], cshm[5dp], score[5dp], syms, coords[5dp]).\n")
    fh.write("#\n")

    pass_count = 0
    fail_count = 0
    skip_count = 0  # SMILES that don't decompose / emit nothing

    t0 = time.monotonic()
    for i, e in enumerate(sample, 1):
        kwargs = dict(
            method="lm",
            max_per_isomer=args.max_per_isomer,
            max_isomers=args.max_isomers,
            max_rotamer_tuples=args.max_rotamer_tuples,
            max_total_time_s=args.budget_s,
            max_per_ring=2,
        )
        try:
            r1 = grace_enumerate(e.smiles, **kwargs)
        except Exception as exc:
            fh.write(f"  [{i:>3}/{len(sample)}] EXC1  {e.name[:32]:<32}  {exc!r}\n")
            fail_count += 1
            continue

        if args.clear_cache_between:
            _clear_cache_dirs()

        try:
            r2 = grace_enumerate(e.smiles, **kwargs)
        except Exception as exc:
            fh.write(f"  [{i:>3}/{len(sample)}] EXC2  {e.name[:32]:<32}  {exc!r}\n")
            fail_count += 1
            continue

        t1 = _result_to_tuple(r1)
        t2 = _result_to_tuple(r2)
        h1 = _hash_tuple(t1)
        h2 = _hash_tuple(t2)
        if r1.n_emitted == 0 and r2.n_emitted == 0:
            # Both empty: byte-identical empty tuple is a trivial pass.
            skip_count += 1
            verdict = "SKIP"
        else:
            if t1 == t2:
                pass_count += 1
                verdict = "PASS"
            else:
                fail_count += 1
                verdict = "FAIL"
        fh.write(f"  [{i:>3}/{len(sample)}] {verdict}  "
                 f"{e.name[:32]:<32}  emit={r1.n_emitted:>3}/{r2.n_emitted:<3}  "
                 f"h1={h1} h2={h2}\n")
        fh.flush()

    dt = time.monotonic() - t0
    n_runnable = pass_count + fail_count
    pct = 100.0 * pass_count / max(n_runnable, 1)
    fh.write("# " + "=" * 76 + "\n")
    fh.write(f"# TOTAL    n={len(sample)}  PASS={pass_count}  FAIL={fail_count}  "
             f"SKIP={skip_count}  (runnable={n_runnable})\n")
    fh.write(f"# PCT_PASS {pct:.1f}%  (over runnable subset)\n")
    fh.write(f"# T_TOTAL  {dt:.1f}s\n")
    fh.close()

    print(f"# Determinism proof: PASS={pass_count}/{n_runnable} "
          f"({pct:.1f}%)  FAIL={fail_count}  SKIP={skip_count}  "
          f"t={dt:.1f}s  out={args.out}")


if __name__ == "__main__":
    main()
