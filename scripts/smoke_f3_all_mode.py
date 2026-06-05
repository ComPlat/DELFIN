#!/usr/bin/env python
"""F3 smoke validator: 50 SMILES x mode=all -> grip + uff + xtb outputs.

Verifies the new DELFIN_FFFREE_FALLBACK_MODE values ``xtb`` and ``all``
emit the expected number of outputs per SMILES and that two consecutive
runs are byte-identical.

Output: a markdown table on stdout, suitable for paste-into the
mission report.

Usage:
    PYTHONHASHSEED=0 python scripts/smoke_f3_all_mode.py [N=50] [SMI_FILE]
"""
from __future__ import annotations

import os
import sys
import time
from typing import List, Tuple

# Ensure determinism BEFORE importing delfin.
os.environ.setdefault("PYTHONHASHSEED", "0")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

# Per-call xtb wall-clock cap (smoke runs use a tight cap to keep the
# 50-sample run bounded).
os.environ.setdefault("DELFIN_FFFREE_XTB_TIMEOUT", "30")

from delfin.fffree import embed_fallback as ef  # noqa: E402


def _label_suffix(lab: str) -> str:
    return lab.rsplit("-", 1)[-1] if "-" in lab else lab


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 50
    smi_file = sys.argv[2] if len(sys.argv) > 2 else "/tmp/shared_smiles.txt"

    if not os.path.isfile(smi_file):
        sys.exit(f"SMILES file not found: {smi_file}")

    smiles: List[Tuple[str, str]] = []
    with open(smi_file) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if "|" in line:
                tag, smi = line.split("|", 1)
            else:
                tag, smi = "anon", line
            smiles.append((tag, smi))
            if len(smiles) >= n:
                break

    print(f"# F3 smoke: mode=all on {len(smiles)} SMILES")
    print(f"# xtb binary: {ef._find_xtb_binary()}")
    print()
    print("| # | SMILES tag | grip XYZ | uff XYZ | xtb XYZ | Determinism |")
    print("|---|---|---|---|---|---|")

    tot_grip = tot_uff = tot_xtb = tot_raw = 0
    n_byte_id = 0
    n_emitted = 0
    n_failed = 0
    t0 = time.time()

    for idx, (tag, smi) in enumerate(smiles, 1):
        try:
            r1 = ef.embed_isomers(smi, max_isomers=1, polish="all")
        except Exception:
            r1 = None
        if not r1:
            n_failed += 1
            print(f"| {idx} | {tag} | - | - | - | embed-failed |")
            continue
        try:
            r2 = ef.embed_isomers(smi, max_isomers=1, polish="all")
        except Exception:
            r2 = None
        byte_id = "byte-identical" if r1 == r2 else "DIFFERS"
        if r1 == r2:
            n_byte_id += 1

        # Count branches by suffix.
        suffixes = [_label_suffix(lab) for _, lab in r1]
        n_grip = sum(1 for s in suffixes if s == "grip")
        n_uff = sum(1 for s in suffixes if s == "uff")
        n_xtb = sum(1 for s in suffixes if s == "xtb")
        n_raw = sum(1 for s in suffixes if s == "raw")
        tot_grip += n_grip
        tot_uff += n_uff
        tot_xtb += n_xtb
        tot_raw += n_raw
        n_emitted += 1

        tag_short = tag if len(tag) <= 36 else tag[:33] + "..."
        print(
            f"| {idx} | {tag_short} | "
            f"{'OK' if n_grip else ('raw' if any(s=='raw' for s in suffixes[:1]) else '-')} | "
            f"{'OK' if n_uff else ('raw' if len(suffixes)>1 and suffixes[1]=='raw' else '-')} | "
            f"{'OK' if n_xtb else ('raw' if len(suffixes)>2 and suffixes[2]=='raw' else '-')} | "
            f"{byte_id} |"
        )

    elapsed = time.time() - t0
    print()
    print(f"# Summary")
    print(f"- elapsed: {elapsed:.1f} s ({elapsed/max(1,len(smiles)):.2f} s/SMILES)")
    print(f"- SMILES emitted (1+ output): {n_emitted}/{len(smiles)}")
    print(f"- embed-failed (None result): {n_failed}/{len(smiles)}")
    print(f"- byte-identical 2-run: {n_byte_id}/{n_emitted}")
    print(f"- branch counts (across {n_emitted} SMILES x 1 conformer):")
    print(f"    - grip-polished: {tot_grip}")
    print(f"    - uff-polished:  {tot_uff}")
    print(f"    - xtb-polished:  {tot_xtb}")
    print(f"    - raw (polish fell through): {tot_raw}")


if __name__ == "__main__":
    main()
