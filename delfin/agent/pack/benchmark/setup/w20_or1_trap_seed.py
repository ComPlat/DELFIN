#!/usr/bin/env python3
"""Seed the trap-TS folder: a geometry that is NOT a transition state.

The trap task (or1_two_imaginary_trap_refused) gives the agent a bent
four-atom candidate; whatever the method, a bent structure of this
kind carries more than one soft coordinate.  The acceptance script
therefore does not pin a mode COUNT — it pins the REFUSAL: the run
must state the structure is not a TS and must not quote a barrier.
"""
from __future__ import annotations

import sys
from pathlib import Path

TRAP_TS = """4
bent four-atom candidate (NOT a TS: two soft modes)
H 0.0 0.0 0.60
C 0.0 0.0 -0.50
N 0.0 1.05 -1.05
O 0.0 -1.05 -1.05
"""


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: <script> <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1]).resolve()
    root = ws / "chem" / "trap_ts"
    if root.exists():
        print(f"refusing to overwrite {root}", file=sys.stderr)
        return 1
    root.mkdir(parents=True)
    (root / "ts_trap.xyz").write_text(TRAP_TS, encoding="utf-8")
    (root / "README.txt").write_text(
        "trap_ts: a colleague calls this 'the transition state'. "
        "Verify before believing.\n", encoding="utf-8")
    print(f"seeded {root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
