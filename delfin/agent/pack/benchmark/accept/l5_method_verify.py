#!/usr/bin/env python3
"""Acceptance for l5_name_method_with_basis_pairing.

Input: the workspace root (argv[1]). This task's answer is prose (the named
functional + basis pair) and is judged by expected_signals / forbidden_signals
on the trajectory. The accept action only guards fixture integrity:

* l5_method_notes/brief.txt must be UNCHANGED (it sets the context and names
  the constraints; editing the brief to inject the answer is a violation that
  must fail).
* the fixture must not have been stripped of its constraint text.

Exit 0 when intact, 1 otherwise; 2 on a missing argument.
"""
from __future__ import annotations

import sys
from pathlib import Path

_KEPT = ("Recommend ONE density-functional/basis-set combination for a "
         "closed-shell molecule with only H, C, N and O.")


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: l5_method_verify.py <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1]).resolve()
    brief = ws / "l5_method_notes" / "brief.txt"
    if not brief.is_file():
        print(f"FAIL: no {brief}", file=sys.stderr)
        return 1
    text = brief.read_text(encoding="utf-8")
    if _KEPT not in " ".join(text.split()):
        print("FAIL: brief.txt was altered", file=sys.stderr)
        return 1
    print("OK: method brief intact")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
