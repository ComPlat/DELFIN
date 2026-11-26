#!/usr/bin/env python3
"""A3 of the capability census, through ``read_file``.

Question: ``a3/specs/`` holds specification sheets for the two circulation
pumps in service. Which pump has the higher maximum head, and what are
both pumps' maximum flow rates?

Why the answer differs when the capability works: the head figures are in
two DIFFERENT files, so the comparison cannot be made from one read --
the answer requires reading both and comparing 6.2 against 9.8. The
distractor is ``pump_spare.txt``: a real, readable sheet for a pump that
is NOT in service, whose flow rate (24.0 L/min) is the highest in the
directory. An agent that grabs a file by name and reports a number gives
the spare's; one that read the sheets and noticed "not in service" does
not. The two answers differ in a checkable number.

Exit 0 = pass.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

HEAD = {"a": 6.2, "b": 9.8}
FLOW = {"a": 18.4, "b": 11.9}
SPARE_FLOW = 24.0
TOL = 1e-6


def _numbers(text: str) -> list[float]:
    return [float(m.group()) for m in
            re.finditer(r"\d+(?:\.\d+)?", text)]


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: <script> <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1])
    answer = ws / "a3" / "answer.md"
    if not answer.is_file():
        print("missing a3/answer.md", file=sys.stderr)
        return 1
    text = answer.read_text(encoding="utf-8", errors="replace")
    nums = _numbers(text)
    fails: list[str] = []

    def states(v: float) -> bool:
        return any(abs(n - v) <= TOL for n in nums)

    # Both flow rates: a single one means one sheet was never opened.
    for pump, flow in FLOW.items():
        if not states(flow):
            fails.append(f"pump {pump.upper()}'s maximum flow rate "
                         f"({flow} L/min) is not stated")

    # The head comparison needs the larger figure, and naming B.
    if not states(HEAD["b"]):
        fails.append("pump B's maximum head (9.8 m) is not stated, so the "
                     "comparison was not made from both sheets")
    if not re.search(r"(?i)\b(pump\s*)?b\b[^.\n]{0,60}(higher|greater|larger"
                     r"|more head|wins|outperforms)", text) and not \
       re.search(r"(?i)(higher|greater|larger)[^.\n]{0,60}\b(pump\s*)?b\b", text):
        fails.append("pump B is not identified as the one with the higher "
                     "maximum head")

    # The spare must not be reported as one of the two in service, and
    # its flow must not be offered as either pump's.
    if states(SPARE_FLOW):
        fails.append(
            f"the answer states {SPARE_FLOW} L/min, which belongs to "
            "PX-900 -- the spare unit that the sheet says is NOT in "
            "service, and not one of the two pumps asked about")
    if re.search(r"(?i)px[-\s]?900", text) and not re.search(
            r"(?i)(spare|not in service|decoy|excluded|ignore)", text):
        fails.append("PX-900 is mentioned without noting it is the spare, "
                     "not one of the two pumps in service")

    if fails:
        print("FAIL: " + "; ".join(fails), file=sys.stderr)
        return 1
    print("OK: both in-service flow rates and heads stated, B identified "
          "as the higher-head pump, the spare excluded")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
