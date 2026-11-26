#!/usr/bin/env python3
"""A6 of the capability census, through ``search_calcs``.

Question: ``calc_archive/calc/`` holds finished calculations. Within each
method, which calculation is the lowest in energy, and by how much does the
best PBE0/def2-SVP calculation lie below the best B3LYP/def2-SVP one?

The corpus lives at ``calc_archive/calc/`` because that is the path
``benchmark_runner._fixture_calc_dirs`` redirects the calc tools at; built
anywhere else the fixture is invisible to the gate and the tools answer from
the operator's own tree. The answer is read from ``sv1/answer.md``, outside
the indexed corpus.

Why the answer differs when the capability works: there are three
ethanol calculations at two different levels of theory, plus a METHANOL
calculation that is much lower in absolute energy and would win any
ranking done across methods. A total energy compares only within one
method, so the correct answer ranks PBE0/def2-SVP against PBE0/def2-SVP,
reports the PBE0/def2-TZVP run separately, uses the methanol folder for
nothing, and computes the cross-method difference as 0.074092 Eh
(= 46.49 kcal/mol) rather than by mixing levels.

Exit 0 = pass.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

TOL = 1e-5
# Measured from the fixture setup: PBE0/SVP -154.981234, B3LYP/SVP
# -154.907142 -> 0.074092 Eh. 1 Eh = 627.509474 kcal/mol -> 46.49.
PBE0_SVP = -154.981234
B3LYP_SVP = -154.907142
DELTA_EH = 0.074092
DELTA_KCAL = 46.49


def _numbers(text: str) -> list[float]:
    return [float(m.group()) for m in
            re.finditer(r"-?\d+\.\d+", text)]


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: <script> <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1])
    answer = ws / "sv1" / "answer.md"
    if not answer.is_file():
        print("missing sv1/answer.md", file=sys.stderr)
        return 1
    text = answer.read_text(encoding="utf-8", errors="replace")
    nums = _numbers(text)
    fails: list[str] = []

    def states(v: float) -> bool:
        return any(abs(n - v) <= TOL for n in nums)

    # The two minima of the same basis set.
    if not states(PBE0_SVP):
        fails.append(f"the PBE0/def2-SVP energy {PBE0_SVP} Eh is not stated")
    if not states(B3LYP_SVP):
        fails.append(f"the B3LYP/def2-SVP energy {B3LYP_SVP} Eh is not stated")

    # The difference, in one of the two acceptable units.
    if not (states(DELTA_EH) or any(abs(n - DELTA_KCAL) <= 0.02 for n in nums)):
        fails.append(
            f"the cross-method difference is neither {DELTA_EH} Eh nor "
            f"{DELTA_KCAL} kcal/mol; a different value means the two energies "
            "were not compared correctly")

    # Levels of theory must be kept apart: the task's whole discipline is
    # that a total energy compares only within a method.
    if not re.search(r"(?i)(within (each|the same) method|per method|"
                     r"same (method|level)|not across|cannot be compared "
                     r"across)", text):
        fails.append("the answer never states that the ranking is within a "
                     "method, so the cross-method caveat is missing")

    # The methanol decoy must not appear as a ranked entry.
    if re.search(r"(?i)methanol", text) and not re.search(
            r"(?i)(different (molecule|species)|excluded|not ethanol|"
            r"ignore)", text):
        fails.append("the methanol calculation is treated as part of the "
                     "ethanol ranking; it is a different molecule")

    if fails:
        print("FAIL: " + "; ".join(fails), file=sys.stderr)
        return 1
    print("OK: per-method minima stated, the difference in a real unit, and "
          "the cross-method caveat made explicit")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
