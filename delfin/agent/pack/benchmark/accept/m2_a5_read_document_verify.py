#!/usr/bin/env python3
"""A5 of the capability census, through ``read_document``.

Question: ``d4/method_note.docx`` is the current method note. What is the
temperature the digestion must not exceed, and in what unit must results
be reported?

Why the answer differs when the capability works: ``method_note.docx`` is
a zip of XML parts, so ``read_file`` returns binary noise -- the limit and
the reporting unit are recoverable only by a reader that understands the
format. ``d4/method_note_old.txt`` is a sitting distractor: a real text
file, easy to read, stating a SUPERSEDED limit of 150 C. An agent that
cannot open the .docx and answers from the .txt gives 150, which this
script fails. That is the difference the capability makes, made visible
as two different numbers rather than as a judgment about tool use.

Grading is on the numbers and the unit in the answer document, never on
which tool was called: reaching the same facts another way still passes,
because the property under test is that the document was understood.

Exit 0 = pass.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

# The current note's limit and reporting unit. The distractor says 150.
CORRECT_LIMIT_C = 120
STALE_LIMIT_C = 150


def _numbers(text: str) -> list[float]:
    return [float(m.group()) for m in
            re.finditer(r"\d+(?:\.\d+)?", text)]


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: <script> <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1])
    answer = ws / "d4" / "answer.md"
    if not answer.is_file():
        print("missing d4/answer.md", file=sys.stderr)
        return 1
    text = answer.read_text(encoding="utf-8", errors="replace")
    fails: list[str] = []
    nums = _numbers(text)

    has_correct = any(abs(n - CORRECT_LIMIT_C) < 1e-6 for n in nums)
    # The stale figure is only a defect when it is ASSERTED as the limit.
    # A correct answer very often mentions 150 C in order to dismiss it
    # ("the superseded note's 150 C does not apply"), and failing that
    # would punish the better answer. So the check is on the assertion
    # frame around it, not on the bare number.
    asserts_stale = re.search(
        r"(?i)(do not exceed|limit(?:ed)? (?:to|is)|maximum(?: is)?|not "
        r"exceed|exceed)\s*(?:of\s*)?150", text) is not None

    if not has_correct:
        print(
            f"FAIL: the note's limit of {CORRECT_LIMIT_C} C is not stated. "
            f"The superseded .txt says {STALE_LIMIT_C} C; answering that is "
            "reading the wrong document, and it is the answer an agent that "
            "cannot open the .docx produces.", file=sys.stderr)
        return 1
    if asserts_stale:
        fails.append(
            f"the answer gives {STALE_LIMIT_C} C AS the limit, which is the "
            "SUPERSEDED figure in method_note_old.txt; the current note "
            f"says {CORRECT_LIMIT_C} C")

    # The reporting unit is the second half of the note and is prose, not
    # a number: it proves the whole document was read rather than a
    # single line pattern-matched out of it.
    if not re.search(r"(?i)mg\s*/\s*kg", text):
        fails.append("the reporting unit 'mg/kg' (dry mass) is not stated")
    # mg/L is wrong only when it is the unit being REPORTED IN. The note
    # itself says "not as mg/L", so a good answer repeats that phrase --
    # measured by trying to catch it, which failed a correct answer.
    asserts_mg_per_l = re.search(
        r"(?i)(report(?:ed|s)?)[^.\n]{0,30}\bmg\s*/\s*l\b", text)
    if asserts_mg_per_l and not re.search(
            r"(?i)(not|never|instead of)[^.\n]{0,20}\bmg\s*/\s*l\b", text):
        fails.append("the answer reports results in mg/L, which the note "
                     "explicitly says NOT to do")

    if fails:
        print("FAIL: " + "; ".join(fails), file=sys.stderr)
        return 1
    print(f"OK: limit {CORRECT_LIMIT_C} C from the current note, reporting "
          "unit mg/kg dry mass")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
