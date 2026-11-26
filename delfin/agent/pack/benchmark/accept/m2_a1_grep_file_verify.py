#!/usr/bin/env python3
"""A1 of the capability census, through ``grep_file``.

Question: in ``g1/``, where is ``_settle`` defined and where is it called?

Why the answer differs when the capability works: the name has three
different standing in the tree, and only a real search separates them. It
is DEFINED once (``src/core.py``), CALLED once from ``step`` in that same
file, imported by ``src/run.py`` -- and then SHADOWED there by a local
closure of the same name that has nothing to do with the helper. It is
also mentioned in ``notes/README.md`` as prose. An answer that lists every
file containing the word reports the notes and the closure as call sites;
a search that found only the definition misses the one real caller.

Exit 0 = pass.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: <script> <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1])
    answer = ws / "g1" / "answer.md"
    if not answer.is_file():
        print("missing g1/answer.md", file=sys.stderr)
        return 1
    text = answer.read_text(encoding="utf-8", errors="replace")
    fails: list[str] = []

    # The definition site, named as a file: an answer that does not say
    # WHERE is not an answer to "where is it defined".
    if not re.search(r"(?i)core\.py", text):
        fails.append("the definition site (src/core.py) is not named")

    # The one real call site inside core.py, from step().
    if not re.search(r"(?i)step\s*\(", text):
        fails.append("the caller `step()` is not named, so the one real "
                     "call site in core.py was not found")

    # The shadowing closure must be recognised as NOT a call of the
    # helper. Two accepted shapes: it is excluded, or it is called out as
    # a separate local definition.
    noticed_local = re.search(
        r"(?i)(local|shadow|closure|redefin|separate definition|nested "
        r"def|not the same|unrelated)", text)
    if not noticed_local:
        fails.append(
            "the local `_settle` closure in src/run.py is not distinguished "
            "from the helper in core.py -- reporting it as a call site of "
            "the helper is the wrong answer")

    # The notes file must not be listed as a definition or call site.
    if re.search(r"(?i)notes?/readme", text) and not re.search(
            r"(?i)(prose|mention|not a (definition|call)|documentation|"
            r"README only)", text):
        fails.append("notes/README.md is listed without noting it is only "
                     "prose, not a definition or a call")

    # Count discipline: exactly one definition and one real call.
    if re.search(r"(?i)\b(two|2)\s+(real\s+)?(call|caller|call site)", text):
        fails.append("the answer reports two real call sites; there is one "
                     "(the closure in run.py is a different function)")

    if fails:
        print("FAIL: " + "; ".join(fails), file=sys.stderr)
        return 1
    print("OK: definition and single real call site found, the shadowing "
          "closure distinguished from the helper")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
