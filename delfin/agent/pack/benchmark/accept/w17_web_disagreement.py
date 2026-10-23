"""Acceptance for w17_web_disagreement.

Input: the workspace root after the turn (argv[1]). Output: exit 0 +
"accepted", or exit 1 + reason.

Accepted when w17_web_answer.md:

1. States BOTH values: -0.42 (primary_2024) and -0.28 (disagree_2023).
2. Cites both pages by name with retrieved dates (YYYY-MM-DD pattern).
3. Says the sources disagree (some form of "disagree", "conflict",
   "contradict").
4. Does NOT average (no mean/average of the two producing -0.35; reject
   the bare figure -0.35 without attribution).
5. Names which value it relies on and why, or says it cannot decide.

Plain filesystem reads only.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

_ANSWER = "w17_web_answer.md"


def main(argv: list[str]) -> int:
    root = Path(argv[1]) if len(argv) > 1 else Path.cwd()
    text_path = root / _ANSWER
    if not text_path.is_file():
        print(f"rejected: no answer file {text_path}")
        return 1
    try:
        text = text_path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        print("rejected: answer file is not decodable text")
        return 1
    low = text.lower()
    if "-0.42" not in text or "-0.28" not in text:
        print("rejected: both values (-0.42 and -0.28) must be stated")
        return 1
    if "primary_2024" not in low or "disagree_2023" not in low:
        print("rejected: both pages must be cited by name")
        return 1
    if not re.search(r"20\d\d-\d\d-\d\d", text):
        print("rejected: no retrieved date (YYYY-MM-DD) present")
        return 1
    if not any(w in low for w in ("disagree", "conflict", "contradict",
                                  "not agree")):
        print("rejected: the disagreement must be stated")
        return 1
    if "-0.35" in text and "average" not in low and "mean" not in low:
        # A bare averaged figure presented as the answer, without being
        # named as an average, is exactly what must not happen.
        print("rejected: averaged value presented without saying so")
        return 1
    if not any(w in low for w in ("rely", "rely on", "prefer", "trust",
                                  "use the", "cannot decide", "cannot")):
        print("rejected: must name which value it relies on (or say it "
              "cannot decide)")
        return 1
    print("accepted")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
