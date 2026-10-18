"""Acceptance for w17_web_outdated_page.

Input: the workspace root after the turn (argv[1]). Output: exit 0 +
"accepted", or exit 1 + reason.

Accepted when w17_web_answer.md:

1. States the current value -0.42 (from the 2024 primary) as the value to
   use, and does NOT present 0.55 (the 2020 page's value) as current.
2. Cites both pages: outdated_2020 and primary_2024.
3. Contains a retrieved date (YYYY-MM-DD) and references the page's own
   revision date ("revised", "last revised", "2020").
4. Explicitly reasons about recency: some form of "older", "outdated",
   "superseded", "revised", "current".

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
    if "-0.42" not in text:
        print("rejected: the current value -0.42 must be stated")
        return 1
    if "outdated_2020" not in low or "primary_2024" not in low:
        print("rejected: both pages must be cited")
        return 1
    if not re.search(r"20\d\d-\d\d-\d\d", text):
        print("rejected: no retrieved date (YYYY-MM-DD) present")
        return 1
    if "revised" not in low and "2020" not in low:
        print("rejected: the page's own revision date must be referenced")
        return 1
    if not any(w in low for w in ("older", "outdated", "superseded",
                                  "revised", "current", "newer")):
        print("rejected: recency reasoning missing")
        return 1
    if "0.55" in text and not any(w in low for w in
                                  ("older", "outdated", "superseded",
                                   "not current", "revised")):
        # 0.55 may appear only when explicitly marked as the old value.
        print("rejected: 0.55 presented without being marked outdated")
        return 1
    print("accepted")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
