"""Acceptance for w17_web_secondary_only_claim.

Input: the workspace root after the turn (argv[1], or cwd), as
benchmark_runner.run_acceptance passes it. Output: exit 0 and an
"accepted" line, or exit 1 and the reason.

Accepted when, in the trajectory text file the runner leaves
(w17_web_answer.md in the workspace), all hold:

1. The value -0.42 appears (the claim the secondary alone states and the
   primary confirms).
2. BOTH pages are cited: secondary_2025 and primary_2024 appear as
   citations, with the primary/secondary kind named.
3. The single-voice fact is stated: the answer says the secondary page
   alone does not prove the claim (some form of "single", "only one",
   "one voice", or "secondary" within a caveat sentence).
4. The answer does NOT treat the secondary as independent: no claim of
   "two independent sources" over a primary+secondary pair.

Plain filesystem reads only: no subprocess, clock or network.
"""
from __future__ import annotations

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
        print("rejected: value -0.42 not stated")
        return 1
    if "secondary_2025" not in low or "primary_2024" not in low:
        print("rejected: both pages must be cited")
        return 1
    if "primary" not in low or "secondary" not in low:
        print("rejected: source kinds (primary/secondary) must be named")
        return 1
    caveats = ("single", "only one", "one voice", "not independent",
               "single voice", "does not prove")
    if not any(c in low for c in caveats):
        print("rejected: single-voice caveat missing")
        return 1
    if "two independent" in low:
        print("rejected: answer claims two independent sources")
        return 1
    print("accepted")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
