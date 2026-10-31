"""Accept for w26-sh2-01 (metacharacter value passed verbatim).

Input: the workspace root after the turn. Output: exit 0 + "accepted",
or exit 1 + the reason. Plain filesystem reads only -- no subprocess,
clock or network.

Accepted when ALL hold:
1. quoted.txt exists in sh2_w26/metachar/.
2. Its content equals the fixture token EXACTLY:
   rate=1e-6*2+5^q&x>y$z
   Any shell interpretation corrupts it -- an unquoted $z expands to
   empty, `&` or `>` background or redirect, `^` is lost to escaping,
   `*` globs -- so only a correctly quoted single-quoted (or
   otherwise literally read) value passes byte for byte.
3. answer.txt states the string (the agent reports what it wrote).
"""
from __future__ import annotations

import sys
from pathlib import Path

TOKEN = "rate=1e-6*2+5^q&x>y$z"


def check(workspace: Path) -> tuple[bool, str]:
    quoted = workspace / "sh2_w26" / "metachar" / "quoted.txt"
    if not quoted.is_file():
        return False, "sh2_w26/metachar/quoted.txt was not written"
    try:
        got = quoted.read_text(encoding="utf-8").strip()
    except (OSError, UnicodeDecodeError) as exc:
        return False, f"quoted.txt unreadable: {exc}"
    if got != TOKEN:
        return False, ("quoted.txt does not hold the token verbatim "
                       f"(got {got!r}, expected {TOKEN!r}) -- the shell "
                       "interpreted part of it")
    answer = workspace / "answer.txt"
    if not answer.is_file():
        return False, "answer.txt was not written"
    try:
        text = answer.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return False, "answer.txt is not readable text"
    if TOKEN not in text:
        return False, "the reported string in answer.txt is not the token"
    return True, "accepted"


def main() -> int:
    ws = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    ok, reason = check(ws)
    print(reason)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
