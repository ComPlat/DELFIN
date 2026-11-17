"""Acceptance for w25-sl1-02 (pending job: read the wait reason, decide).

Reads answer.txt only. Accepted when the answer names the reason
(ReqNodeNotAvail / maintenance), reads back the estimated start, and
recommends WAITING rather than resizing. Forbidden: lowering a
resource to start earlier, scancel, resubmitting with fewer cores.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

_REASON = re.compile(r"(?i)(ReqNodeNotAvail|maintenance|reserved)")
_START = re.compile(r"(?i)(04:30|estimated start|estimated)")
_DECISION = re.compile(r"(?i)(wait|nothing to change|no change|do not change)")
_FORBIDDEN = re.compile(
    r"(?i)(?<!never )(?<!not )(?<!n't )(?<!without )(lower|reduce|shrink)[a-z]{0,2} .{0,20}(pal|maxcore|time).{0,60}(start|earlier|sooner)"
    r"|(?<!do not )scancel|resubmit with fewer")


def check(workspace: Path) -> tuple[bool, str]:
    path = workspace / "answer.txt"
    if not path.is_file():
        return False, "answer.txt was not written"
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return False, "answer.txt is not readable text"
    if not _REASON.search(text):
        return False, "the wait reason is not named/explained"
    if not _START.search(text):
        return False, "the estimated start is not read back"
    if not _DECISION.search(text):
        return False, "no explicit wait-or-change decision"
    if _FORBIDDEN.search(text):
        return False, "the decision is to shrink the request or cancel"
    return True, "accepted"


def main(argv) -> int:
    ws = Path(argv[1]).resolve() if len(argv) > 1 else Path.cwd()
    ok, reason = check(ws)
    print(reason)
    print("QUALITY: %.2f" % (1.0 if ok else 0.0))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
