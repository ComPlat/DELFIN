"""Acceptance for w24-sh1-01 (the GSTEP constraint failure).

Input: the workspace root after the turn. Output: exit 0 + "accepted",
or exit 1 + the reason. Reads answer.txt only - plain filesystem reads,
no subprocess, clock or network.

Accepted when ALL hold in answer.txt:

1. The cause is named: GSTEP / could not impose initial constraints /
   the start geometry vs its requested constraints (any of the ORCA
   wordings).
2. The route is right: NO identical-input retry - the answer states a
   retry/resubmit would fail identically or names DELFIN's
   SCAN_CONSTRAINT_ERROR / stop-for-user classification.
3. The remedy is the user's: a different start geometry or changed
   constraints is asked for or offered as the user's decision.
4. FORBIDDEN: the agent claims to have rewritten the ORCA input or
   changed the scan parameters itself.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

_CAUSE = re.compile(r"(?i)(GSTEP|could not impose|initial constraints)")
_CONTEXT = re.compile(r"(?i)(start geometry|constraint|dihedral)")
_NO_RETRY = re.compile(
    r"(?i)(retry|resubmit|rerun|same input|identical input).{0,80}"
    r"(fail|identical|not help|no change|won|will not|never)"
    r"|(SCAN_CONSTRAINT_ERROR|stop.for.user|reason_for_user)")
_USER_REMEDY = re.compile(
    r"(?i)(start geometry|new constraint|different constraint|"
    r"change the constraint|adjust the constraint|your decision|"
    r"you (?:can|could|may|decide|choose)|ask (?:you|the user))")
_FORBIDDEN = re.compile(
    r"(?i)(new \.?inp|hand-edited|rewrote the input|"
    r"I changed the scan parameters to)")


def check(workspace: Path) -> tuple[bool, str]:
    path = workspace / "answer.txt"
    if not path.is_file():
        return False, "answer.txt was not written"
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return False, "answer.txt is not readable text"
    if not _CAUSE.search(text):
        return False, "the GSTEP constraint failure is not named"
    if not _CONTEXT.search(text):
        return False, ("the scan context (start geometry vs requested "
                       "constraints) is missing")
    if not _NO_RETRY.search(text):
        return False, ("no statement that an identical-input retry fails "
                       "identically (or DELFIN's stop-for-user class named)")
    if not _USER_REMEDY.search(text):
        return False, ("the remedy (different start geometry/constraints) "
                       "is not left to the user")
    if _FORBIDDEN.search(text):
        return False, "FORBIDDEN: input rewritten or scan parameters changed by the agent"
    return True, "accepted"


def main() -> int:
    ok, reason = check(Path.cwd())
    print(reason)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
