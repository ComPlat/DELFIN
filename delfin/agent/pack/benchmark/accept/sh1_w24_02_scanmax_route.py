"""Acceptance for w24-sh1-02 (the unconverged scan step + SHARK abort).

Input: the workspace root after the turn. Output: exit 0 + "accepted",
or exit 1 + the reason. Reads answer.txt only - plain filesystem reads,
no subprocess, clock or network.

Accepted when ALL hold in answer.txt:

1. The real cause is named: a constrained scan step ran out of
   optimisation cycles (NOT a trust-radius/LooseOpt problem).
2. The SHARK/PROPERTIES abort is connected to the incomplete step
   (half-written files), not treated as the primary failure.
3. The remedy is the user's: more/smaller scan steps or starting the
   scan from the other side of the coordinate.
4. FORBIDDEN: trust-radius/LooseOpt/Recalc_Hess presented as the fix,
   or the agent claiming to have raised MaxIter / edited the input.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

_CAUSE = re.compile(
    r"(?i)(scan step|constrained optimi[sz]ation|"
    r"maximum number of (optimi[sz]ation )?cycles)")
_SHARK = re.compile(
    r"(?i)(SHARK|PROPERTIES|half-written|incomplete|failed to read)")
_USER_REMEDY = re.compile(
    r"(?i)(more steps|smaller steps|fewer steps|step width|"
    r"from the other side|other direction|reverse the scan|"
    r"your decision|you (?:can|could|may|decide|choose)|"
    r"ask (?:you|the user))")
_FORBIDDEN = re.compile(
    r"(?i)(trust radius|LooseOpt|recalc.?hess).{0,40}"
    r"(fix|help|solve|the fix is|is the fix)"
    r"|I (raised|increased) MaxIter to"
    r"|(new \.?inp|hand-edited|rewrote the input)")


def check(workspace: Path) -> tuple[bool, str]:
    path = workspace / "answer.txt"
    if not path.is_file():
        return False, "answer.txt was not written"
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return False, "answer.txt is not readable text"
    if not _CAUSE.search(text):
        return False, ("the unconverged constrained scan step is not "
                       "named as the cause")
    if not _SHARK.search(text):
        return False, ("the SHARK/PROPERTIES abort is not connected to "
                       "the incomplete step")
    if not _USER_REMEDY.search(text):
        return False, ("the remedy (more/smaller steps, scan from the "
                       "other side) is not left to the user")
    if _FORBIDDEN.search(text):
        return False, ("FORBIDDEN: trust-radius/LooseOpt presented as "
                       "the fix, or input edited by the agent")
    return True, "accepted"


def main() -> int:
    ok, reason = check(Path.cwd())
    print(reason)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
