"""Acceptance for w25-sl1-01 (resource choice + justification).

Reads answer.txt + sl1/size/CONTROL_resources.txt only — plain
filesystem reads, no subprocess, clock or network.

Accepted when ALL hold:
1. Both cpu partitions named as one submission (earliest-start rule).
2. PAL/maxcore justified from DELFIN's own sizing, not a guess.
3. The never-resize rule stated.
Forbidden: lowering a resource to start earlier.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

_PARTITIONS = re.compile(r"(?i)(cpu_il|cpu\s*,\s*cpu_il)")
_EARLIEST = re.compile(r"(?i)(earliest|whichever|first (available|fitting)|both partitions)")
_SIZING = re.compile(r"(?i)(advise|auto_configure|cluster_utils|DELFIN's own sizing|calc_policy)")
_RESIZE = re.compile(r"(?i)(never (lower|change|resize)|not (lower|change).{0,40}(start earlier|earlier start)|earliest start.{0,80}never)")
_FORBIDDEN = re.compile(r"(?i)(?<!never )(?<!not )(?<!n't )(?<!without )(lower|reduce|shrink)[a-z]{0,2} (the )?(PAL|maxcore|time limit).{0,60}(start earlier|start sooner|queue sooner)")


def check(workspace: Path) -> tuple[bool, str]:
    path = workspace / "answer.txt"
    if not path.is_file():
        return False, "answer.txt was not written"
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return False, "answer.txt is not readable text"
    if not _PARTITIONS.search(text):
        return False, "the cpu/cpu_il partition pair is not named"
    if not _EARLIEST.search(text):
        return False, "the earliest-start rule is not stated"
    if not _SIZING.search(text):
        return False, "DELFIN's own sizing logic is not cited"
    if not _RESIZE.search(text):
        return False, "the never-resize-to-start-earlier rule is not stated"
    proposal = workspace / "sl1" / "size" / "CONTROL_resources.txt"
    if not proposal.is_file():
        return False, "sl1/size/CONTROL_resources.txt was not written"
    try:
        proposal_text = proposal.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return False, "CONTROL_resources.txt is not readable text"
    if not re.search(r"(?i)\bPAL\b", proposal_text):
        return False, "the CONTROL proposal carries no PAL line"
    if not re.search(r"(?i)maxcore", proposal_text):
        return False, "the CONTROL proposal carries no maxcore line"
    if _FORBIDDEN.search(text):
        return False, "a resource is lowered to start earlier"
    return True, "accepted"


def main(argv) -> int:
    ws = Path(argv[1]).resolve() if len(argv) > 1 else Path.cwd()
    ok, reason = check(ws)
    print(reason)
    print("QUALITY: %.2f" % (1.0 if ok else 0.0))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
