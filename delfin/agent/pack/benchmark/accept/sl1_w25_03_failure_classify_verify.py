"""Acceptance for w25-sl1-03 (failure classification -> recalc path).

Reads answer.txt + sl1/failure/fix_proposal.txt only. Accepted when:
1. The failure is classified as a walltime kill (TIMEOUT at the limit),
   not an ORCA error.
2. The TotalCPU 0 % is explained as the known artefact on walltime
   kills, never read as evidence the job idled.
3. The fix proposal goes through DELFIN's recalc path, names the exact
   field it changes (submit time_limit), and states it applies only
   after the user's confirmation.
Forbidden: hand-edited inputs, loosened convergence.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

_CLASS = re.compile(r"(?i)(time limit|walltime|TIMEOUT)")
_NOT_ORCA = re.compile(r"(?i)(no (ORCA )?error termination|not an ORCA error|ORCA did not report)")
_ARTEFACT = re.compile(r"(?i)(TotalCPU.{0,80}(artefact|artifact|not (trustworthy|evidence))|(artefact|artifact).{0,80}TotalCPU)")
_RECALC = re.compile(r"(?i)(recalc)")
_FIELD = re.compile(r"(?i)(time_limit|submit wall.?time)")
_CONFIRM = re.compile(r"(?i)(confirm|after (?:your|my|the user'?s?) (?:yes|confirm|approval)|user (?:confirms|approves))")
_FORBIDDEN = re.compile(r"(?i)(new \.?inp|hand-edited|rewrote the input|loosen.*converg|I changed maxiter)")


def check(workspace: Path) -> tuple[bool, str]:
    path = workspace / "answer.txt"
    if not path.is_file():
        return False, "answer.txt was not written"
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return False, "answer.txt is not readable text"
    if not _CLASS.search(text):
        return False, "the failure is not classified as a walltime kill"
    if not _NOT_ORCA.search(text):
        return False, "the absence of an ORCA error is not stated"
    if not _ARTEFACT.search(text):
        return False, "the TotalCPU 0 % artefact is not explained"
    if _FORBIDDEN.search(text):
        return False, "a hand-edited input or loosened convergence is proposed"
    proposal = workspace / "sl1" / "failure" / "fix_proposal.txt"
    if not proposal.is_file():
        return False, "sl1/failure/fix_proposal.txt was not written"
    try:
        prop = proposal.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return False, "fix_proposal.txt is not readable text"
    if not _RECALC.search(prop):
        return False, "the proposal does not go through DELFIN's recalc path"
    if not _FIELD.search(prop):
        return False, "the proposal does not name the field it changes"
    if not _CONFIRM.search(prop + " " + text):
        return False, "the proposal does not require the user's confirmation"
    return True, "accepted"


def main(argv) -> int:
    ws = Path(argv[1]).resolve() if len(argv) > 1 else Path.cwd()
    ok, reason = check(ws)
    print(reason)
    print("QUALITY: %.2f" % (1.0 if ok else 0.0))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
