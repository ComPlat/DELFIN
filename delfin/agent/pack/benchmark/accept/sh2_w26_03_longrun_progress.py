"""Accept for w26-sh2-03 (long command: progress and the result).

Input: the workspace root after the turn. Output: exit 0 + "accepted",
or exit 1 + the reason. Plain filesystem reads only -- no subprocess,
clock or network.

The transcript facts come as transcript_signals.txt, written by the
benchmark caller from the transcript (same convention as the MN1
family); missing means unmeasured, never accepted.

Accepted when ALL hold:
1. job.log holds the DONE line with the fixed numbers -- the job
   really ran to the end.
2. answer.txt reports the three numbers WITH their units: 84 frames,
   rmsd 0.183 A, energy -2141.77 kJ/mol.
3. The transcript record shows a background start (bash_background)
   and a collect; busy polling -- four or more status/output/log-reread
   calls -- is FORBIDDEN (same rule the MN1 family applies).
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

_DONE = re.compile(r"(?i)DONE: 84 frames, rmsd 0\.183 A, "
                   r"energy -2141\.77 kJ/mol")
_FRAMES = re.compile(r"(?i)84 frames")
_RMSD = re.compile(r"(?i)rmsd\D{0,4}0\.183")
_ENERGY = re.compile(r"(?i)energy\D{0,8}-2141\.77\s*kJ/mol")
_BG = re.compile(r"(?i)\b(bash_background|run_in_background)\b")
_POLL = re.compile(r"(?i)\b(bash_status|bash_output)\b")
_LOG_REREAD = re.compile(
    r"(?i)\b(cat|tail|head|less|more)\b[^\n]*job\.log")


def check(workspace: Path) -> tuple[bool, str]:
    log = workspace / "sh2_w26" / "longrun" / "job.log"
    if not log.is_file():
        return False, "job.log is missing -- the job never ran"
    try:
        log_text = log.read_text(encoding="utf-8")
    except OSError as exc:
        return False, f"job.log unreadable: {exc}"
    if not _DONE.search(log_text):
        return False, ("the DONE line never reached job.log -- the job "
                       "was not collected to the end")
    answer = workspace / "answer.txt"
    if not answer.is_file():
        return False, "answer.txt was not written"
    try:
        text = answer.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return False, "answer.txt is not readable text"
    for name, pat in (("frames", _FRAMES), ("rmsd", _RMSD),
                      ("energy", _ENERGY)):
        if not pat.search(text):
            return False, f"the {name} value (with its unit) is not reported"
    ts = workspace / "transcript_signals.txt"
    if not ts.is_file():
        return False, ("transcript_signals.txt was not provided -- the "
                       "background-start and polling rubric is unmeasured")
    try:
        signals = ts.read_text(encoding="utf-8")
    except OSError as exc:
        return False, f"transcript_signals.txt unreadable: {exc}"
    if not _BG.search(signals):
        return False, "no background start in the transcript record"
    polls = len(_POLL.findall(signals)) + len(_LOG_REREAD.findall(signals))
    if polls >= 4:
        return False, ("FORBIDDEN: %d status/output/log-reread calls -- busy "
                       "polling instead of waiting on the condition" % polls)
    return True, "accepted"


def main() -> int:
    ws = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    ok, reason = check(ws)
    print(reason)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
