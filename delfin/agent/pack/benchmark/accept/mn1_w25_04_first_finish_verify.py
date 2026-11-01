"""Accept for w25-mn1-04 (two background jobs: act on whichever finishes first).

Input: the workspace root after the turn. Output: exit 0 + "accepted",
or exit 1 + the reason. Plain filesystem reads only.

Accepted when ALL hold:
1. Both output logs exist (fast.log AND slow.log): both jobs were
   started and ran.
2. answer.txt names the FAST job as the one that finished first and
   reports its output (3 iterations, residual 1.2e-6).
3. answer.txt states the other job's status at that moment (still
   running / not finished).
4. transcript_signals.txt shows TWO bash_background calls (both jobs
   started in the background); FORBIDDEN: >=4 bash_status calls.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

_FAST_OUT = re.compile(r"(?i)(3 iterations|1\.2e-?6|residual)")
_FIRST = re.compile(r"(?i)(fast|first|fast_job|fast\.log)")
_OTHER = re.compile(r"(?i)(still running|not finished|not yet|slower|"
                    r"other job|slow job .*(running|unfinished)|running)")
_BG = re.compile(r"(?i)bash_background")
# Polling channels: status + output tools, and repeated re-reads of the
# watched logs. >=4 status calls, or >=4 log rereads, is polling.
_POLL = re.compile(r"(?i)\b(bash_status|bash_output)\b")
_LOG_REREAD = re.compile(
    r"(?i)\b(cat|tail|head|less|more)\b[^\n]*(fast\.log|slow\.log|two_jobs)")


def check(workspace: Path) -> tuple[bool, str]:
    folder = workspace / "mn1_w25" / "two_jobs"
    fast = folder / "fast.log"
    slow = folder / "slow.log"
    if not fast.is_file():
        return False, "fast.log is missing -- the fast job was not started/run"
    if not slow.is_file():
        return False, ("slow.log is missing -- the slow job was not started; "
                       "only acting on the first finisher does not drop the "
                       "second one")
    try:
        fast_text = fast.read_text(encoding="utf-8")
    except OSError as exc:
        return False, f"fast.log unreadable: {exc}"
    if not _FAST_OUT.search(fast_text):
        return False, "fast.log does not hold the fast job's output"
    answer = workspace / "answer.txt"
    if not answer.is_file():
        return False, "answer.txt was not written"
    try:
        text = answer.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return False, "answer.txt is not readable text"
    if not _FIRST.search(text):
        return False, "answer.txt does not name the job that finished first"
    if not _FAST_OUT.search(text):
        return False, "the fast job's output (3 iterations, residual) is not summarised"
    if not _OTHER.search(text):
        return False, ("the other job's status at that moment (still "
                       "running) is not stated")
    ts = workspace / "transcript_signals.txt"
    if not ts.is_file():
        return False, ("transcript_signals.txt was not provided -- the "
                       "background-start and polling rubric is unmeasured")
    try:
        signals = ts.read_text(encoding="utf-8")
    except OSError as exc:
        return False, f"transcript_signals.txt unreadable: {exc}"
    if len(_BG.findall(signals)) < 2:
        return False, ("fewer than two bash_background calls -- not both "
                       "jobs were started in the background")
    polls = len(_POLL.findall(signals)) + len(_LOG_REREAD.findall(signals))
    if polls >= 4:
        return False, ("FORBIDDEN: %d status/output/log-reread calls -- busy "
                       "polling" % polls)
    return True, "accepted"


def main() -> int:
    ws = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    ok, reason = check(ws)
    print(reason)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
