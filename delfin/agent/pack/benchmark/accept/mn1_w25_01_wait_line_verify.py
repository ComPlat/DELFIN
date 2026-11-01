"""Accept for w25-mn1-01 (start a background job, wait on a log LINE).

Input: the workspace root after the turn. Output: exit 0 + "accepted",
or exit 1 + the reason. Plain filesystem reads only -- no subprocess,
clock or network.

Accepted when ALL hold:
1. job.log next to the started script holds the PHASE TWO line -- the
   job really ran (started, not merely described).
2. answer.txt states the phase-two reason the log gives (the optimizer
   leaving the first basin).
3. The job was started as a background job: the transcript facts
   (tool-call names) come as transcript_signals.txt, written by the
   benchmark caller from the transcript; missing means
   unmeasured, never accepted.
4. FORBIDDEN: four or more bash_status calls in the transcript record
   (busy polling) -- enforced on the transcript record here, not by a
   wall clock.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

_PHASE_TWO = re.compile(r"(?i)PHASE TWO REACHED")
_REASON = re.compile(r"(?i)(first basin|optimizer|left the)")
_BG = re.compile(r"(?i)bash_background")
# Polling channels, all of them: the status tool, the output reader, and
# any shell re-read of the watched log (cat/tail/head/less). Four or more
# of these in the transcript record is busy polling; three (start, one
# check, one collect) always stay legal.
_POLL = re.compile(r"(?i)\b(bash_status|bash_output)\b")
_LOG_REREAD = re.compile(
    r"(?i)\b(cat|tail|head|less|more)\b[^\n]*job\.log")


def check(workspace: Path) -> tuple[bool, str]:
    folder = workspace / "mn1_w25" / "wait_line"
    log = folder / "job.log"
    if not log.is_file():
        return False, "job.log is missing -- the background job never ran"
    try:
        log_text = log.read_text(encoding="utf-8")
    except OSError as exc:
        return False, f"job.log unreadable: {exc}"
    if not _PHASE_TWO.search(log_text):
        return False, "the PHASE TWO line never reached job.log -- job not collected"
    answer = workspace / "answer.txt"
    if not answer.is_file():
        return False, "answer.txt was not written"
    try:
        text = answer.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return False, "answer.txt is not readable text"
    if not _REASON.search(text):
        return False, ("the phase-two reason (optimizer leaving the first "
                       "basin) is not stated in answer.txt")
    # The transcript facts come from the caller, as a file, so no wall
    # clock and no transcript parsing happen here.
    ts = workspace / "transcript_signals.txt"
    if not ts.is_file():
        return False, ("transcript_signals.txt was not provided -- the "
                       "background-start and polling rubric is unmeasured")
    try:
        signals = ts.read_text(encoding="utf-8")
    except OSError as exc:
        return False, f"transcript_signals.txt unreadable: {exc}"
    if not _BG.search(signals):
        return False, "no bash_background call in the transcript record"
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
