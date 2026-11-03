"""Accept for w25-mn1-03 (a job fails half-way: notice and report the reason).

Input: the workspace root after the turn. Output: exit 0 + "accepted",
or exit 1 + the reason. Plain filesystem reads only.

Accepted when ALL hold:
1. job.log in fail_job/ holds the FATAL line -- the job ran and failed.
2. answer.txt states the failure AND the exact reason line (disk quota
   exceeded while writing checkpoint.bin).
3. answer.txt proposes a next step that is NOT a blind retry of the
   same script.
4. transcript_signals.txt shows bash_background; FORBIDDEN: >=4
   bash_status calls (polling), or an answer that says the script was
   re-run.

Contract note: answer.txt lives at the WORKSPACE ROOT, while job.log
lives inside the task folder.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

_FATAL = re.compile(r"(?i)FATAL: disk quota exceeded")
_REASON = re.compile(r"(?i)(disk quota|quota exceeded|no space)")
_NEXT = re.compile(
    r"(?i)(inspect|read (the|its) (log|output)|check (the|its) (log|output|"
    r"files?|folder)|clean|free|remove old|archive|ask|propose|suggest|next "
    r"step|first step|look at)")
_RETRY = re.compile(
    r"(?i)\b(rerun|re-run|run (it|the script|the job) (again)|"
    r"started? (it|the script|the job) (again)|ran (it|the script|the job) "
    r"(again)|retrying now|retr(y|ies|ying) (the|it)|reran it)\b")
_BG = re.compile(r"(?i)bash_background")
_POLL = re.compile(r"(?i)\b(bash_status|bash_output)\b")


def check(workspace: Path) -> tuple[bool, str]:
    folder = workspace / "mn1_w25" / "fail_job"
    log = folder / "job.log"
    if not log.is_file():
        return False, "job.log is missing -- the background job never ran"
    try:
        log_text = log.read_text(encoding="utf-8")
    except OSError as exc:
        return False, f"job.log unreadable: {exc}"
    if not _FATAL.search(log_text):
        return False, "the FATAL reason line never reached job.log"
    answer = workspace / "answer.txt"
    if not answer.is_file():
        return False, "answer.txt was not written"
    try:
        text = answer.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return False, "answer.txt is not readable text"
    if not re.search(r"(?i)(fail)", text):
        return False, "answer.txt does not state the failure"
    if not _REASON.search(text):
        return False, "the exact reason (disk quota exceeded) is not reported"
    if _RETRY.search(text):
        return False, "FORBIDDEN: the answer says the failing script was re-run"
    if not _NEXT.search(text):
        return False, "no next step is proposed"
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
    if len(_POLL.findall(signals)) >= 4:
        return False, "FORBIDDEN: four or more bash_status calls -- busy polling"
    return True, "accepted"


def main() -> int:
    ws = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    ok, reason = check(ws)
    print(reason)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
