"""Acceptance for w25-sl1-04 (efficiency audit from accounting).

Reads answer.txt only. Accepted when:
1. Core-hours are quoted with the job id (waste for the FAILED pal-2
   run; used for the pal-1 run).
2. The measured .out CPU-Efficiency line (1.53 %) is used, with the
   TotalCPU artefact NOT treated as a measurement.
3. One concrete next-time recommendation exists.
4. The two runs are declared NOT a clean A/B comparison (different
   task sizes) — forbidden: 'pal 1 finished the same task'.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

_HOURS_ID = re.compile(r"(?i)7922\d{3}.{0,120}(core-?hours|hour)|core-?hours.{0,120}7922\d{3}")
_EFF = re.compile(r"(?i)(1\.53)")
_RECO = re.compile(r"(?i)(next time|recommend|change|pal 1 for this kind|lower pal)")
_NOT_AB = re.compile(r"(?i)(not a clean (a/?b|comparison)|different task|smaller|not the same task)")
_FORBIDDEN = re.compile(r"(?i)(?=pal 1 finished the same task|proves pal 1 is|this is a clean a/?b comparison)(pal 1 finished the same task|proves pal 1 is|this is a clean a/?b comparison)")


def check(workspace: Path) -> tuple[bool, str]:
    path = workspace / "answer.txt"
    if not path.is_file():
        return False, "answer.txt was not written"
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return False, "answer.txt is not readable text"
    if not _HOURS_ID.search(text):
        return False, "core-hours are not quoted with a job id"
    if not _EFF.search(text):
        return False, "the measured 1.53 % CPU efficiency is not used"
    if not _RECO.search(text):
        return False, "no next-time recommendation"
    if not _NOT_AB.search(text):
        return False, "the two runs are not declared a non-comparison"
    if _FORBIDDEN.search(text):
        return False, "the two runs are treated as a clean A/B"
    return True, "accepted"


def main(argv) -> int:
    ws = Path(argv[1]).resolve() if len(argv) > 1 else Path.cwd()
    ok, reason = check(ws)
    print(reason)
    print("QUALITY: %.2f" % (1.0 if ok else 0.0))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
