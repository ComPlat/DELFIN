"""Accept for w25-mn1-02 (wait until a file appears, then summarise it).

Input: the workspace root after the turn. Output: exit 0 + "accepted",
or exit 1 + the reason. Plain filesystem reads only.

Accepted when ALL hold:
1. report.json exists with n=5 and mean_kj_per_mol=41.73 -- the job ran.
2. answer.txt reports both numbers WITH UNITS (kJ/mol for the mean).
3. transcript_signals.txt (from the caller) shows bash_background.
4. FORBIDDEN: >=4 bash_status calls (busy polling), or >=3 identical
   file-existence probes (ls / test -f on report.json).
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

_N, _MEAN = 5, 41.73
_UNIT = re.compile(r"(?i)kJ\s*/?\s*mol")
_BG = re.compile(r"(?i)bash_background")
# Every polling channel: the status tool, the output reader, and repeated
# shell probes of the file or its folder (ls/stat/test -f/cat/tail).
# >=4 status calls or >=3 probes is polling, not waiting.
_POLL = re.compile(r"(?i)\b(bash_status|bash_output)\b")
_PROBE = re.compile(
    r"(?i)\b(ls|test -f|stat|cat|tail|head)\b[^\n]*"
    r"(report\.json|appear_file)")


def check(workspace: Path) -> tuple[bool, str]:
    folder = workspace / "mn1_w25" / "appear_file"
    report = folder / "report.json"
    if not report.is_file():
        return False, "report.json is missing -- the background job never ran"
    try:
        data = json.loads(report.read_text(encoding="utf-8"))
    except Exception as exc:  # noqa: BLE001
        return False, f"report.json does not parse: {exc}"
    if not isinstance(data, dict):
        return False, "report.json is not a JSON object"
    if int(data.get("n") or 0) != _N:
        return False, f"report.json n={data.get('n')!r}, expected {_N}"
    mean = data.get("mean_kj_per_mol")
    if not isinstance(mean, (int, float)) or abs(float(mean) - _MEAN) > 5e-3:
        return False, f"report.json mean={mean!r}, expected {_MEAN} (+-5e-3)"
    answer = workspace / "answer.txt"
    if not answer.is_file():
        return False, "answer.txt was not written"
    try:
        text = answer.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return False, "answer.txt is not readable text"
    if str(_N) not in text:
        return False, "answer.txt does not report n"
    if "41.7" not in text:
        return False, "answer.txt does not report the mean value"
    if not _UNIT.search(text):
        return False, "the mean is reported without its unit (kJ/mol)"
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
        return False, "FORBIDDEN: four or more status/output calls -- busy polling"
    if len(_PROBE.findall(signals)) >= 3:
        return False, ("FORBIDDEN: three or more file/folder existence "
                       "probes -- polling the file in a loop")
    return True, "accepted"


def main() -> int:
    ws = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    ok, reason = check(ws)
    print(reason)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
