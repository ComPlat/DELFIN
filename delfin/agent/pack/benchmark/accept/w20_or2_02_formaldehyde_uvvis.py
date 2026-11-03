"""Acceptance for w20-or2-02 (formaldehyde UV-Vis via the _esd workflow).

Input: the workspace root after the turn. Output: exit 0 + "accepted",
or exit 1 + the reason.

Accepted when ALL hold in or2/formaldehyde/result_uv.txt:

1. An S1 excitation energy in eV inside 4.2 +/- 0.4, with the unit
   stated.
2. An oscillator strength reported (any decimal number next to
   oscillator/fosc/f).
3. The workflow is named: _esd or ESD or TD-DFT (tddft) appears in the
   text — the number must come from DELFIN's excited-state workflow,
   not from memory.
4. The state character (n->pi* or equivalent) is named.

Plain filesystem reads only: no subprocess, clock or network.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

_RESULT = Path("or2") / "formaldehyde" / "result_uv.txt"
_EV = re.compile(r"([-+]?\d+(?:\.\d+)?)\s*eV", re.IGNORECASE)
_FOSC = re.compile(r"(?i)(oscillator|fosc)\D*?(\d+\.\d+)")
_WORKFLOW = re.compile(r"(?i)(_esd|ESD|TD.?DFT|TDDFT)")
_CHARACTER = re.compile(r"(?i)(n.?pi|lone.?pair|character)")


def check(workspace: Path) -> tuple[bool, str]:
    path = workspace / _RESULT
    if not path.is_file():
        return False, f"{_RESULT} was not written"
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return False, f"{_RESULT} is not readable text"
    evs = [float(m.group(1)) for m in _EV.finditer(text)]
    if not any(3.8 <= e <= 4.6 for e in evs):
        return False, f"no S1 energy in 3.8-4.6 eV; found {evs}"
    if not _FOSC.search(text):
        return False, "no oscillator strength reported"
    if not _WORKFLOW.search(text):
        return (False, "workflow not named: the energy must come from "
                       "the _esd / TD-DFT workflow, not from memory")
    if not _CHARACTER.search(text):
        return False, "state character (n->pi*) not named"
    return True, "accepted: S1 in band, fosc, workflow named, character"


def main() -> None:
    ws = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    ok, msg = check(ws)
    print(msg)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
