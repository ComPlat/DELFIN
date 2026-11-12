"""Acceptance for w20-or2-01 (formaldehyde IR, scaled frequencies).

Input: the workspace root after the turn. Output: exit 0 + "accepted",
or exit 1 + the reason.

Accepted when ALL hold in or2/formaldehyde/result.txt:

1. A C=O stretch is reported in cm-1 inside the band 1725 +/- 60
   (scaled value; the raw band is 1665-1785 as stated in the task).
2. A scaling factor is named next to a source (any of: manual, paper,
   citation, DOI, ORCA manual, IUPAC, database keyword).
3. The point group C2v is stated.

Plain filesystem reads only: no subprocess, clock or network.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

_RESULT = Path("or2") / "formaldehyde" / "result.txt"
_BAND = (1725.0 - 60.0, 1725.0 + 60.0)
_SOURCE = re.compile(r"(?i)(source|cited|reference|doi|manual|paper|"
                     r"literature|database|iupac)")
_SCALE = re.compile(r"(?i)scal\w*\s*(factor)?\s*[=:≈]?\s*0?\.\d{2,4}|"
                    r"0\.\d{3}\b")
_C2V = re.compile(r"(?i)C2v")
_NUMBER = re.compile(r"([-+]?\d+(?:\.\d+)?)\s*cm-?1", re.IGNORECASE)


def check(workspace: Path) -> tuple[bool, str]:
    path = workspace / _RESULT
    if not path.is_file():
        return False, f"{_RESULT} was not written"
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return False, f"{_RESULT} is not readable text"
    freqs = [float(m.group(1)) for m in _NUMBER.finditer(text)]
    if not any(_BAND[0] <= n <= _BAND[1] for n in freqs):
        return (False, f"no C=O stretch in {_BAND[0]}-{_BAND[1]} cm-1; "
                       f"found {freqs}")
    if not _SCALE.search(text):
        return False, "no scaling factor stated"
    if not _SOURCE.search(text):
        return (False, "scaling factor without a source — the factor "
                       "must be stated with where it comes from")
    if not _C2V.search(text):
        return False, "point group C2v not stated"
    return True, "accepted: scaled C=O in band, factor + source, C2v"


def main() -> None:
    ws = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    ok, msg = check(ws)
    print(msg)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
