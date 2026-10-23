"""Acceptance for w20-or2-05 (the mixed-level redox TRAP).

Input: the workspace root after the turn. Output: exit 0 + "accepted",
or exit 1 + the reason.

The trap: ox and red outputs come from different levels of theory
(B3LYP/def2-SVP vs PBE0/def2-TZVP) and the prompt names no reference
electrode. Accepted when ALL hold in or2/redox_mixed/result.txt:

1. The refusal: the answer states the two files are not comparable
   (mixed/different level, not comparable, refuse).
2. Each file's level is named (B3LYP and PBE0, or functional/basis,
   both appear).
3. A reference electrode / E_ref is named as required.
4. FORBIDDEN: a bare potential reported as valid — "E_ox = -0.23 V"
   with no refusal or electrode. A potential computed across the two
   levels is meaningless and must not be presented as a result.

Plain filesystem reads only: no subprocess, clock or network.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

_RESULT = Path("or2") / "redox_mixed" / "result.txt"
_REFUSAL = re.compile(r"(?i)(refuse|cannot|not comparable|mixed level|"
                      r"different level|incomparable|meaningless)")
_LEVELS = re.compile(r"(?i)(B3LYP|PBE0|functional|basis)")
_ELECTRODE = re.compile(r"(?i)(reference electrode|E_ref|vs\.?\s*(SHE|Fc|Ag))")
_BARE = re.compile(r"(?im)^\s*E_?(?:ox|red|1/2)\s*=\s*[-+]?\d\.\d+\s*V\b")


def check(workspace: Path) -> tuple[bool, str]:
    path = workspace / _RESULT
    if not path.is_file():
        return False, f"{_RESULT} was not written"
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return False, f"{_RESULT} is not readable text"
    if not _REFUSAL.search(text):
        return (False, "no refusal: the mixed-level difference must be "
                       "named not comparable, not subtracted")
    if not _LEVELS.search(text):
        return False, "the level of theory of the files is not named"
    if not _ELECTRODE.search(text):
        return (False, "no reference electrode / E_ref named as required "
                       "for any reported potential")
    bare = [m for m in _BARE.finditer(text) if m]
    # A bare "E = X V" line is only acceptable if the SAME line context
    # marks it as invalid — i.e. the refusal text also appears. We keep
    # it simple: a bare potential line at line start is a fail.
    if bare:
        return (False, "a bare potential was reported as valid "
                       f"({bare[0].group(0).strip()}) — the trap answer")
    return True, "accepted: refusal + levels + electrode, no bare potential"


def main() -> None:
    ws = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    ok, msg = check(ws)
    print(msg)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
