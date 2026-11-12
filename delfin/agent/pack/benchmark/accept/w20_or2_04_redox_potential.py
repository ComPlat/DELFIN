"""Acceptance for w20-or2-04 (oxidation potential, electrode + cycle).

Input: the workspace root after the turn. Output: exit 0 + "accepted",
or exit 1 + the reason.

Accepted when ALL hold in or2/redox/result.txt:

1. An oxidation potential in V inside 0.55 +/- 0.45 (0.10-1.00).
2. A reference electrode is named next to the number (vs. SHE, Fc/Fc+,
   Ag/Ag+, NHE...). A potential with no electrode is a rejection — that
   is the core of this package's rubric.
3. The thermodynamic cycle / free-energy relation is stated (G_ox,
   G_red, free energy, Delta G).
4. The redox workflow is named (redox_potential / DELFIN redox).

Plain filesystem reads only: no subprocess, clock or network.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

_RESULT = Path("or2") / "redox" / "result.txt"
_POTENTIAL = re.compile(r"(?i)E_?(?:ox|1/2|red)?\D{0,20}?([-+]?\d\.\d+)\s*V\b")
_ELECTRODE = re.compile(r"(?i)(vs\.?\s*(SHE|NHE|Fc|Ag|RHE)|reference "
                        r"electrode|E_ref)")
_CYCLE = re.compile(r"(?i)(thermodynamic cycle|G_?ox|G_?red|Delta.?G|"
                    r"free energy)")
_WORKFLOW = re.compile(r"(?i)(redox_potential|redox workflow|DELFIN'?s? "
                       r"redox)")


def check(workspace: Path) -> tuple[bool, str]:
    path = workspace / _RESULT
    if not path.is_file():
        return False, f"{_RESULT} was not written"
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return False, f"{_RESULT} is not readable text"
    volts = [float(m.group(1)) for m in _POTENTIAL.finditer(text)]
    if not any(0.10 <= v <= 1.00 for v in volts):
        return False, f"no oxidation potential in 0.10-1.00 V; found {volts}"
    if not _ELECTRODE.search(text):
        return (False, "no reference electrode next to the potential — "
                       "an E value without its electrode is not a result")
    if not _CYCLE.search(text):
        return False, "thermodynamic cycle / free-energy relation not stated"
    if not _WORKFLOW.search(text):
        return False, "DELFIN redox workflow not named"
    return True, "accepted: E in band, electrode, cycle, workflow"


def main() -> None:
    ws = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    ok, msg = check(ws)
    print(msg)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
