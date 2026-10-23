"""Acceptance for w20-or2-03 (spin states via OCCUPIER).

Input: the workspace root after the turn. Output: exit 0 + "accepted",
or exit 1 + the reason.

Accepted when ALL hold in or2/spin/result.txt:

1. The OCCUPIER workflow is named (occupier / OCCUPIER) — the spin
   comparison must come from DELFIN's preferred-spin workflow.
2. More than one multiplicity is reported (the comparison actually
   happened), and 5 is among them or named the ground state.
3. Relative energies in kcal/mol are stated (a bare multiplicities
   list is not a comparison).
4. An S2 / spin-contamination value is reported for the chosen state.

Plain filesystem reads only: no subprocess, clock or network.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

_RESULT = Path("or2") / "spin" / "result.txt"
_WORKFLOW = re.compile(r"(?i)occupier")
_MULTI = re.compile(r"(?i)(multiplicity|mult\.?|2S\+1)\D{0,12}(\d)")
_KCAL = re.compile(r"kcal/mol", re.IGNORECASE)
_S2 = re.compile(r"(?i)(S\^?2|spin contamination)\D{0,15}\d")


def check(workspace: Path) -> tuple[bool, str]:
    path = workspace / _RESULT
    if not path.is_file():
        return False, f"{_RESULT} was not written"
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return False, f"{_RESULT} is not readable text"
    if not _WORKFLOW.search(text):
        return (False, "OCCUPIER workflow not named — the comparison "
                       "must come from DELFIN's preferred-spin workflow")
    multis = {int(m.group(2)) for m in _MULTI.finditer(text)}
    if len(multis) < 2:
        return False, f"fewer than two multiplicities compared: {multis}"
    if 5 not in multis and not re.search(r"(?i)ground\D{0,10}5", text):
        return False, "multiplicity 5 not among the compared/named states"
    if not _KCAL.search(text):
        return False, "relative energies not given in kcal/mol"
    if not _S2.search(text):
        return False, "no S2 / spin-contamination value reported"
    return True, "accepted: OCCUPIER named, >=2 multiplicities, kcal/mol, S2"


def main() -> None:
    ws = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    ok, msg = check(ws)
    print(msg)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
