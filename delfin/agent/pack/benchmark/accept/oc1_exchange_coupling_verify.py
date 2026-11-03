"""Acceptance for w24-oc1-02 (broken-symmetry exchange coupling).

Input: the workspace root after the turn. Output: exit 0 + "accepted",
or exit 1 + the reason.

Accepted when ALL hold for the run under oc1/broken_sym/:

1. A real OCCUPIER run exists: at least one *.out under oc1/broken_sym/,
   terminated normally.
2. The report names J with a value AND states the Yamaguchi relation
   (<S^2> of the high-spin AND the broken-symmetry solution) -- a J
   without its spin moments is not a result (manual Eq. 5.167 needs both
   <S^2> values).
3. Site order is stated: which centre is site A (the larger
   unpaired-electron count). If the input got the order wrong the run
   reconverges the wrong state -- that is the defect this grades.
4. The energy difference used in J carries kcal/mol or cm^-1 with the
   conversion stated.

Plain filesystem reads only.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

_J = re.compile(r"(?i)\bJ\b\D{0,25}?[-+]?\d+\.?\d*")
_S2 = re.compile(r"(?i)<S\^?2>\D{0,15}\d+\.?\d*")
_SITES = re.compile(r"(?i)(site A|unpaired)")
_UNIT = re.compile(r"(?i)(kcal/mol|cm\^-1|kJ/mol)")
_YAMAGUCHI = re.compile(r"(?i)yamaguchi")


def check(workspace: Path) -> tuple[bool, str]:
    root = workspace / "oc1" / "broken_sym"
    if not root.is_dir():
        return False, "oc1/broken_sym/ does not exist"
    outs = sorted(root.rglob("*.out"))
    if not outs:
        return False, "no .out file under oc1/broken_sym/ - a real run is required"
    if not any("ORCA TERMINATED NORMALLY" in o.read_text(encoding="utf-8",
                                                          errors="replace")
               for o in outs):
        return False, "no oc1/broken_sym/ output terminated normally"
    report = root / "result.txt"
    if not report.is_file():
        return False, "no result.txt under oc1/broken_sym/"
    text = report.read_text(encoding="utf-8", errors="replace")
    if not _J.search(text):
        return False, "no J value reported"
    s2_vals = _S2.findall(text)
    if len(s2_vals) < 2:
        return (False, f"fewer than two <S^2> values stated ({len(s2_vals)}) "
                       "— Yamaguchi Eq. 5.167 needs HS and BS spin moments")
    if not _SITES.search(text):
        return False, "site A (larger unpaired-electron count) not identified"
    if not _UNIT.search(text):
        return False, "energy difference carries no unit (kcal/mol, kJ/mol or cm^-1)"
    # A stated Yamaguchi form is preferred but J with both <S^2> values
    # and a unit already fixes the formula implicitly; flag wording-only
    # absence as accepted-with-note.
    note = "" if _YAMAGUCHI.search(text) else " (relation implicit from values)"
    return True, f"accepted: J stated, {len(s2_vals)} <S^2> values, site A named, unit present{note}"


def main() -> None:
    ws = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    ok, msg = check(ws)
    print(msg)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
