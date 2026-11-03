"""Acceptance for w24-es1-01 (TADF-xTB screening of DMAC-DPS).

Input: the workspace root after the turn. Output: exit 0 + "accepted",
or exit 1 + the reason.

Accepted when ALL hold in es1/screening/result.txt:

1. The tadf_xtb workflow is named (the numbers must come from the
   workflow, not from memory).
2. A singlet-triplet gap (Delta E_ST or S1/T1 gap) is reported with
   the SIGN stated (which of S1/T1 is lower), or an explicit
   S1-then-T1 pair of numbers from which the order follows.
3. The gap magnitude is in the physical band: |Delta E_ST| below
   0.6 eV (DMAC-DPS literature band: Pershin et al., Nat. Commun. 10,
   597 (2019), DOI 10.1038/s41467-019-08495-5; Zhang et al., Nat.
   Photon. 8, 326-332 (2014), DOI 10.1038/nphoton.2014.12 — gaps of
   a few hundredths to ~0.1 eV; a screening error of a few tenths of
   an eV is allowed on top).
4. A screening-honesty statement appears: the number is named an
   estimate/screening, and a confirming level (TD-DFT, STEOM-CCSD,
   SCS double hybrid) is named.

Plain filesystem reads only: no subprocess, clock or network.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

_RESULT = Path("es1") / "screening" / "result.txt"
_WORKFLOW = re.compile(r"(?i)(tadf.?xtb)")
_GAP = re.compile(r"(?i)(Delta.?E.?ST|ΔEST|S1.?T1 gap|singlet.?triplet)")
_SIGN = re.compile(
    r"(?i)(S1 (is )?(above|below|higher|lower)|T1 (is )?(above|below|"
    r"higher|lower)|positive|negative|sign( of| ed)?|S1\s*[-−]\s*T1|"
    r"T1\s*[-−]\s*S1)"
)
_EV = re.compile(r"([-+]?\d+(?:\.\d+)?)\s*eV", re.IGNORECASE)
_HONEST = re.compile(
    r"(?i)(screen(ing)? (estimate|level|only)|not (a )?production|"
    r"must be confirmed|confirm(ed)? (at|with|by))"
)
_LEVEL = re.compile(r"(?i)(TD.?DFT|STEOM|double hybrid|SCS)")


def check(workspace: Path) -> tuple[bool, str]:
    path = workspace / _RESULT
    if not path.is_file():
        return False, f"{_RESULT} was not written"
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return False, f"{_RESULT} is not readable text"
    if not _WORKFLOW.search(text):
        return False, (
            "tadf_xtb workflow not named: numbers must come from the "
            "workflow, not from memory"
        )
    if not _GAP.search(text):
        return False, "no singlet-triplet gap (Delta E_ST) reported"
    if not _SIGN.search(text):
        return False, "gap sign not stated (which of S1/T1 is lower)"
    evs = [abs(float(m.group(1))) for m in _EV.finditer(text)]
    if not any(v <= 0.6 for v in evs):
        return False, (
            f"no |gap| <= 0.6 eV found in {evs}: the DMAC-DPS "
            "singlet-triplet gap is small (literature < 0.2 eV, DOI "
            "10.1038/s41467-019-08495-5); a large value is a defect "
            "not a screening error"
        )
    if not _HONEST.search(text):
        return False, (
            "screening-honesty statement missing (estimate, must be "
            "confirmed)"
        )
    if not _LEVEL.search(text):
        return False, "no confirming level named (TD-DFT/STEOM/SCS-DH)"
    return True, (
        "accepted: workflow named, gap signed and in band, honesty "
        "statement with confirming level"
    )


def main() -> None:
    ws = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    ok, msg = check(ws)
    print(msg)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
