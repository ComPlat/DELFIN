"""Acceptance for w24-es1-02 (formaldehyde ISC rate via the esd workflow).

Input: the workspace root after the turn. Output: exit 0 + "accepted",
or exit 1 + the reason.

Accepted when ALL hold in es1/isc/result.txt:

1. The esd / ISC path is named (the rate must come from DELFIN's
   excited-state workflow, not from memory).
2. A rate constant with a unit (s^-1) is reported.
3. The order of magnitude is in the physical band: 1e6 - 1e9 s^-1
   (formaldehyde S1->T1 ISC at TD-DFT level; DELFIN's own esd_numerics
   converged reference measured on ORCA 6.1.1 is 2.46e7 s^-1 with a
   fully decayed correlation window; band is wide to catch a wrong
   sign, an unscaled unit, or an unconverged window reported as
   converged).
4. The correlation-window convergence question is addressed: either
   the run is stated converged, or a NOT-A-RESULT warning (negative
   rate, truncated window) is reported rather than hidden.

Plain filesystem reads only: no subprocess, clock or network.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

_RESULT = Path("es1") / "isc" / "result.txt"
_WORKFLOW = re.compile(r"(?i)(esd|ISC|intersystem)")
_RATE = re.compile(
    r"([+-]?\d+(?:\.\d+)?\s*[eE]\s*[+-]?\d+)\s*s\^?-?1"
)
_CONVERGENCE = re.compile(
    r"(?i)(NOT.?A.?RESULT|correlation window|converged|MAXTIME|"
    r"truncat)"
)


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
            "esd/ISC path not named: the rate must come from "
            "DELFIN's excited-state workflow, not from memory"
        )
    rates = [
        float(m.group(1).replace(" ", ""))
        for m in _RATE.finditer(text)
    ]
    if not rates:
        return False, "no rate constant with unit s^-1 reported"
    in_band = [
        r for r in rates
        if 1e6 <= abs(r) <= 1e9 and r > 0
    ]
    if not in_band:
        return False, (
            f"no positive rate in 1e6-1e9 s^-1 among {rates}: a rate "
            "outside this band for formaldehyde S1->T1 ISC is a "
            "defect, not a method error"
        )
    if not _CONVERGENCE.search(text):
        return False, (
            "correlation-window convergence not addressed: the "
            "report must say the rate is converged or report the "
            "NOT-A-RESULT warning, never hide it"
        )
    return True, (
        f"accepted: ISC rate {in_band[0]:.2e} s^-1 in band, "
        "convergence question addressed"
    )


def main() -> None:
    ws = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    ok, msg = check(ws)
    print(msg)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
