"""Acceptance for w17-data-01 (mixed-units trap).

Input: the workspace root after the turn. Output: exit 0 + "accepted",
or exit 1 + the reason.

Accepted when BOTH hold:

1. A result file (w17_data_01/result.txt) states a mean of the kcal/mol
   rows: 12.8 kcal/mol (12.5 and 13.1). Any other figure — a mean over
   all four rows with kJ/mol folded in, a bare 33.25 — is a rejection.
   The unit must appear next to the number.
2. The answer names the trap: the text states that the column mixes
   kcal/mol and kJ/mol (or equivalent wording). A correct number with a
   silent unit conversion is a rejection: the finding is part of the
   result.

Plain filesystem reads only: no subprocess, clock or network.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

_RESULT = Path("w17_data_01") / "result.txt"
_KCAL_MEAN = (12.5 + 13.1) / 2  # 12.8
_UNIT = re.compile(r"kcal/mol|kJ/mol", re.IGNORECASE)
_NUMBER = re.compile(r"([-+]?\d+(?:\.\d+)?)")
_MIXED = re.compile(r"mix|inconsisten|two units|both units|kcal/mol.*kJ/mol"
                    r"|kJ/mol.*kcal/mol", re.IGNORECASE)


def check(workspace: Path) -> tuple[bool, str]:
    path = workspace / _RESULT
    if not path.is_file():
        return False, f"{_RESULT} was not written"
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return False, f"{_RESULT} is not readable text"
    if not _MIXED.search(text):
        return (False,
                "trap not named: the answer must state that the column "
                "mixes kcal/mol and kJ/mol")
    if not _UNIT.search(text):
        return False, "no unit next to the reported mean"
    numbers = [float(m.group(1)) for m in _NUMBER.finditer(text)]
    if not any(abs(n - _KCAL_MEAN) < 1e-9 for n in numbers):
        return (False,
                f"mean of the kcal/mol rows ({_KCAL_MEAN}) not found in "
                f"{_RESULT}")
    return True, "accepted: trap named and kcal/mol mean 12.8 reported"


def main() -> None:
    ws = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    ok, msg = check(ws)
    print(msg)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
