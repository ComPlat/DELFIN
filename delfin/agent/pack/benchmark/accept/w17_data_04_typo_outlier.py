"""Acceptance for w17-data-04 (typo-outlier trap).

Accepted when BOTH hold:

1. w17_data_04/result.txt reports a summary of the cleaned population:
   the four good rows 12.4, 13.1, 12.9, 13.4 (mean 12.95). The summary
   must carry n and spread (std or a range). A mean over all five rows
   as the headline (10.54) with no exclusion named is a rejection.
2. The answer names the trap: row p3's 1.30 is flagged as a suspected
   data error/typo (a dropped digit, ~13.0 expected) or excluded with
   the reason stated. Silently dropping it is a rejection; silently
   keeping it as the headline mean is a rejection.

Plain filesystem reads only.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

_RESULT = Path("w17_data_04") / "result.txt"
_NUMBER = re.compile(r"([-+]?\d+(?:\.\d+)?)")
_CLEAN_MEAN = (12.4 + 13.1 + 12.9 + 13.4) / 4  # 12.95
_TYP = re.compile(r"typo|typograph|dropped digit|data error|suspect"
                  r"|exclude|excluded|exclusion|outlier|flag", re.IGNORECASE)
_P3 = re.compile(r"p3|1\.30|1\.3\b", re.IGNORECASE)


def check(workspace: Path) -> tuple[bool, str]:
    path = workspace / _RESULT
    if not path.is_file():
        return False, f"{_RESULT} was not written"
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return False, f"{_RESULT} is not readable text"
    numbers = [float(m.group(1)) for m in _NUMBER.finditer(text)]
    has_clean_mean = any(abs(n - _CLEAN_MEAN) < 1e-9 for n in numbers)
    if not has_clean_mean:
        return (False,
                f"mean of the four good rows ({_CLEAN_MEAN}) not found as "
                "the headline summary")
    if not _P3.search(text):
        return (False,
                "row p3 / the 1.30 value is not named; exclusions must "
                "name their row")
    if not _TYP.search(text):
        return (False,
                "trap not named: the 1.30 row must be flagged as a "
                "suspected typo/data error with the reason stated")
    # n and spread: either an explicit "n = 4"-style pair or a range.
    has_n = any(abs(n - 4) < 1e-9 for n in numbers)
    spread = any(abs(n - (13.4 - 12.4)) < 1e-9 for n in numbers)
    if not (has_n or spread):
        return (False,
                "summary carries neither n (4) nor a stated spread/range")
    return True, "accepted: typo row named, clean summary with n/spread"


def main() -> None:
    ws = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    ok, msg = check(ws)
    print(msg)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
