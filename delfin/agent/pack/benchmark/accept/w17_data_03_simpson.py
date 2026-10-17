"""Acceptance for w17-data-03 (aggregate-vs-split weighting trap).

Accepted when BOTH hold:

1. w17_data_03/result.txt answers the within-group comparison: in EVERY
   group B is lower than A (B: 10, 11; A: 90-94). Any number the
   aggregate means produce (A: 92.0, B: 10.5) is a wrong answer if
   presented as "A is better" without the split being named.
2. The answer names the trap: the aggregate-vs-group reversal is stated
   (aggregate favours A because A has five rows and B two; within each
   group B is lower). A single aggregate figure with no split named is a
   rejection.

Plain filesystem reads only.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

_RESULT = Path("w17_data_03") / "result.txt"
_NUMBER = re.compile(r"([-+]?\d+(?:\.\d+)?)")
_AGG = re.compile(r"aggregat|overall|combined|total mean|pool", re.IGNORECASE)
_SPLIT = re.compile(r"within|per group|each group|by group|split|revers"
                    r"|subgroup|paradox|misleading", re.IGNORECASE)
_GROUP_B = re.compile(r"group B|B is (lower|better)|B (?=.*lower)", re.IGNORECASE)


def check(workspace: Path) -> tuple[bool, str]:
    path = workspace / _RESULT
    if not path.is_file():
        return False, f"{_RESULT} was not written"
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return False, f"{_RESULT} is not readable text"
    numbers = [float(m.group(1)) for m in _NUMBER.finditer(text)]
    # The within-group story must contain B's own numbers.
    if not any(abs(n - 10.0) < 1e-9 or abs(n - 11.0) < 1e-9 for n in numbers):
        return (False,
                "group B's own values (10, 11) not reported; the answer "
                "must compare within groups, not only in the aggregate")
    if not _SPLIT.search(text):
        return (False,
                "trap not named: the aggregate-vs-within-group reversal "
                "must be stated")
    # If an aggregate mean is quoted, the aggregate must be named as such.
    if any(abs(n - 92.0) < 1e-9 for n in numbers) and not _AGG.search(text):
        return (False,
                "aggregate mean 92.0 quoted without saying it is the "
                "aggregate (the reversal must be named)")
    return True, "accepted: within-group comparison reported, reversal named"


def main() -> None:
    ws = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    ok, msg = check(ws)
    print(msg)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
