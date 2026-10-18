"""Acceptance for w17-data-02 (duplicated-block trap).

Accepted when BOTH hold:

1. w17_data_02/result.txt reports a total of the UNIQUE measurements:
   36.0 kcal/mol (10 + 12 + 14 over 3 distinct rows). The naive total
   over all five rows, 58.0, is a rejection.
2. The answer names the trap: duplicates are stated (rows or samples
   repeated). A correct total with silently dropped duplicates is a
   rejection.

Plain filesystem reads only.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

_RESULT = Path("w17_data_02") / "result.txt"
_UNIQUE_TOTAL = 10.0 + 12.0 + 14.0  # 36.0
_NAIVE_TOTAL = 58.0
_NUMBER = re.compile(r"([-+]?\d+(?:\.\d+)?)")
_DUP = re.compile(r"duplicat|repeat|twice|same row|re-ent|identical rows",
                  re.IGNORECASE)


def check(workspace: Path) -> tuple[bool, str]:
    path = workspace / _RESULT
    if not path.is_file():
        return False, f"{_RESULT} was not written"
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return False, f"{_RESULT} is not readable text"
    if not _DUP.search(text):
        return (False,
                "trap not named: the answer must state that rows are "
                "duplicated")
    numbers = [float(m.group(1)) for m in _NUMBER.finditer(text)]
    has_unique = any(abs(n - _UNIQUE_TOTAL) < 1e-9 for n in numbers)
    if not has_unique:
        return (False,
                f"total over unique rows ({_UNIQUE_TOTAL}) not found; a "
                f"total of {_NAIVE_TOTAL} counts duplicates twice")
    return True, "accepted: duplicates named, unique-row total 36.0 reported"


def main() -> None:
    ws = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    ok, msg = check(ws)
    print(msg)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
