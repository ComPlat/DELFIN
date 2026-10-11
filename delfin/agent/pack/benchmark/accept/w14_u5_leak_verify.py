"""Acceptance for w14_u5_dsml_leak_is_reissued_as_a_real_call.

Input: the workspace root after the turn (argv[1], or cwd), as
benchmark_runner.run_acceptance passes it.
Output: exit 0 and an "accepted" line, or exit 1 and the reason.

Accepted when both hold:

1. w14_u5_readings/level.txt holds one number equal to 41.7, optionally
   followed by the unit "m". The prompt does not state the value, so only a
   read of readings.txt that ran, followed by a write that ran, puts it
   there; a call written as text runs nothing. Formatting (unit, trailing
   zero, whitespace) is not judged.
2. w14_u5_readings/readings.txt is unchanged from the fixture, so the
   value was not planted by rewriting the source.

An unreadable or undecodable file is a rejection with a message, not a
traceback. Plain filesystem reads only: no subprocess, clock or network.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

from delfin.agent.pack.benchmark.setup.w14_u5_leak_setup import (
    _READINGS,
    _VALUE,
)

_FOLDER = "w14_u5_readings"
_LEVEL = "level.txt"
_SOURCE = "readings.txt"
_NUMBER_WITH_UNIT = re.compile(r"\s*([-+]?\d+(?:\.\d+)?)\s*(?:m)?\s*")


def _read(path: Path) -> str | None:
    try:
        return path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return None


def check(workspace: Path) -> tuple[bool, str]:
    """``(accepted, message)`` for one workspace."""
    d = workspace / _FOLDER
    source = d / _SOURCE
    level = d / _LEVEL

    if not source.is_file():
        return False, "readings.txt is missing (fixture not built)"
    if _read(source) != _READINGS:
        return (False,
                "readings.txt differs from the fixture; the source must stay "
                "unchanged")

    if not level.is_file():
        return (False,
                "level.txt was not written; no call that ran produced it")
    body = _read(level)
    if body is None:
        return False, "level.txt could not be read as UTF-8 text"
    m = _NUMBER_WITH_UNIT.fullmatch(body)
    if not m or abs(float(m.group(1)) - float(_VALUE)) > 1e-9:
        return (False,
                f"level.txt holds {body.strip()[:80]!r}, not the number "
                f"{_VALUE}")

    return True, ("accepted: level.txt holds 41.7 and readings.txt is "
                  "unchanged")


def main() -> None:
    ws = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    ok, msg = check(ws)
    print(msg)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
