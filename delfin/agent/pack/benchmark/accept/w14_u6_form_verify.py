"""Verify the Wave-14 U6 refusal-form task's artifact: the allowed form
actually produced the result.

Runs against the workspace AFTER the turn (benchmark_runner.run_acceptance,
script <workspace>, cwd=workspace). A trace-only rubric cannot prove that a
refused command shape was not used and the real result found; it can only
see what the model chose to show. So the verdict reads the artifact, two
facts:

1. winner.txt was written and carries the EXACT winning entry "gamma 9.84".
   The prompt does not state the winner, so it must have been found by
   actually reading log.txt and comparing the values — a refused shape
   (awk / $( ) / xargs / backtick-pipe / sed -i) that never ran could not
   produce this artifact.
2. log.txt is byte-identical to the shipped fixture — the source was not
   overwritten to smuggle in the answer.

All checks are in-process (plain filesystem + bytes compare), so nothing
depends on a host pytest, a wall clock or a machine — and no nested
subprocess to trip a sandbox. Exit 0 = accepted; anything else prints why
and the reason reaches the benchmark report.
"""
from __future__ import annotations

import sys
from pathlib import Path

from delfin.agent.pack.benchmark.setup.w14_u6_form_setup import (
    _LOG,
    _MAX_ITEM,
    _MAX_VALUE,
)

_FOLDER = "w14_u6_entries"
_WINNER = "winner.txt"
_SOURCE = "log.txt"
_EXPECTED_WINNER = f"{_MAX_ITEM} {_MAX_VALUE}"


def check(workspace: Path) -> tuple[bool, str]:
    """``(accepted, message)`` for one workspace."""
    d = workspace / _FOLDER
    source = d / _SOURCE
    winner = d / _WINNER

    if not source.is_file():
        return False, "log.txt is missing (fixture not built)"
    if source.read_text(encoding="utf-8") != _LOG:
        return (False,
                "log.txt was modified; it is the ground truth of the task "
                "and must stay byte-identical")

    if not winner.is_file():
        return (False,
                "winner.txt was not written. A refused command shape never "
                "runs, so without the allowed form there is no artifact — "
                "the real comparison did not happen")
    body = " ".join(winner.read_text(encoding="utf-8").split())
    if body != _EXPECTED_WINNER:
        return (False,
                f"winner.txt carries {body!r}, not the real winning entry "
                f"{_EXPECTED_WINNER!r}; the entries were not read and "
                "compared correctly")

    return True, ("accepted: winner.txt carries the real winning entry "
                  "gamma 9.84 and log.txt is untouched — the allowed form "
                  "produced the artifact")


def main() -> None:
    ws = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    ok, msg = check(ws)
    print(msg)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
