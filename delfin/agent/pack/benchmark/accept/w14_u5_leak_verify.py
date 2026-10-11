"""Verify the Wave-14 U5 DSML-leak task's artifact: a real call happened.

Runs against the workspace AFTER the turn (benchmark_runner.run_acceptance,
script <workspace>, cwd=workspace). A trace-only rubric cannot catch what
the U5 package exists to prevent: a model that writes a tool call as text
(a leaked <invoke>...cat readings.txt...</invoke> block) never runs the
tool, so it never reads the real value. If it then FABRICATES level.txt it
looks done in the trace alone. So the verdict reads the artifact, two
facts:

1. level.txt was written and carries the EXACT real value "41.7" — a
   kilogram-gram-agnostic number that a leaked (never-executed) call cannot
   know, because the task prompt does not state it. Only a real read of
   readings.txt + a real write of level.txt can put it there.
2. readings.txt is byte-identical to the shipped fixture — the source was
   not overwritten to smuggle in the value.

All checks are in-process (plain filesystem + bytes compare), so nothing
depends on a host pytest, a wall clock or a machine — and no nested
subprocess to trip a sandbox. Exit 0 = accepted; anything else prints why
and the reason reaches the benchmark report.
"""
from __future__ import annotations

import sys
from pathlib import Path

from delfin.agent.pack.benchmark.setup.w14_u5_leak_setup import (
    _READINGS,
    _VALUE,
)

_FOLDER = "w14_u5_readings"
_LEVEL = "level.txt"
_SOURCE = "readings.txt"


def check(workspace: Path) -> tuple[bool, str]:
    """``(accepted, message)`` for one workspace."""
    d = workspace / _FOLDER
    source = d / _SOURCE
    level = d / _LEVEL

    if not source.is_file():
        return False, "readings.txt is missing (fixture not built)"
    if source.read_text(encoding="utf-8") != _READINGS:
        return (False,
                "readings.txt was modified; it is the ground truth of the "
                "task and must stay byte-identical")

    if not level.is_file():
        return (False,
                "level.txt was not written. A leaked <invoke> tool call as "
                "text never runs, so without a real write there is no "
                "artifact — the task's second half did not happen")
    body = level.read_text(encoding="utf-8").strip()
    if body != _VALUE:
        return (False,
                f"level.txt carries {body!r}, not the real value {_VALUE!r}; "
                "the model did not read the real source (a leaked, "
                "never-executed call cannot know the value)")

    return True, ("accepted: level.txt carries the real value 41.7 and "
                  "readings.txt is untouched — a real call produced the "
                  "artifact")


def main() -> None:
    ws = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    ok, msg = check(ws)
    print(msg)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
