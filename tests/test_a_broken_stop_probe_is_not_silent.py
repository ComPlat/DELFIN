"""A stop probe that raises is said out loud, once.

Found from outside: a session testing today's Stop mechanism wrote a
probe that defined ``__bool__`` and not ``__call__``. Calling it raised
``TypeError``, the wait swallowed the exception and read it as "do not
stop" — at every poll, for the whole command. The one failure a stop
mechanism must not have is the quiet one.

It is not turned into a stop: a flaky probe would then end every command
that runs, which is worse than one that cannot be stopped. It is
reported once, and not asked again — the command runs unstoppable rather
than unstoppable and unremarked.

And a test run's verdict is kept whole in the trace: at the ordinary
2000-character cap about half of the recorded ones were cut mid-JSON,
which is what the flake report had to count as unreadable.
"""

from __future__ import annotations

import logging
import threading
import time

import pytest

from delfin.agent import contained_run, tool_trace


class _NotCallable:
    """A probe of the shape that caused this: truthy, not callable."""

    def __bool__(self) -> bool:
        return True


def test_a_probe_that_raises_is_reported(caplog):
    with caplog.at_level(logging.WARNING, logger="delfin.agent.contained_run"):
        result = contained_run.run(["sleep", "0.6"], timeout=10,
                                   should_stop=_NotCallable())
    said = [r for r in caplog.records if "stop probe raised" in r.getMessage()]
    assert len(said) == 1, "said once, not at every poll"
    assert "can no longer be stopped" in said[0].getMessage()
    assert result.returncode == 0, "a broken probe must not end the command"


def test_a_probe_that_raises_does_not_end_the_command():
    """The safety direction does not flip: a probe nobody can ask is not
    an instruction to kill what is running."""
    calls = {"n": 0}

    def _raises():
        calls["n"] += 1
        raise RuntimeError("no")

    result = contained_run.run(["echo", "still here"], timeout=10,
                               should_stop=_raises)
    assert result.returncode == 0
    assert result.stdout.strip() == "still here"
    assert calls["n"] <= 1, "asked once, then left alone"


def test_a_working_probe_still_stops_the_command():
    flag = {"stop": False}

    def _flip():
        time.sleep(0.3)
        flag["stop"] = True

    threading.Thread(target=_flip, daemon=True).start()
    started = time.monotonic()
    result = contained_run.run(["sleep", "20"], timeout=20,
                               should_stop=lambda: flag["stop"])
    assert result.returncode == contained_run.STOPPED_RETURNCODE
    assert time.monotonic() - started < 5


# -- the verdict a later reader needs -------------------------------------

def test_a_test_runs_verdict_is_kept_whole(tmp_path, monkeypatch):
    """Its result is a verdict with a failure list, and a diagnostic that
    reads it back cannot use half a JSON object."""
    monkeypatch.setattr(tool_trace, "_DIR", tmp_path)
    long_result = {"status": "failed", "failures": [f"tests/t{i}.py::x"
                                                    for i in range(400)]}
    import json
    tool_trace.record("s1", tool="run_tests", tool_input="{}",
                      output=json.dumps(long_result), ok=False)
    entries = tool_trace.read("s1", root=tmp_path)
    assert entries, "the call was recorded"
    assert json.loads(entries[0]["output"])["failures"], "and it parses"


def test_an_ordinary_tool_is_still_capped(tmp_path, monkeypatch):
    monkeypatch.setattr(tool_trace, "_DIR", tmp_path)
    tool_trace.record("s2", tool="bash", tool_input="{}",
                      output="x" * 50_000, ok=True)
    entries = tool_trace.read("s2", root=tmp_path)
    assert 0 < len(entries[0]["output"]) <= 2000
