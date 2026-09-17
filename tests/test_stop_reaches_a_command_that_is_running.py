"""Stop reaches a command while it runs, not after it finishes.

The engine checks its stop flag between stream events, and the tool loop
checks it between rounds — but a single tool call that lasts half an
hour is inside neither. A session pressed Stop during a full suite run
and then "did not react any more": the worker was sitting in
``proc.wait(timeout=1800)`` and nothing in that wait ever looked
(2026-09-17).

  a running command is ended      the wait asks every quarter second
  and says it was stopped         exit 130, the shell's own convention
  a finished one is untouched     the ordinary path is unchanged
  the group goes with it          nothing survives the stop
  a test run says what it is      "stopped", not "failed": the suite's
                                  state does not follow from this
"""

from __future__ import annotations

import subprocess
import threading
import time

import pytest

from delfin.agent import contained_run


def _stop_after(seconds: float) -> "callable":
    flag = {"stop": False}

    def _flip():
        time.sleep(seconds)
        flag["stop"] = True

    threading.Thread(target=_flip, daemon=True).start()
    return lambda: flag["stop"]


def test_a_command_that_is_asked_to_stop_ends_at_once():
    started = time.monotonic()
    result = contained_run.run(["sleep", "30"], timeout=30,
                               should_stop=_stop_after(0.3))
    waited = time.monotonic() - started
    assert result.returncode == contained_run.STOPPED_RETURNCODE
    assert waited < 5, "the wait must not run to the command's own end"


def test_a_command_nobody_stops_runs_as_before():
    result = contained_run.run(["echo", "hi"], timeout=10,
                               should_stop=lambda: False)
    assert result.returncode == 0
    assert result.stdout.strip() == "hi"


def test_the_timeout_still_ends_a_command_with_a_probe():
    with pytest.raises(subprocess.TimeoutExpired):
        contained_run.run(["sleep", "10"], timeout=0.5,
                          should_stop=lambda: False)


def test_a_probe_that_raises_does_not_stop_anything():
    def _bad():
        raise RuntimeError("no")

    result = contained_run.run(["echo", "ok"], timeout=10, should_stop=_bad)
    assert result.returncode == 0


def test_whatever_the_command_started_ends_with_it(tmp_path):
    """The group, not just the process: a stop that leaves children
    behind is the containment hole this runner exists to close."""
    marker = tmp_path / "child-was-here"
    script = (f"sh -c 'sleep 30; touch {marker}' & sleep 30")
    result = contained_run.run(script, shell=True, timeout=30,
                               should_stop=_stop_after(0.3))
    assert result.returncode == contained_run.STOPPED_RETURNCODE
    time.sleep(0.5)
    assert not marker.exists()


def test_a_test_run_that_was_stopped_says_so(tmp_path):
    """Not "failed": nothing about the suite follows from a stop."""
    from delfin.agent import test_runner

    (tmp_path / "test_slow.py").write_text(
        "import time\n\n\ndef test_slow():\n    time.sleep(30)\n",
        encoding="utf-8")
    out = test_runner.run_tests(
        workspace=tmp_path, target="test_slow.py", timeout_s=30,
        should_stop=_stop_after(1.0))
    assert out["status"] == "stopped"
    assert "on request" in out.get("note", "")
