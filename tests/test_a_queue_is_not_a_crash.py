"""A benchmark that could not start printed a traceback, not a sentence.

`fixture_run_lock` takes the fixture workspaces exclusively and refuses
rather than waits, on a sound ground: two runs in one checkout measure
each other's writes. Its message says so, names the pid and model
holding the lock, and gives both ways out.

That message was reaching the caller at the bottom of a twenty-line
traceback. In a queued script's log that is where it is buried: three
arms of a GLM A/B died this way on 2026-09-10 because the script was
launched while a calibration run still held the lock, and the log read
as a crash rather than as a queue — so the arms looked like a broken
model instead of a scheduling mistake.

The lock's behaviour is unchanged. Only the rendering is.
"""

from __future__ import annotations

import pytest

from delfin.agent.benchmark_runner import (
    BenchmarkRunInProgress, fixture_run_lock)

def test_the_lock_still_refuses_a_second_holder(tmp_path):
    """The premise. If this ever starts waiting instead, the rendering
    below is answering a question nobody asks."""
    with fixture_run_lock(tmp_path, owner="first"):
        with pytest.raises(BenchmarkRunInProgress):
            with fixture_run_lock(tmp_path, owner="second"):
                pass


def test_the_refusal_names_the_holder_and_both_ways_out(tmp_path):
    with fixture_run_lock(tmp_path, owner="benchmark run, model X"):
        try:
            with fixture_run_lock(tmp_path, owner="second"):
                pass
        except BenchmarkRunInProgress as exc:
            msg = str(exc)
        else:
            pytest.fail("the lock did not refuse")
    assert "model X" in msg
    assert "wait" in msg and "separate checkout" in msg


def test_the_lock_is_released_again(tmp_path):
    with fixture_run_lock(tmp_path, owner="first"):
        pass
    with fixture_run_lock(tmp_path, owner="second"):
        pass


def test_a_collision_is_one_line_and_its_own_exit_code(monkeypatch, capsys):
    """Through the CLI, which is where it was ugly.

    Not by spawning a second process: tests/conftest.py redirects
    ``_RUN_LOCK_DIR`` per test so a suite never touches the user's real
    lock, so a parent and a child cannot contend for the same file by
    construction. What changed here is the RENDERING of an exception, so
    the exception is what gets raised.
    """
    from delfin.agent import benchmark_runner as _br
    from delfin.agent import cli as _cli

    held = ("the fixture workspaces are held by pid 4242: benchmark run, "
            "model kit.deepseek-v4-flash. Two runs in one checkout measure "
            "each other's writes -- wait for it, or use a separate checkout.")

    def _refuse(*a, **k):
        raise _br.BenchmarkRunInProgress(held)

    monkeypatch.setattr(_br, "run_suite", _refuse)
    rc = _cli.main(["bench", "run", "--model", "kit.test", "--provider", "kit",
                    "--task", "science_a_number_carries_its_unit_and_its_factor",
                    "--repeats", "1"])
    out = capsys.readouterr()
    both = out.out + out.err
    assert rc == 3, (rc, both[-400:])
    assert "Traceback" not in both
    assert "benchmark not started" in both
    assert held in both, both[-400:]


def test_the_run_is_not_recorded_when_it_never_started(monkeypatch, capsys):
    """A refused run must not write a result file: a baseline built from
    one would compare against a suite that did not execute."""
    from delfin.agent import benchmark as _bm
    from delfin.agent import benchmark_runner as _br
    from delfin.agent import cli as _cli

    def _refuse(*a, **k):
        raise _br.BenchmarkRunInProgress("held by pid 1: x")

    wrote = []
    monkeypatch.setattr(_br, "run_suite", _refuse)
    monkeypatch.setattr(_bm, "write_run",
                        lambda *a, **k: wrote.append(a) or "never")
    _cli.main(["bench", "run", "--model", "kit.test", "--provider", "kit",
               "--task", "science_a_number_carries_its_unit_and_its_factor",
               "--repeats", "1"])
    capsys.readouterr()
    assert not wrote


def test_a_run_that_starts_is_untouched(monkeypatch, capsys):
    """The other half — the guard must catch only this exception."""
    from delfin.agent import benchmark_runner as _br
    from delfin.agent import cli as _cli

    monkeypatch.setattr(_br, "run_suite",
                        lambda *a, **k: (_ for _ in ()).throw(
                            RuntimeError("something else entirely")))
    with pytest.raises(RuntimeError, match="something else entirely"):
        _cli.main(["bench", "run", "--model", "kit.test", "--provider", "kit",
                   "--task",
                   "science_a_number_carries_its_unit_and_its_factor",
                   "--repeats", "1"])
    capsys.readouterr()
