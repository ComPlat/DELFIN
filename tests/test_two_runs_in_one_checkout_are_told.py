"""Two runs in one checkout measure each other's writes.

Demonstrated twice, not argued.

2026-08-12: a full test-suite run failed two office fixture tests while a
benchmark was running, and ``git status`` caught it mid-attempt --
``D tests/fixtures/office_workspace/buchungen.csv``, deleted by the agent
under test and not yet restored by the fixture guard. Neither side was
broken. The measurement was, and it reported a product defect that did
not exist.

2026-08-13, the other direction: a suite and a live benchmark in one
checkout, and the benchmark's own "this run changed N file(s)" report
blamed itself for 569 files the suite had written.

The lock stops nobody. It names who is holding the fixtures, which is
the whole difference between a run that is wrong and a run that says why.
"""

from __future__ import annotations

import os
import threading

import pytest

from delfin.agent import benchmark_runner as br


def test_a_second_run_is_told_who_holds_the_fixtures(tmp_path):
    with br.fixture_run_lock(tmp_path, owner="the first run"):
        with pytest.raises(br.BenchmarkRunInProgress) as caught:
            with br.fixture_run_lock(tmp_path, owner="the second run"):
                pass
    message = str(caught.value)
    assert "the first run" in message
    assert str(os.getpid()) in message
    # ... and it says what to do about it.
    assert "separate checkout" in message


def test_the_lock_is_released_when_the_run_ends(tmp_path):
    with br.fixture_run_lock(tmp_path, owner="first"):
        pass
    # A held lock that outlived its run would make the next one impossible.
    with br.fixture_run_lock(tmp_path, owner="second"):
        pass


def test_the_lock_is_released_when_the_run_raises(tmp_path):
    with pytest.raises(ValueError):
        with br.fixture_run_lock(tmp_path, owner="first"):
            raise ValueError("the run blew up")
    with br.fixture_run_lock(tmp_path, owner="second"):
        pass


def test_exactly_one_of_many_contenders_gets_in(tmp_path):
    """Four processes' worth of threads, one winner, and everyone else
    told rather than silently interleaved."""
    got_in: list[int] = []
    refused: list[str] = []
    ready = threading.Barrier(4)

    def _contend(n: int) -> None:
        ready.wait()
        try:
            with br.fixture_run_lock(tmp_path, owner=f"run {n}"):
                got_in.append(n)
                import time
                time.sleep(0.3)
        except br.BenchmarkRunInProgress as exc:
            refused.append(str(exc))

    threads = [threading.Thread(target=_contend, args=(n,))
               for n in range(4)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    assert len(got_in) == 1, got_in
    assert len(refused) == 3
    assert all("held by" in r for r in refused)


def test_the_lock_file_lives_beside_the_fixtures(tmp_path):
    """Not in /tmp and not in the user's home: it belongs to the checkout
    whose fixtures it protects, because that is what two runs share."""
    with br.fixture_run_lock(tmp_path, owner="first"):
        pass
    assert (tmp_path / "tests" / "fixtures" / ".benchmark-run.lock").exists()


def test_a_filesystem_without_flock_measures_rather_than_refuses(tmp_path,
                                                                 monkeypatch):
    """A measurement aid must not become a reason not to measure."""
    import fcntl

    def _no_flock(*a, **k):
        raise OSError(95, "operation not supported")

    monkeypatch.setattr(fcntl, "flock", _no_flock)
    ran = False
    with br.fixture_run_lock(tmp_path, owner="first"):
        ran = True
    assert ran
