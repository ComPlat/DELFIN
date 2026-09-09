"""SIGTERM does not run a `finally`.

The pristine guard snapshots the fixture directories before an attempt
and restores them after, so a benchmark task that edits, moves or deletes
the workbooks leaves nothing behind. It restores in `__exit__`, which
Python runs when the block is left — including on KeyboardInterrupt,
because that is an exception.

SIGTERM is not. `timeout 3600 delfin-agent bench run` sends it, the
interpreter exits where it stands, `__exit__` never runs, and the
directory keeps whatever the agent had done to it. Observed 2026-09-08:
a run killed by its own timeout during an office task left
buchungen.csv, rechnungen.csv, kostenstellen_roh.csv, inventar.csv and
README.md deleted from the checkout, and a half-written .wb-*.xlsx
beside them.

The same shape reaches the user: a scheduled run killed by its runner, a
container stopped, a `kill` on a long benchmark. One handler turns the
signal into the exception the guard already survives.
"""

from __future__ import annotations

import signal

import pytest

from delfin.agent.benchmark_runner import _PristineWorkspace


def test_a_term_signal_becomes_an_exception_the_guard_survives(tmp_path):
    """Inside the guard, SIGTERM raises rather than exiting where it
    stands — which is what lets `__exit__` put the files back."""
    before = signal.getsignal(signal.SIGTERM)
    raised = False
    try:
        with _PristineWorkspace(tmp_path):
            handler = signal.getsignal(signal.SIGTERM)
            assert callable(handler), "no handler installed"
            try:
                handler(signal.SIGTERM, None)
            except KeyboardInterrupt:
                raised = True
    except KeyboardInterrupt:
        raised = True
    assert raised, "SIGTERM did not become an exception"
    assert signal.getsignal(signal.SIGTERM) == before, (
        "the previous handler was not put back")


def test_the_handler_is_removed_again_afterwards(tmp_path):
    before = signal.getsignal(signal.SIGTERM)
    with _PristineWorkspace(tmp_path):
        pass
    assert signal.getsignal(signal.SIGTERM) == before


def test_a_nested_guard_stands_down_instead_of_deadlocking(tmp_path):
    """The lock is an exclusive flock, and a second guard opened a second
    descriptor on the same file and waited for a lock its own process
    already held — a silent deadlock, not an error.

    It cost a CI run: the job sat there until the 25-minute limit
    cancelled it, with two tests passed and nothing to say why. Nesting is
    a programming mistake either way; hanging is the worst way to report
    one, so the inner guard stands down and the outer keeps the lock and
    the restore.
    """
    ws = tmp_path / "tests" / "fixtures" / "office_workspace"
    ws.mkdir(parents=True)
    (ws / "buchungen.csv").write_text("id;betrag\n1;10\n", encoding="utf-8")
    before = signal.getsignal(signal.SIGTERM)
    with _PristineWorkspace(tmp_path):
        with _PristineWorkspace(tmp_path):
            (ws / "buchungen.csv").unlink()
        # The inner guard restored nothing: that is the outer one's job.
        assert not (ws / "buchungen.csv").exists()
    assert (ws / "buchungen.csv").is_file(), "the outer guard did not restore"
    assert signal.getsignal(signal.SIGTERM) == before


def test_the_files_come_back_when_the_block_is_left_by_an_exception(tmp_path):
    """The property the handler exists to reach, exercised through the
    exception it raises."""
    ws = tmp_path / "tests" / "fixtures" / "office_workspace"
    ws.mkdir(parents=True)
    (ws / "buchungen.csv").write_text("id;betrag\n1;10\n", encoding="utf-8")
    with pytest.raises(KeyboardInterrupt):
        with _PristineWorkspace(tmp_path):
            (ws / "buchungen.csv").unlink()
            raise KeyboardInterrupt
    assert (ws / "buchungen.csv").is_file(), "the fixture was not restored"
