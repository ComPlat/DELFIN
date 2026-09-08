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


def test_a_nested_guard_does_not_lose_the_original_handler(tmp_path):
    before = signal.getsignal(signal.SIGTERM)
    with _PristineWorkspace(tmp_path):
        with _PristineWorkspace(tmp_path):
            pass
    assert signal.getsignal(signal.SIGTERM) == before


def test_the_files_come_back_when_the_block_is_left_by_an_exception(tmp_path):
    """The property the handler exists to reach, exercised through the
    exception it raises."""
    ws = tmp_path / "office_workspace"
    ws.mkdir()
    (ws / "buchungen.csv").write_text("id;betrag\n1;10\n", encoding="utf-8")
    with pytest.raises(KeyboardInterrupt):
        with _PristineWorkspace(tmp_path):
            (ws / "buchungen.csv").unlink()
            raise KeyboardInterrupt
    assert (ws / "buchungen.csv").is_file(), "the fixture was not restored"
