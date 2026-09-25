"""A neighbour restoring the fixture workspaces is not a leak.

SLURM job 7188718 (2026-09-25): the parallel suite run, one pytest
process per test file, ended eleven files with "N passed, 1 error" --
the error always in the teardown of the checkout guard in conftest,
reporting paths under ``tests/fixtures/office_workspace`` as "appeared
in the checkout during the run". The same files were green
individually (116 passed).

What the guard's git-status scan saw was not a leak. The benchmark
workspace guard (``_PristineWorkspace`` in
delfin/agent/benchmark_runner.py) restores the office workspace by
removing the directory and copying a snapshot back -- and during an
attempt the tracked CSVs can be gone for the length of that attempt,
because the agent under test deletes them as part of its task and the
restore happens only when the attempt ends. A neighbouring pytest
process whose teardown scan falls into that window sees `` D`` entries
for tracked files and booked them as paths that had appeared.

The fix keeps the finding and questions the moment: a deletion or a
modification of a tracked file is real damage and stays a finding --
but the scan is CONFIRMED by rescanning with growing gaps before it
reports, so a change that a neighbour's restore window put there for
an instant is waited out and a change that is still there after ~2 s
is not a race. A genuinely new untracked path is reported as before.
"""

from __future__ import annotations

import subprocess
import sys
import threading
from pathlib import Path

import pytest

import conftest as cf

_WS_REL = "tests/fixtures/office_workspace"
_TRACKED = ("buchungen.csv", "inventar.csv", "kostenstellen_roh.csv",
            "rechnungen.csv")


def _make_repo(tmp_path: Path) -> Path:
    """A checkout shaped like the real one: a generated root holding
    tracked files plus a local .gitignore for the workbooks."""
    root = tmp_path / "repo"
    ws = root / _WS_REL
    ws.mkdir(parents=True)
    (ws / ".gitignore").write_text("*.xlsx\n.fixture-stamp\n",
                                   encoding="utf-8")
    for name in _TRACKED:
        (ws / name).write_text("a;b\n", encoding="utf-8")
    for cmd in (["git", "init", "-q"],
                ["git", "add", "-A"],
                ["git", "-c", "user.email=t@t", "-c", "user.name=t",
                 "commit", "-q", "-m", "fixtures"]):
        subprocess.run(cmd, cwd=root, check=True,
                       capture_output=True)
    return root


def _scan(monkeypatch, root: Path) -> frozenset:
    """The guard's own raw git-status scan, pointed at ``root``.

    The function reads its checkout root from a module global at call
    time, so pointing it at a scratch repo needs no signature change
    and runs the very code the incident ran.
    """
    monkeypatch.setattr(cf, "_CHECKOUT_ROOT", root)
    return cf._unexpected_under_a_generated_root()


def _confirmed(monkeypatch, root: Path) -> frozenset:
    """The confirmed scan the teardown uses, pointed at ``root``."""
    monkeypatch.setattr(cf, "_CHECKOUT_ROOT", root)
    return cf._confirmed_under_a_generated_root()


def test_a_deletion_inside_a_neighbours_restore_window_is_waited_out(
        tmp_path, monkeypatch):
    """The incident, as a race: the tracked file is gone when the scan
    runs and back shortly after, the way a neighbouring process
    restoring a snapshot makes it.

    Red on the previous commit: there was no confirmation, so the raw
    scan's word was final and the deletion was reported -- eleven
    teardown errors in the parallel run.
    """
    root = _make_repo(tmp_path)
    victim = root / _WS_REL / "buchungen.csv"
    gone = victim.read_bytes()          # what the neighbour restores
    victim.unlink()                     # the neighbour's window opens
    # Evidence the window is real: the raw scan sees the deletion.
    assert f"{_WS_REL}/buchungen.csv" in _scan(monkeypatch, root)
    # The guard's confirmed scan runs while the window is open; the
    # neighbour puts the file back 0.15 s in, well inside the
    # confirmation gaps (0.05 + 0.2 s), the way a copytree restore
    # closes its window.
    seen: list = []

    def guard():
        seen.append(_confirmed(monkeypatch, root))

    thread = threading.Thread(target=guard)
    thread.start()
    monkeypatch.undo()
    monkeypatch.setattr(cf, "_CHECKOUT_ROOT", root)
    threading.Timer(0.15, lambda: victim.write_bytes(gone)).start()
    thread.join(timeout=30)
    assert victim.read_bytes() == gone
    assert seen == [frozenset()], (
        f"a transient deletion was booked as a leak: {seen}")


def test_a_deletion_that_stays_is_still_a_leak(tmp_path, monkeypatch):
    """The control the fix must not weaken: a tracked file deleted for
    good -- a test that damages versioned files -- is reported even
    after the full confirmation window. Red is the correct answer
    here, before the fix and after it."""
    root = _make_repo(tmp_path)
    (root / _WS_REL / "buchungen.csv").unlink()
    reported = _confirmed(monkeypatch, root)
    assert f"{_WS_REL}/buchungen.csv" in reported


def test_a_file_no_ignore_rule_describes_is_still_reported(
        tmp_path, monkeypatch):
    """A genuinely new, non-ignored path inside a generated root stays
    a finding through the confirmation. Green before the fix and after
    it -- this is the control that the confirmation narrows the scan's
    timing, not its scope."""
    root = _make_repo(tmp_path)
    (root / _WS_REL / "not_described_by_any_rule.txt").write_text(
        "leak\n", encoding="utf-8")
    reported = _confirmed(monkeypatch, root)
    assert f"{_WS_REL}/not_described_by_any_rule.txt" in reported
