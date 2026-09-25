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

import os
import subprocess
import sys
import threading
import time
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


# ---------------------------------------------------------------------------
# The incident as it actually happened: two processes, not a timer
# ---------------------------------------------------------------------------

_NEIGHBOUR = r"""
import os
import sys
import time
from pathlib import Path

# A process shaped like the benchmark runner: hold the workspace in
# its mid-attempt state (tracked files removed, as the agent under
# test removes them) and restore it after the guard below has scanned.
root = Path(sys.argv[1])
ws = root / "tests" / "fixtures" / "office_workspace"
kept = {p.name: p.read_bytes() for p in ws.iterdir()
        if p.name.endswith(".csv")}

from delfin.agent import benchmark_runner as br

with br._PristineWorkspace(root):
    for p in ws.glob("*.csv"):
        p.unlink()
    # Hold the window open until the guard has taken its raw scan;
    # the file this process waits for says so.
    (ws / ".neighbour-was-here").write_text("", encoding="utf-8")
    while not (root / ".guard-scan-started").exists():
        time.sleep(0.02)
    # A window the guard's confirmation gaps can outwait: 0.1 s after
    # the raw scan, the restore closes it -- the shape a snapshot
    # restore has when the attempt ends.
    time.sleep(0.1)
    for name, blob in kept.items():
        (ws / name).write_bytes(blob)
    (ws / ".neighbour-was-here").unlink()
os._exit(0)
"""


def test_the_incident_two_processes_no_false_error(tmp_path, monkeypatch):
    """The report's own reproduction: a parallel suite run, one pytest
    process per file. One child drives the real workspace guard the
    way a benchmark attempt does -- tracked CSVs gone for the length
    of the attempt, restored at the end -- while THIS process runs the
    guard's teardown scan, unmodified, over the shared checkout.

    Red on the previous commit: the raw scan booked the neighbour's
    transient deletions as appearances. Green with the confirmation:
    the neighbour's window is waited out, and no false error.
    """
    root = _make_repo(tmp_path)
    signal_scan = root / ".guard-scan-started"
    neighbour_log = tmp_path / "neighbour.log"
    # Extra paths the neighbour legitimately creates inside the
    # workspace during its attempt: the guard's generated-root scan
    # reports them only if they SURVIVE, which they do not.
    script = tmp_path / "neighbour.py"
    script.write_text(_NEIGHBOUR, encoding="utf-8")
    env = dict(os.environ)
    env["PYTHONPATH"] = str(Path(__file__).resolve().parents[1])
    # A private HOME for the child, the way conftest.child_env gives
    # one: _PristineWorkspace redirects its user state, and none of it
    # belongs in the real home.
    child_home = tmp_path / "neighbour_home"
    child_home.mkdir()
    env["HOME"] = str(child_home)
    neighbour = subprocess.Popen(
        [sys.executable, str(script), str(root)],
        stdout=subprocess.PIPE,
        stderr=open(neighbour_log, "w"), text=True,
        env=env)
    try:
        # Wait until the neighbour's window is open (its marker exists
        # and the tracked CSVs are gone), then run the guard's own
        # teardown path -- the raw scan first, for evidence, then the
        # confirmed one the fixture now uses.
        marker = root / _WS_REL / ".neighbour-was-here"
        for _ in range(500):            # up to ~10 s
            if marker.exists() and not (root / _WS_REL /
                                        "buchungen.csv").exists():
                break
            time.sleep(0.02)
        else:
            raise AssertionError("neighbour never opened its window: "
                                 + neighbour.poll().__repr__())
        assert f"{_WS_REL}/buchungen.csv" in _scan(monkeypatch, root), (
            "the raw scan no longer sees the neighbour's window -- "
            "the reproduction would be vacuous")
        # Tell the neighbour its window has been seen; the confirmed
        # scan runs while the window is still open and must wait it
        # out, reporting nothing that survives.
        signal_scan.write_text("", encoding="utf-8")
        reported = _confirmed(monkeypatch, root)
        rc = neighbour.wait(timeout=120)
        assert rc == 0, ("the neighbour crashed: "
                         + neighbour_log.read_text(encoding="utf-8"))
    finally:
        if signal_scan.exists():
            signal_scan.unlink()
        if neighbour.poll() is None:
            neighbour.kill()
            neighbour.wait(timeout=30)
    # The checkout is whole again and the guard reported nothing.
    assert reported == frozenset(), (
        f"a neighbour's restore window was booked as a leak: "
        f"{sorted(reported)}")
    assert (root / _WS_REL / "buchungen.csv").is_file()
    # The real teardown assert shape: nothing appeared AND nothing
    # under the generated root survived the neighbour.
    assert not [p for p in (root / _WS_REL).iterdir()
                if p.name == ".neighbour-was-here"]
