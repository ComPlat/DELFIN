"""The teardown guard must not book a parallel neighbour's work.

The SLURM suite runs one pytest process per test file, sixteen at a
time, all in ONE checkout. The session guard in conftest compares the
whole tree before and after ITS OWN run -- and everything a neighbour
does to the shared tree in that window arrives as a "path appeared in
the checkout" error, booked on the last test of this file:

* runs 7199892/7214159/7210626/7225552: tracked files under
  tests/fixtures/{behavior,user_project,science}_workspace booked as
  new paths. The benchmark guard (_PristineWorkspace) removes those
  workspaces and copies a snapshot back per attempt, so a walk taken
  mid-attempt missed the files and a later walk "found" them. The
  office workspace got the generated-root treatment for exactly this
  after SLURM 7188718; the other three never did.
* run 7214160: ``.delfin_leak_probe`` -- a transient directory another
  process's test creates at the checkout top and removes in its
  finally. Booked because it existed for the instant of this run's
  teardown walk.
* coordinator reproduction (s13, 2026-09-26): the operator drops new
  SLURM logs into ``.gate/rot/logs/`` while tests run.

What must NOT change: a path that persists is reported exactly as
before -- a real leak in a suite process, or leftover damage inside a
fixture workspace, still fails the run.
"""

from __future__ import annotations

import sys
import threading
import time
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

import conftest as cf  # noqa: E402


# The three fixture workspaces the benchmark guard snapshots and
# restores per attempt -- the reason they cannot be plain tree.
_RESTORED_WS = ("tests/fixtures/behavior_workspace",
                "tests/fixtures/user_project_workspace",
                "tests/fixtures/science_workspace")


def _make_repo(tmp_path: Path) -> Path:
    """A checkout shaped like the real one: tracked files inside the
    restored fixture workspaces, a git repository around them."""
    import subprocess
    root = tmp_path / "repo"
    (root / "delfin").mkdir(parents=True)
    (root / "delfin" / "__init__.py").write_text("", encoding="utf-8")
    for rel in _RESTORED_WS:
        (root / rel).mkdir(parents=True)
        (root / rel / "tracked.txt").write_text("a\n", encoding="utf-8")
    for cmd in (["git", "init", "-q"], ["git", "add", "-A"],
                ["git", "-c", "user.email=t@t", "-c", "user.name=t",
                 "commit", "-q", "-m", "fixtures"]):
        subprocess.run(cmd, cwd=root, check=True, capture_output=True)
    return root


@pytest.fixture()
def repo(tmp_path, monkeypatch):
    root = _make_repo(tmp_path)
    monkeypatch.setattr(cf, "_CHECKOUT_ROOT", root)
    # A short confirmation ladder: the gaps only pace how long a
    # transient is waited out, not what is reported.
    monkeypatch.setattr(cf, "_CONFIRM_GAPS_S", (0.05, 0.1, 0.2))
    return root


def _teardown_diff(repo, before):
    """The computation the session fixture's assert runs, against the
    CURRENT conftest: with the factored function when the fix is in,
    and with the pre-fix shape (raw diff plus the generated-root
    confirmation) when it is not, so the control runs red on the
    previous commit for the behaviour and not for a missing name."""
    if hasattr(cf, "_teardown_new_paths"):
        return cf._teardown_new_paths(before, frozenset())
    new = sorted(cf._checkout_entries() - before)
    new += sorted(
        str(cf._CHECKOUT_ROOT / p)
        for p in cf._confirmed_under_a_generated_root() - frozenset())
    return new


def test_a_restored_fixture_file_is_not_booked(repo, monkeypatch):
    """The incident: this run's before-walk falls into a neighbour's
    mid-attempt window (the tracked file is gone), the neighbour
    restores it, and the teardown walk finds it -- as a "new" path,
    because the before-walk never saw it."""
    tracked = repo / _RESTORED_WS[0] / "tracked.txt"
    before = cf._checkout_entries()
    # The neighbour's attempt holds the workspace: the file is away.
    tracked.unlink()
    before = cf._checkout_entries()
    assert str(tracked) not in before, "precondition: the window is open"
    # The neighbour's attempt ends: the snapshot goes back.
    tracked.write_text("a\n", encoding="utf-8")
    reported = _teardown_diff(repo, before)
    assert not any("tracked.txt" in p for p in reported), sorted(reported)


def test_operator_writes_under_gate_are_not_booked(repo):
    """Logs the operator drops into .gate/ while a run is going are
    coordination data, not a suite leak."""
    before = cf._checkout_entries()
    (repo / ".gate" / "rot" / "logs").mkdir(parents=True)
    (repo / ".gate" / "rot" / "logs" / "7234567_base_test_x.log").write_text(
        "ERROR ...\n", encoding="utf-8")
    reported = _teardown_diff(repo, before)
    assert not any(".gate" in p for p in reported), sorted(reported)


def test_a_transient_entry_at_the_top_is_waited_out(repo):
    """The .delfin_leak_probe incident: a directory that exists for an
    instant -- another process creates it and removes it in a finally
    -- must not survive the confirmation window."""
    probe = repo / ".delfin_leak_probe"
    before = cf._checkout_entries()
    probe.mkdir()
    assert str(probe) in cf._checkout_entries() - before, (
        "precondition: the walk sees the probe")
    vanish = threading.Timer(0.05, probe.rmdir)
    vanish.start()
    reported = _teardown_diff(repo, before)
    vanish.join()
    assert not any(".delfin_leak_probe" in p for p in reported), (
        sorted(reported))


def test_a_persistent_entry_at_the_top_is_still_booked(repo):
    """The teeth: a path that is still there after the confirmation
    window is not a race and is reported exactly as before."""
    (repo / "delfin" / "agent").mkdir()
    (repo / "delfin" / "agent" / "leaked.txt").write_text("x\n",
                                                          encoding="utf-8")
    before = cf._checkout_entries()
    (repo / "delfin" / "agent" / "leaked2.txt").write_text("x\n",
                                                           encoding="utf-8")
    reported = _teardown_diff(repo, before)
    assert any(p.endswith("leaked2.txt") for p in reported), sorted(reported)


def test_leftover_damage_inside_a_restored_workspace_is_still_booked(repo):
    """Inside a declared root git's rules decide: an attempt's leftover
    that no rule describes survives the confirmation and is reported --
    exempting the restore race must not exempt real damage."""
    ws = repo / _RESTORED_WS[1]
    before = cf._checkout_entries()
    before_generated = cf._unexpected_under_a_generated_root()
    (ws / "attempt_leftover.out").write_text("junk\n", encoding="utf-8")
    reported = cf._teardown_new_paths(before, before_generated)
    assert any("attempt_leftover.out" in p for p in reported), (
        sorted(reported))
