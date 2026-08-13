"""Nothing under ~/.delfin ever expired.

``cleanup_old_sessions`` existed and had no production caller. There was no
prune at all for sub-agent sessions, running markers, tool traces, turn
metrics, transcript archives, handoffs, bundles or bug reports. Observed on
the audited machine: session files back to April (143 of 154 from one
month), 1748 tool traces, 2076 turn-metric files.

Retention is the amplifier. A directory that is briefly group-readable is a
one-day problem; a directory that is group-readable AND keeps everything is
a four-month archive of transcripts and raw tool output. Fixing the modes
without fixing the growth leaves the second half standing.

Everything here runs against a FIXTURE tree. The sweep deletes, and a test
that pointed it at the real state would remove the user's own sessions --
which is why the three session-store directories became module constants
the suite can redirect.
"""

from __future__ import annotations

import os
import stat
import time

import pytest

from delfin.agent import state_paths


_DAY = 86400.0


def _aged(path, days_old: float, *, now: float):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("{}\n", encoding="utf-8")
    stamp = now - days_old * _DAY
    os.utime(path, (stamp, stamp))
    return path


def _mode(p):
    return stat.S_IMODE(os.stat(p).st_mode)


@pytest.fixture(autouse=True)
def _let_the_sweep_run():
    state_paths.reset_maintenance_flag()
    yield
    state_paths.reset_maintenance_flag()


# ---------------------------------------------------------------------------
# The prune itself
# ---------------------------------------------------------------------------

def test_an_old_file_is_removed_and_a_recent_one_is_not(tmp_path):
    now = time.time()
    old = _aged(tmp_path / "traces" / "april.jsonl", 120, now=now)
    new = _aged(tmp_path / "traces" / "today.jsonl", 1, now=now)

    removed = state_paths.prune_old(tmp_path / "traces",
                                    max_age_days=30, now=now)

    assert removed == 1
    assert not old.exists()
    assert new.exists()


def test_retention_can_be_switched_off(tmp_path):
    """Somebody keeping a four-month archive on purpose must be able to."""
    now = time.time()
    old = _aged(tmp_path / "traces" / "april.jsonl", 120, now=now)

    assert state_paths.prune_old(tmp_path / "traces",
                                 max_age_days=0, now=now) == 0
    assert old.exists()


def test_the_prune_leaves_the_directory_in_place(tmp_path):
    """Removing it would race the writer about to use it."""
    now = time.time()
    _aged(tmp_path / "traces" / "old.jsonl", 90, now=now)
    state_paths.prune_old(tmp_path / "traces", max_age_days=30, now=now)
    assert (tmp_path / "traces").is_dir()


def test_the_prune_does_not_follow_a_symlink_out_of_the_tree(tmp_path):
    now = time.time()
    outside = _aged(tmp_path / "precious.txt", 400, now=now)
    root = tmp_path / "traces"
    root.mkdir()
    (root / "link.jsonl").symlink_to(outside)

    state_paths.prune_old(root, max_age_days=30, now=now)

    assert outside.exists(), "the prune deleted through a symlink"


def test_a_glob_bounds_the_prune_to_what_it_owns(tmp_path):
    """The rotated audit logs share ~/.delfin with the settings file and
    the credential store. An unbounded prune there would delete both."""
    now = time.time()
    rotated = _aged(tmp_path / "audit-2026-W20.log", 200, now=now)
    creds = _aged(tmp_path / "credentials.json", 200, now=now)

    state_paths.prune_old(tmp_path, max_age_days=30,
                          glob="audit-*.log", now=now)

    assert not rotated.exists()
    assert creds.exists(), "the prune reached past its glob"


# ---------------------------------------------------------------------------
# The two tiers
# ---------------------------------------------------------------------------

def test_a_resumable_session_outlives_a_tool_trace(tmp_path):
    """Telemetry about finished work and work the user can still resume
    are not the same thing, and a single number would have to be wrong for
    one of them."""
    now = time.time()
    trace = _aged(tmp_path / "traces" / "t.jsonl", 45, now=now)
    session = _aged(tmp_path / "sessions" / "s.json", 45, now=now)

    state_paths.run_startup_maintenance(
        dirs=[state_paths.StateDir(tmp_path / "traces", "derived"),
              state_paths.StateDir(tmp_path / "sessions", "session")],
        retention_days=30, session_retention_days=90, now=now)

    assert not trace.exists()
    assert session.exists()


def test_a_bug_report_is_pruned_but_its_modes_are_left_alone(tmp_path):
    """A filed report is made group-readable ON PURPOSE so a maintainer
    can read it. A sweep that tightened it back would fight the writer and
    leave the mode depending on when the last start-up happened."""
    now = time.time()
    fresh = tmp_path / "agent_bugs" / "20260810_x"
    old = tmp_path / "agent_bugs" / "20260401_y"
    _aged(fresh / "report.json", 2, now=now)
    _aged(old / "report.json", 120, now=now)
    os.chmod(fresh / "report.json", 0o640)

    state_paths.run_startup_maintenance(
        dirs=[state_paths.StateDir(tmp_path / "agent_bugs", "prune")],
        retention_days=30, now=now)

    assert not old.exists(), "the stale report directory was left behind"
    assert (fresh / "report.json").exists()
    assert _mode(fresh / "report.json") == 0o640, "the sweep re-tightened it"


def test_a_keep_target_is_repaired_but_never_pruned(tmp_path):
    """The active audit log and the failure log are the running record and
    both already self-trim."""
    now = time.time()
    log = _aged(tmp_path / "audit.log", 400, now=now)
    os.chmod(log, 0o664)

    state_paths.run_startup_maintenance(
        dirs=[state_paths.StateDir(tmp_path, "keep", "audit.log")],
        retention_days=30, now=now)

    assert log.exists()
    assert _mode(log) == 0o600


# ---------------------------------------------------------------------------
# The sweep at start-up
# ---------------------------------------------------------------------------

def test_the_sweep_prunes_and_repairs_in_one_pass(tmp_path):
    now = time.time()
    old = _aged(tmp_path / "traces" / "old.jsonl", 90, now=now)
    keep = _aged(tmp_path / "traces" / "new.jsonl", 2, now=now)
    os.chmod(keep, 0o664)
    os.chmod(tmp_path / "traces", 0o775)

    result = state_paths.run_startup_maintenance(
        dirs=[state_paths.StateDir(tmp_path / "traces", "derived")],
        retention_days=30, now=now)

    assert result["ran"] is True
    assert result["pruned"] == 1
    assert not old.exists()
    assert _mode(keep) == 0o600
    assert _mode(tmp_path / "traces") == 0o700


def test_the_sweep_runs_once_per_process(tmp_path):
    """It hangs off the first state write, so a second pass on every
    session save would put a tree walk on the hot path."""
    now = time.time()
    _aged(tmp_path / "traces" / "old.jsonl", 90, now=now)
    dirs = [state_paths.StateDir(tmp_path / "traces", "derived")]

    first = state_paths.run_startup_maintenance(dirs=dirs,
                                                retention_days=30, now=now)
    second = state_paths.run_startup_maintenance(dirs=dirs,
                                                 retention_days=30, now=now)

    assert first["ran"] is True
    assert second["ran"] is False
    assert second["pruned"] == 0


def test_the_first_session_write_triggers_the_sweep(tmp_path, monkeypatch):
    """The helper existing with no caller is the bug this closes -- that is
    exactly what happened to cleanup_old_sessions."""
    from delfin.agent import session_store

    now = time.time()
    monkeypatch.setattr(session_store, "_SESSIONS_DIR", tmp_path / "sessions")
    monkeypatch.setattr(session_store, "_transcript_archive_path",
                        lambda: tmp_path / "transcript_archive")
    stale = _aged(tmp_path / "transcript_archive" / "old.jsonl", 200, now=now)
    os.chmod(stale, 0o664)

    session_store._ensure_dir()

    assert not stale.exists(), "the sweep never ran"


def test_the_sweep_never_raises_on_a_directory_it_cannot_read(tmp_path):
    """Start-up maintenance that can abort start-up is worse than the
    permissions it fixes."""
    result = state_paths.run_startup_maintenance(
        dirs=[state_paths.StateDir(tmp_path / "does_not_exist", "derived")],
        retention_days=30)
    assert result["ran"] is True
    assert result["pruned"] == 0


# ---------------------------------------------------------------------------
# The directory table
# ---------------------------------------------------------------------------

def test_every_state_directory_the_suite_redirected_is_swept(tmp_path):
    """A directory missing from the table silently keeps growing, which is
    the failure mode this whole finding is about."""
    names = {mod.rsplit(".", 1)[-1] + "." + attr
             for mod, attr, _k, _g, _f in state_paths._STATE_DIR_SOURCES}
    for expected in (
        "session_store._SESSIONS_DIR",
        "session_store._transcript_archive_path",
        "session_store._handoffs_path",
        "session_store._bundles_path",
        "tool_trace._DIR",
        "turn_metrics._DIR",
        "subagents._SESSIONS_DIR",
        "subagents._RUNNING_DIR",
        "subagents._PENDING_DIR",
        "bug_report._FALLBACK_DIR",
    ):
        assert expected in names, expected


def test_the_sweep_resolves_directories_live_not_at_import(tmp_path,
                                                           monkeypatch):
    """A snapshot taken at import time would point the sweep at the real
    home from inside a test run -- and its prune deletes."""
    from delfin.agent import tool_trace

    monkeypatch.setattr(tool_trace, "_DIR", tmp_path / "redirected")
    (tmp_path / "redirected").mkdir()

    found = {str(d.path) for d in state_paths.state_dirs()}

    assert str(tmp_path / "redirected") in found


def test_no_swept_directory_is_the_users_real_state():
    """The guard that makes the rest of this file safe to run."""
    import pathlib

    real = (pathlib.Path.home() / ".delfin").resolve()
    for entry in state_paths.state_dirs():
        try:
            entry.path.resolve().relative_to(real)
        except ValueError:
            continue
        # ~/.delfin itself is allowed only for glob-bounded, repair-only
        # targets; a recursive prune there would delete the user's
        # settings, credentials and memory.
        assert entry.glob, f"unbounded sweep of real state: {entry.path}"
