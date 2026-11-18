"""The retention sweep can be run by a person, and the service is bounded.

Two findings of the 2026-11-17 state audit, which belong together because
the first is what makes the second matter:

* ``state_paths.run_startup_maintenance`` has done the repair and the
  prune since it was written, and its only caller is
  ``session_store._startup_maintenance`` -- once per process, behind a
  module flag. Nothing a person can invoke calls it, and
  ``cleanup_old_sessions`` has no production caller at all. A retention
  policy nobody can run is a setting, not a policy.
* ``~/.delfin/agent_service/`` was in no prune table. Every turn of every
  session appends to a ``<sid>.jsonl`` under ``sessions/``, and the
  supervisor's whole stdout and stderr go to one ``logs/service.log``.
  The only limits were in memory (``RecentEvents._trim``,
  ``REPLAY_BUDGET_BYTES``), which bound what is replayed, not what is
  kept. A full disk is how the state files the agent depends on become
  unreadable in the first place.

What this file pins, beyond "the numbers come out right":

* the prune's default is the DRY RUN, and the listing comes from
  ``prune_old`` itself rather than from a second implementation of the
  age rule in the caller;
* the ACTIVE ``service.log`` is not prunable -- the running supervisor
  holds it open as its stdout, so only rotated copies are swept;
* rotation happens where the log is opened, and a second rotation does
  not overwrite the first copy.
"""

from __future__ import annotations

import os
import time

from delfin.agent import state_paths as SP
from delfin.agent.service import eventlog, launcher


def _aged(path, days):
    path.write_text("x", encoding="utf-8")
    old = time.time() - days * 86400.0
    os.utime(path, (old, old))
    return path


# --- the dry run ----------------------------------------------------------

def test_a_dry_run_reports_and_deletes_nothing(tmp_path):
    seen = []
    _aged(tmp_path / "old.jsonl", 90)
    _aged(tmp_path / "fresh.jsonl", 1)

    n = SP.prune_old(tmp_path, max_age_days=30, dry_run=True,
                     on_remove=seen.append)

    assert n == 1, "it counts what it would delete"
    assert [q.name for q in seen] == ["old.jsonl"]
    assert (tmp_path / "old.jsonl").exists(), "nothing may be removed"
    assert (tmp_path / "fresh.jsonl").exists()


def test_the_real_prune_removes_exactly_what_the_dry_run_named(tmp_path):
    """One implementation of the age rule, so the two cannot disagree."""
    _aged(tmp_path / "old.jsonl", 90)
    _aged(tmp_path / "fresh.jsonl", 1)
    named = []
    SP.prune_old(tmp_path, max_age_days=30, dry_run=True,
                 on_remove=named.append)

    removed = SP.prune_old(tmp_path, max_age_days=30)

    assert removed == len(named)
    assert not (tmp_path / "old.jsonl").exists()
    assert (tmp_path / "fresh.jsonl").exists()


def test_a_dry_run_leaves_empty_subdirectories_alone(tmp_path):
    """The real prune drops them; a dry run that did would change state."""
    sub = tmp_path / "one"
    sub.mkdir()
    _aged(sub / "old.jsonl", 90)

    SP.prune_old(tmp_path, max_age_days=30, dry_run=True)

    assert sub.is_dir(), "a dry run must not remove a directory either"


def test_the_sweep_still_repairs_in_a_dry_run(tmp_path):
    """Repair destroys nothing, so it is not withheld from a dry run."""
    f = _aged(tmp_path / "loose.jsonl", 1)
    os.chmod(f, 0o644)

    result = SP.run_startup_maintenance(
        dirs=[SP.StateDir(path=tmp_path, kind="derived", glob=None)],
        force=True, dry_run=True, retention_days=30,
        session_retention_days=30)

    assert result["repaired"] >= 1
    assert f.exists()


# --- the service directories are in the tree ------------------------------

def test_the_state_tree_and_the_service_agree_on_where_it_keeps_state(
        tmp_path, monkeypatch):
    """The layout is written twice, so it is pinned once.

    state_paths resolves the service's directories itself instead of
    importing the service layer: the table's resolver runs on the first
    session write of a process, and naming service.eventlog would pull
    those modules into that path -- a coupling the session store did not
    have. The copy is the price, and this is what keeps it honest: both
    answers must agree under every redirection either of them honours.
    """
    from delfin.agent.service import eventlog as el

    for env, value in (
        ("DELFIN_AGENT_SERVICE_DIR", str(tmp_path / "explicit")),
        ("DELFIN_SCRATCH_STATE", str(tmp_path / "scratch")),
    ):
        monkeypatch.delenv("DELFIN_AGENT_SERVICE_DIR", raising=False)
        monkeypatch.delenv("DELFIN_SCRATCH_STATE", raising=False)
        monkeypatch.setenv(env, value)
        assert SP._service_dir() == el.service_dir(), env
        assert SP._service_sessions_dir() == el.sessions_dir(), env
        assert SP._service_logs_dir() == launcher.service_log().parent, env

    # And with neither set, both fall back to the same home.
    monkeypatch.delenv("DELFIN_AGENT_SERVICE_DIR", raising=False)
    monkeypatch.delenv("DELFIN_SCRATCH_STATE", raising=False)
    monkeypatch.setattr("pathlib.Path.home", lambda: tmp_path / "home")
    assert SP._service_dir() == el.service_dir()


def test_the_sweep_does_not_import_the_service_layer(monkeypatch):
    """The resolver runs on the first session write. It must not be the
    thing that loads the service modules."""
    import inspect

    src = inspect.getsource(SP)
    table = src[src.index("_STATE_DIR_SOURCES"):src.index("def state_dirs")]
    assert "delfin.agent.service" not in table, (
        "no table entry may name a module of the service layer")


def test_the_per_session_event_logs_are_a_named_state_directory(
        tmp_path, monkeypatch):
    monkeypatch.setenv(eventlog.SERVICE_DIR_ENV, str(tmp_path / "svc"))
    sessions = eventlog.sessions_dir()
    sessions.mkdir(parents=True)
    _aged(sessions / "s1.jsonl", 90)

    entries = [e for e in SP.state_dirs()
               if e.path.resolve() == sessions.resolve()]
    assert entries, (
        "the per-session event logs must be in the state tree: nothing "
        "else removed them")
    assert entries[0].kind == "derived"
    assert entries[0].glob == "*.jsonl", (
        "bounded by the glob, so the sweep touches nothing else in there")


def test_the_log_path_still_resolves_under_the_named_directory(
        tmp_path, monkeypatch):
    monkeypatch.setenv(eventlog.SERVICE_DIR_ENV, str(tmp_path / "svc"))
    assert eventlog.log_path("abc").parent == eventlog.sessions_dir()


def test_only_rotated_service_logs_are_prunable(tmp_path, monkeypatch):
    """The active log is the running supervisor's stdout."""
    monkeypatch.setenv(eventlog.SERVICE_DIR_ENV, str(tmp_path / "svc"))
    logs = launcher.service_log().parent
    logs.mkdir(parents=True)
    _aged(launcher.service_log(), 90)
    _aged(logs / "service.log.20260101-000000", 90)

    entries = [e for e in SP.state_dirs()
               if e.path.resolve() == logs.resolve()]
    assert entries, "the service log directory must be in the tree"
    assert entries[0].glob == "service.log.*", (
        "the ACTIVE service.log must not match: pruning it would delete "
        "the log of a live service")

    removed = SP.prune_old(logs, max_age_days=30, glob=entries[0].glob)
    assert removed == 1
    assert launcher.service_log().exists(), "the active log survives"


# --- rotation -------------------------------------------------------------

def test_a_small_log_is_not_rotated(tmp_path):
    log = tmp_path / "service.log"
    log.write_bytes(b"k" * 1000)
    assert launcher._rotate_service_log(log) is None
    assert log.exists()


def test_an_oversized_log_is_moved_aside_whole(tmp_path):
    log = tmp_path / "service.log"
    log.write_bytes(b"g" * (launcher.SERVICE_LOG_MAX_BYTES + 1))

    aside = launcher._rotate_service_log(log)

    assert aside is not None and aside.name.startswith("service.log.")
    assert aside.stat().st_size == launcher.SERVICE_LOG_MAX_BYTES + 1, (
        "a rename, so nothing written so far is lost")
    assert not log.exists(), "the next open recreates it"


def test_a_second_rotation_keeps_the_first_copy(tmp_path):
    """os.replace overwrites its destination without a word."""
    log = tmp_path / "service.log"
    for marker in (b"a", b"b"):
        log.write_bytes(marker * (launcher.SERVICE_LOG_MAX_BYTES + 1))
        launcher._rotate_service_log(log)

    kept = sorted(q for q in tmp_path.iterdir()
                  if q.name.startswith("service.log."))
    assert len(kept) == 2, "both rotated logs are kept"
    assert {q.read_bytes()[:1] for q in kept} == {b"a", b"b"}


def test_rotation_happens_where_the_log_is_opened():
    """Not a style point: the descriptor is handed to the supervisor as its
    stdout, and a rename follows the inode. Rotating under a running
    service would leave it writing into the rotated file while the live
    log stayed empty."""
    import inspect

    src = inspect.getsource(launcher.start_service)
    rotate = src.index("_rotate_service_log(")
    open_fd = src.index("os.open(log")
    assert rotate < open_fd, (
        "the rotation must precede the open of the log it rotates")
