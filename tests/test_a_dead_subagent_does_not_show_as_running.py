"""A subagent whose process was killed showed as running forever.

The live registry is one file per subagent under
``~/.delfin/subagent_running/``. It is created when the run starts and
removed in a ``finally``. That covers every way a RUN can end -- and none
of the ways the PROCESS can end.

Kill the dashboard, lose the ssh session, let the machine reboot mid-run:
the ``finally`` never executes, the file stays, and ``read_running()``
simply globs whatever files are there. The panel keeps showing the
subagent as working, and ``get_subagent_result`` keeps answering
``status: "running"`` -- so a parent agent polling a backgrounded
delegate waits for a report that can no longer arrive. Nothing ages the
entry out, so it survives every later session too.

The mechanism to fix it already exists a module away: ``bash_jobs``
records the owning pid plus its start ticks and checks both, the ticks
being the guard against a recycled pid pointing at an unrelated process.
The same stamp is applied here rather than reinvented.

Entries written before the field existed are treated as alive -- an
unknown owner is not evidence of a dead one, and reaping them would
delete live work on the first run after an upgrade.
"""

from __future__ import annotations

import json
import os

import pytest

from delfin.agent import subagents as sa


@pytest.fixture
def running_dir(tmp_path, monkeypatch):
    d = tmp_path / "subagent_running"
    d.mkdir()
    monkeypatch.setattr(sa, "_RUNNING_DIR", d)
    return d


def _write(d, sa_id: str, **extra) -> None:
    entry = {"type": "explore", "description": "d", "started_at": 0.0,
             "actions": [], "last_action": "", "transcript": []}
    entry.update(extra)
    (d / f"{sa_id}.json").write_text(json.dumps(entry), encoding="utf-8")


def _dead_pid() -> int:
    """A pid that is certainly not running: fork a child and reap it."""
    pid = os.fork()
    if pid == 0:                       # pragma: no cover - child
        os._exit(0)
    os.waitpid(pid, 0)
    return pid


# ---------------------------------------------------------------------------
# A dead owner means a dead subagent
# ---------------------------------------------------------------------------

def test_an_entry_whose_owner_is_gone_is_not_running(running_dir):
    _write(running_dir, "gone", owner_pid=_dead_pid())
    assert "gone" not in sa.read_running()


def test_an_entry_owned_by_this_process_is_running(running_dir):
    _write(running_dir, "live", **sa._owner_stamp())
    assert "live" in sa.read_running()


def test_a_recycled_pid_does_not_resurrect_it(running_dir):
    """Our own pid, but stamped with a start time that is not ours."""
    _write(running_dir, "recycled", owner_pid=os.getpid(),
           owner_start_ticks=1)
    assert "recycled" not in sa.read_running()


def test_an_entry_without_an_owner_still_counts_as_running(running_dir):
    """Written before the field existed. An unknown owner is not evidence
    of a dead one, and reaping it would delete live work after an upgrade."""
    _write(running_dir, "legacy")
    assert "legacy" in sa.read_running()


# ---------------------------------------------------------------------------
# The caller waiting for the report is told
# ---------------------------------------------------------------------------

def test_polling_a_killed_subagent_says_it_died(running_dir, monkeypatch):
    monkeypatch.setattr(sa, "load_subagent_session", lambda *_a, **_k: None)
    _write(running_dir, "killed", owner_pid=_dead_pid())
    out = sa.get_subagent_result("killed")
    assert out["status"] == "died"
    assert "no longer running" in out["error"]


def test_the_report_is_still_returned_if_the_session_was_written(
    running_dir, monkeypatch,
):
    """The process died AFTER the session was persisted -- the work is not
    lost, and reporting "died" would throw away a usable report."""
    monkeypatch.setattr(
        sa, "load_subagent_session",
        lambda *_a, **_k: {"messages": [
            {"role": "assistant", "content": "the report"}]})
    _write(running_dir, "late", owner_pid=_dead_pid())
    out = sa.get_subagent_result("late")
    assert out["status"] == "finished"
    assert out["final_text"] == "the report"


def test_an_unknown_id_is_still_unknown(running_dir, monkeypatch):
    monkeypatch.setattr(sa, "load_subagent_session", lambda *_a, **_k: None)
    assert sa.get_subagent_result("never-existed")["status"] == "unknown"


def test_a_live_subagent_still_reports_running(running_dir):
    _write(running_dir, "busy", **sa._owner_stamp())
    assert sa.get_subagent_result("busy")["status"] == "running"


# ---------------------------------------------------------------------------
# ...and the file does not pile up forever
# ---------------------------------------------------------------------------

def test_the_stale_file_is_reaped(running_dir):
    _write(running_dir, "gone", owner_pid=_dead_pid())
    _write(running_dir, "live", **sa._owner_stamp())
    assert sa.reap_dead_running() == ["gone"]
    assert not (running_dir / "gone.json").exists()
    assert (running_dir / "live.json").exists()


def test_reaping_an_empty_registry_is_harmless(running_dir):
    assert sa.reap_dead_running() == []


def test_a_new_run_reaps_what_the_last_crash_left(running_dir):
    """The reap has to be triggered by something; starting a run is the
    moment that is guaranteed to happen and cheap enough to carry it."""
    _write(running_dir, "gone", owner_pid=_dead_pid())
    sa._running_update("fresh", {"type": "explore", "description": "d",
                                 "started_at": 0.0})
    assert not (running_dir / "gone.json").exists()


def test_a_started_entry_carries_its_owner(running_dir):
    sa._running_update("fresh", {"type": "explore", "description": "d",
                                 "started_at": 0.0})
    ent = json.loads((running_dir / "fresh.json").read_text(encoding="utf-8"))
    assert ent["owner_pid"] == os.getpid()


def test_removing_an_entry_still_works(running_dir):
    sa._running_update("x", {"type": "explore"})
    sa._running_update("x", None)
    assert not (running_dir / "x.json").exists()
