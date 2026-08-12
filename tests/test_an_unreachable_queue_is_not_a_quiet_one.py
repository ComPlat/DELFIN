"""Every SLURM state query failed open into silence, and into the wrong list.

FAILING OPEN. The query returned ``""`` on any non-zero exit or exception,
so "the scheduler says nothing changed" and "the scheduler could not be
asked" arrived as the same answer. A watch list could then sit for days
reporting nothing while the queue it was watching was unreachable, and the
only exit was a 7-day age prune that popped the entry without a word.

ONE BAD ID DROPPED THE BATCH. ``squeue -j`` exits non-zero as soon as ONE
of the ids has aged out of the queue, and the whole batch went with it: a
watch list of ten jobs stopped reporting because the eleventh had finished.

TWO WATCH LISTS, ONE WATCHER. The daemon polls
``~/.delfin/watched_jobs.json``; everything the agent registers goes to
``<workspace>/.delfin/agent_watched_jobs.json``, which is read only from
inside a turn. So once the session ended, a twelve-hour job the agent had
submitted itself was watched by nobody -- while the tool result that
registered it said, in so many words, that its completion would be
reported without polling.

Every scheduler call here is injected. No real queue is touched.
"""

from __future__ import annotations

import json
import time

import pytest

from delfin.agent import job_monitor as JM


@pytest.fixture(autouse=True)
def _isolated(tmp_path, monkeypatch):
    monkeypatch.setattr(JM, "_AGENT_WATCH_INDEX_PATH",
                        tmp_path / "agent_watch_index.json")
    monkeypatch.setattr(JM, "_WATCHED_PATH", tmp_path / "watched_jobs.json")


def _run(*, squeue=None, sacct=None, squeue_fails=False):
    """An injected scheduler. ``None`` from a command means "could not run"."""
    def _fn(cmd):
        if cmd[0] == "squeue":
            if squeue_fails:
                return None
            return squeue if squeue is not None else ""
        if cmd[0] == "sacct":
            return sacct
        return ""
    return _fn


# ---------------------------------------------------------------------------
# Tri-state: known / absent / unavailable
# ---------------------------------------------------------------------------

def test_a_state_the_scheduler_gave_us_is_reported_as_itself():
    states = JM.query_job_states_detailed(
        ["111"], _run(squeue="111 RUNNING\n"))

    assert states == {"111": "RUNNING"}


def test_an_id_the_scheduler_does_not_know_is_absent_not_unavailable():
    states = JM.query_job_states_detailed(
        ["111"], _run(squeue="", sacct=""))

    assert states == {"111": ""}


def test_a_scheduler_that_cannot_be_asked_says_so():
    """This is the whole point: silence and 'unchanged' are different facts."""
    states = JM.query_job_states_detailed(
        ["111"], _run(squeue_fails=True, sacct=None))

    assert states == {"111": JM.STATE_UNAVAILABLE}


def test_one_aged_out_id_does_not_drop_the_whole_batch():
    """squeue -j fails as a unit; the rest of the list is still knowable."""
    def _fn(cmd):
        if cmd[0] == "squeue":
            ids = cmd[2].split(",")
            if "999" in ids and len(ids) > 1:
                return None             # squeue's all-or-nothing failure
            if ids == ["999"]:
                return None
            return "".join(f"{i} RUNNING\n" for i in ids)
        return "999 COMPLETED\n"

    states = JM.query_job_states_detailed(["111", "222", "999"], _fn)

    assert states["111"] == "RUNNING"
    assert states["222"] == "RUNNING"
    assert states["999"] == "COMPLETED"


def test_the_narrow_view_still_only_returns_placed_jobs():
    """The old signature keeps its old meaning for callers that branch on it."""
    states = JM.query_job_states(
        ["111", "222"], _run(squeue="111 RUNNING\n", sacct=""))

    assert states == {"111": "RUNNING"}


# ---------------------------------------------------------------------------
# The degradation is visible where it matters
# ---------------------------------------------------------------------------

def test_a_watched_job_whose_queue_is_unreachable_is_reported_as_unknown(
        tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    JM.register_agent_job(ws, "555", "twelve-hour opt")

    done = JM.check_agent_jobs(
        ws, _run(squeue_fails=True, sacct=None))

    assert len(done) == 1
    assert done[0]["state"] == JM.STATE_UNAVAILABLE
    assert done[0]["ok"] is False
    assert "unknown, not unchanged" in done[0]["degraded"]


def test_the_entry_is_kept_while_the_queue_is_unreachable(tmp_path):
    """Dropping it would turn a temporary outage into a permanent loss."""
    ws = tmp_path / "ws"
    ws.mkdir()
    JM.register_agent_job(ws, "555", "twelve-hour opt")

    JM.check_agent_jobs(ws, _run(squeue_fails=True, sacct=None))
    data = json.loads(JM._agent_watch_path(ws).read_text())

    assert "555" in data["jobs"]
    assert data["jobs"]["555"]["unavailable_since"] > 0


def test_the_job_is_reported_normally_once_the_queue_answers_again(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    JM.register_agent_job(ws, "555", "twelve-hour opt")
    JM.check_agent_jobs(ws, _run(squeue_fails=True, sacct=None))

    done = JM.check_agent_jobs(ws, _run(squeue="", sacct="555 COMPLETED\n"))

    assert [d["state"] for d in done] == ["COMPLETED"]
    assert done[0]["ok"] is True
    assert JM.check_agent_jobs(ws, _run(squeue="", sacct="")) == []


def test_giving_up_on_an_unresolved_job_is_announced(tmp_path, monkeypatch):
    """The 7-day prune used to pop the entry silently. Nobody was told that
    the thing they were promised a report on would never get one."""
    emitted: list[tuple] = []
    monkeypatch.setattr(
        "delfin.agent.attention.emit_attention",
        lambda kind, **kw: emitted.append((kind, kw.get("title", ""))))
    ws = tmp_path / "ws"
    ws.mkdir()
    JM.register_agent_job(ws, "555", "twelve-hour opt")
    data = JM.load_watched(JM._agent_watch_path(ws))
    data["jobs"]["555"]["added_at"] = time.time() - 8 * 24 * 3600
    JM.save_watched(data, JM._agent_watch_path(ws))

    assert JM.check_agent_jobs(ws, _run(squeue="", sacct="")) == []

    assert emitted, "a job we stopped watching without an answer must be said"
    assert emitted[0][0] == "run_failed"
    assert "555" in emitted[0][1]


@pytest.mark.parametrize("resolution", [
    {"last_state": "COMPLETED"},
    {"last_state": "OUT_OF_MEMORY"},
    {"daemon_notified": True},
])
def test_a_job_that_did_get_an_answer_is_pruned_without_an_alarm(
        tmp_path, monkeypatch, resolution):
    """Only the UNRESOLVED prune is worth an interruption."""
    emitted: list = []
    monkeypatch.setattr(
        "delfin.agent.attention.emit_attention",
        lambda kind, **kw: emitted.append(kind))
    ws = tmp_path / "ws"
    ws.mkdir()
    JM.register_agent_job(ws, "555", "done long ago")
    path = JM._agent_watch_path(ws)
    data = JM.load_watched(path)
    data["jobs"]["555"]["added_at"] = time.time() - 8 * 24 * 3600
    data["jobs"]["555"].update(resolution)
    JM.save_watched(data, path)

    JM.check_agent_jobs(ws, _run(squeue="", sacct=""))

    assert emitted == []


# ---------------------------------------------------------------------------
# The daemon reads the agent's own watch lists too
# ---------------------------------------------------------------------------

def test_registering_a_job_makes_its_workspace_findable(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    JM.register_agent_job(ws, "555", "opt")

    assert JM.agent_watch_workspaces() == [str(ws)]


def test_the_sweep_covers_every_workspace_a_job_was_registered_in(tmp_path):
    first, second = tmp_path / "a", tmp_path / "b"
    for folder in (first, second):
        folder.mkdir()
    JM.register_agent_job(first, "111", "one")
    JM.register_agent_job(second, "222", "two")

    done = JM.check_all_agent_jobs(
        _run(squeue="", sacct="111 COMPLETED\n222 FAILED\n"))

    assert sorted(d["job_id"] for d in done) == ["111", "222"]
    assert {d["workspace"] for d in done} == {str(first), str(second)}


def test_the_daemon_does_not_eat_the_event_the_next_turn_is_owed(tmp_path):
    """Two consumers, two flags: the daemon notifies, the turn consumes."""
    ws = tmp_path / "ws"
    ws.mkdir()
    JM.register_agent_job(ws, "555", "opt")
    run = _run(squeue="", sacct="555 COMPLETED\n")

    assert [d["job_id"] for d in JM.check_all_agent_jobs(run)] == ["555"]
    # The daemon does not report it twice...
    assert JM.check_all_agent_jobs(run) == []
    # ...and the in-turn drain still gets it exactly once.
    assert [d["job_id"] for d in JM.check_agent_jobs(ws, run)] == ["555"]
    assert JM.check_agent_jobs(ws, run) == []


def test_the_daemon_notifies_about_a_job_no_session_is_open_for(
        tmp_path, monkeypatch):
    emitted: list[tuple] = []
    monkeypatch.setattr(
        "delfin.agent.attention.emit_attention",
        lambda kind, **kw: emitted.append((kind, kw.get("title", ""))))
    monkeypatch.setattr(JM, "_default_run",
                        _run(squeue="", sacct="555 COMPLETED\n"))
    monkeypatch.setattr(JM, "check_all_agent_jobs",
                        lambda *a, **kw: [{
                            "job_id": "555", "kind": "slurm",
                            "description": "twelve-hour opt",
                            "state": "COMPLETED", "ok": True,
                            "exit_code": None, "signatures": [],
                            "workspace": str(tmp_path)}])

    reported = JM.check_agent_workspaces_once({})

    assert [r["job_id"] for r in reported] == ["555"]
    assert emitted == [("run_finished", "Job 555 finished")]


def test_a_failure_found_by_the_daemon_goes_through_the_usual_pipeline(
        tmp_path, monkeypatch):
    """Findings are recorded and announced exactly as for the shared list."""
    recorded: list = []
    announced: list = []
    monkeypatch.setattr(JM, "record_finding",
                        lambda finding, path=None: recorded.append(finding))
    monkeypatch.setattr(JM, "announce",
                        lambda finding, settings=None: announced.append(finding))
    monkeypatch.setattr(JM, "diagnose_finding",
                        lambda finding, **kw: finding)
    monkeypatch.setattr(JM, "check_all_agent_jobs",
                        lambda *a, **kw: [{
                            "job_id": "666", "kind": "slurm",
                            "description": "opt", "state": "OUT_OF_MEMORY",
                            "ok": False, "exit_code": None,
                            "signatures": ["out-of-memory"],
                            "workspace": str(tmp_path)}])

    JM.check_agent_workspaces_once({})

    assert [f.job_id for f in recorded] == ["666"]
    assert [f.state for f in announced] == ["OUT_OF_MEMORY"]
    assert recorded[0].signatures == ["out-of-memory"]


def test_the_daemon_loop_polls_both_lists(monkeypatch):
    seen: list[str] = []
    monkeypatch.setattr(JM, "monitor_settings",
                        lambda settings=None: {
                            "enabled": True, "interval_s": 1,
                            "auto_diagnose": False, "webhook_url": "",
                            "provider": "", "model": "", "backend": ""})
    monkeypatch.setattr(JM, "check_once", lambda: seen.append("shared") or [])
    monkeypatch.setattr(JM, "check_agent_workspaces_once",
                        lambda settings=None: seen.append("agent") or [])

    JM.run_loop(max_iterations=1, sleep_fn=lambda _s: None)

    assert seen == ["shared", "agent"]


def test_an_unreachable_queue_does_not_overwrite_the_last_known_state(
        tmp_path):
    """check_once used to record "" over a real state and lose it."""
    path = tmp_path / "watched.json"
    JM.add_watch("777", folder=str(tmp_path), path=path)
    JM.check_once(path, _run(squeue="777 RUNNING\n"))
    assert JM.load_watched(path)["jobs"]["777"]["last_state"] == "RUNNING"

    JM.check_once(path, _run(squeue_fails=True, sacct=None))

    entry = JM.load_watched(path)["jobs"]["777"]
    assert entry["last_state"] == "RUNNING"
    assert entry["unavailable_since"] > 0
