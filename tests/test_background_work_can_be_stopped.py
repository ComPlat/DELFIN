"""Background work can be stopped from the dashboard, one × per row.

Seen live on 2026-09-15: a new session's Background panel listed a CI watch
pending for four hours and a wake-up waiting for that CI result -- both left
by an older session, and neither could be ended from the screen.
"""

from __future__ import annotations

import time

import pytest

from delfin.agent import background_view as bgv
from delfin.agent import bash_jobs as bj
from delfin.agent import job_monitor as jm
from delfin.agent import scheduler as S
from delfin.agent import subagents as sa


@pytest.fixture
def ws(tmp_path, monkeypatch):
    monkeypatch.setattr(jm, "_AGENT_WATCH_INDEX_PATH", tmp_path / "index.json")
    monkeypatch.setattr(S, "_GLOBAL", S.Scheduler(path=tmp_path / "cron.json"))
    return tmp_path


def test_a_watch_is_stopped_and_the_job_left_alone(ws):
    jm.register_agent_job(ws, "ci:ComPlat/DELFIN@8feb1765", "CI for 8feb1765",
                          extra={"session_id": "S"})
    message = bgv.cancel(ws, "watches", "ci:ComPlat/DELFIN@8feb1765")
    assert "keeps running" in message
    assert jm.load_watched(jm._agent_watch_path(ws))["jobs"] == {}


def test_a_wake_up_is_cancelled(ws):
    ent = S.get_scheduler().schedule_once(
        delay_seconds=600, prompt="check CI", session_id="S")
    assert "cancelled" in bgv.cancel(ws, "wakeups", ent.id)
    assert S.get_scheduler().list_entries() == []


def test_a_shell_is_terminated(ws):
    job = bj.get_registry().start(command="sleep 60", cwd=str(ws),
                                  workspace=str(ws), session_id="S")
    bgv.cancel(ws, "shells", job.job_id)
    deadline = time.time() + 10
    while job.poll() is None and time.time() < deadline:
        time.sleep(0.1)
    assert job.poll() is not None


def test_a_background_agent_hears_its_own_stop(monkeypatch):
    monkeypatch.setattr(sa, "read_running", lambda **_: {"a1": {"type": "explore"}})
    assert sa.cancel_background("a1") is True
    assert sa.cancel_requested("a1") is True
    assert sa.cancel_background("nope") is False
    sa._CANCELLED.discard("a1")


def test_a_stopped_run_sees_the_stop_in_its_client(monkeypatch):
    seen = []

    class _Client:
        model, _provider, _base_url = "m", "", ""
        _permissions = None

        def set_permissions(self, perms):
            self._permissions = perms

        def stream_message(self, messages, system, max_tokens):
            seen.append(self.should_stop())
            return iter(())

    sa._CANCELLED.add("run1")
    try:
        sa.run_subagent(subagent_type="explore", description="probe",
                        prompt="a self-contained briefing of some length",
                        parent_client=_Client(), parent_perms=None,
                        sa_id="run1")
    finally:
        sa._CANCELLED.discard("run1")
    assert seen and seen[0] is True


def test_the_panel_offers_a_stop_for_this_sessions_work_only(tmp_path,
                                                             monkeypatch):
    pytest.importorskip("ipywidgets")
    import ipywidgets as widgets

    from delfin.dashboard import tab_agent
    from delfin.dashboard.context import DashboardContext

    monkeypatch.setattr(jm, "_AGENT_WATCH_INDEX_PATH", tmp_path / "index.json")
    monkeypatch.setattr(S, "_GLOBAL", S.Scheduler(path=tmp_path / "cron.json"))
    for name in ("calc", "archive", "office", "project"):
        (tmp_path / name).mkdir()
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda script: None
    ctx.agent_workspace = str(tmp_path / "project")
    _tab, refs = tab_agent.create_tab(ctx)
    project = tmp_path / "project"

    def _stops():
        out = []
        stack = [_tab]
        while stack:
            node = stack.pop()
            if isinstance(node, widgets.Button) and node.description == "×":
                out.append(node)
            stack.extend(getattr(node, "children", ()) or ())
        return out

    jm.register_agent_job(project, "ci:ComPlat/DELFIN@8feb1765",
                          "left by an older session")
    refs["refresh_background"]()
    assert _stops() == [], "a new session shows nothing it does not own"

    refs["state"]["active_session_id"] = "S"
    jm.register_agent_job(project, "ci:ComPlat/DELFIN@269429eb", "mine",
                          extra={"session_id": "S"})
    refs["refresh_background"]()
    (stop,) = _stops()
    stop.click()
    assert list(jm.load_watched(jm._agent_watch_path(project))["jobs"]) == [
        "ci:ComPlat/DELFIN@8feb1765"]
    assert "Stopped watching" in str(refs["state"]["chat_messages"])
    refs["shutdown"]()
