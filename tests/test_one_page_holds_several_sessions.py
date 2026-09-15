"""One page holds several agent sessions.

Each session is its own agent tab -- its own engine, settings, chat and bug
report -- and the session list shows one at a time. That only works if a
tab keeps to itself: its own working directory, page scripts that act on
the session on screen, no second copy of a conversation, and background
work that stops when the session is closed.
"""

from __future__ import annotations

import inspect

import pytest

pytest.importorskip("ipywidgets")

from delfin.agent import scheduler as S                     # noqa: E402
from delfin.dashboard import tab_agent                      # noqa: E402
from delfin.dashboard.context import DashboardContext       # noqa: E402


@pytest.fixture(autouse=True)
def _scheduler(tmp_path, monkeypatch):
    monkeypatch.setattr(S, "_GLOBAL", S.Scheduler(path=tmp_path / "cron.json"))


def _ctx(tmp_path, **own):
    for name in ("calc", "archive", "office"):
        (tmp_path / name).mkdir(exist_ok=True)
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda script: None
    for key, value in own.items():
        setattr(ctx, key, value)
    return ctx


def test_a_session_works_in_the_directory_it_was_given(tmp_path):
    project = tmp_path / "project"
    project.mkdir()
    _tab, refs = tab_agent.create_tab(_ctx(tmp_path, agent_workspace=str(project)))
    assert refs["workspace"]() == str(project)
    refs["shutdown"]()


def test_the_page_scripts_act_on_the_session_on_screen(tmp_path):
    """document.querySelector returns the first match on the page, so Enter
    in the second session clicked the first session's Send."""
    ctx = _ctx(tmp_path)
    _tab, refs = tab_agent.create_tab(ctx)
    refs["shutdown"]()
    scripts = "\n".join(ctx.init_js_parts)
    assert "window.__delfinQ = window.__delfinQ ||" in scripts
    source = inspect.getsource(tab_agent)
    assert "document.querySelector('.delfin-agent-" not in source


def test_a_session_the_list_opens_does_not_reopen_the_last_one(tmp_path,
                                                                monkeypatch):
    monkeypatch.setenv("DELFIN_RESUME_SESSION", "latest")
    asked = []
    monkeypatch.setattr("delfin.agent.session_store.resume_latest",
                        lambda *a, **kw: asked.append(1) or None)
    _tab, refs = tab_agent.create_tab(_ctx(tmp_path, sessions_managed=True))
    refs["shutdown"]()
    assert asked == []


def test_a_conversation_open_in_another_session_is_not_loaded_twice(tmp_path):
    ctx = _ctx(tmp_path, session_open_elsewhere=lambda sid: sid == "taken")
    _tab, refs = tab_agent.create_tab(ctx)
    refs["load_session"]("taken")
    refs["shutdown"]()
    assert refs["state"]["active_session_id"] == ""
    assert "already open in another session" in str(refs["state"]["chat_messages"])


def test_a_closed_session_stops_its_background_work(tmp_path):
    _tab, refs = tab_agent.create_tab(_ctx(tmp_path))
    state = refs["state"]
    sch = S.get_scheduler()
    sch.add_fire_listener(id(state), lambda entry: None, lambda owner: False)

    refs["shutdown"]()

    assert state["_closed"] is True and state["_subagent_live_stop"] is True
    assert all(key != id(state) for key, _cb, _owns in sch._listeners)
