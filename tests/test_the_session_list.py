"""The agent tab lists its sessions on the left and shows one on the right.

Several sessions run side by side, each a whole agent tab with its own
settings. Each works in the directory it was started in, the one on screen
is the one the rest of the dashboard sees, and the sessions that were open
come back when the dashboard starts again.
"""

from __future__ import annotations

import pytest

pytest.importorskip("ipywidgets")

import ipywidgets as widgets                                # noqa: E402

from delfin.agent import scheduler as S                     # noqa: E402
from delfin.dashboard import agent_sessions as AS           # noqa: E402
from delfin.dashboard.context import DashboardContext       # noqa: E402


class _Build:
    """Stands in for tab_agent.create_tab: records what each session got."""

    def __init__(self):
        self.contexts = []
        self.closed = []

    def __call__(self, ctx):
        state = {"active_session_id": ctx.initial_session_id or "",
                 "chat_messages": [], "streaming": False, "engine": None}
        ctx.agent_state = state
        self.contexts.append(ctx)
        return widgets.VBox(), {
            "state": state,
            "shutdown": lambda: self.closed.append(ctx),
            "workspace": lambda: ctx.agent_workspace,
        }


@pytest.fixture
def home(tmp_path, monkeypatch):
    monkeypatch.setattr(AS, "_OPEN_SESSIONS_PATH", tmp_path / "open.json")
    monkeypatch.setattr(S, "_GLOBAL", S.Scheduler(path=tmp_path / "cron.json"))
    monkeypatch.delenv("DELFIN_RESUME_SESSION", raising=False)
    for name in ("launch", "calc", "agent_workspace", "archive", "office"):
        (tmp_path / name).mkdir()
    monkeypatch.setenv("DELFIN_LAUNCH_CWD", str(tmp_path / "launch"))
    saved = {}
    monkeypatch.setattr(AS, "_saved_session", lambda sid: saved.get(sid))
    return tmp_path, saved


def _ctx(tmp_path):
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           agent_dir=tmp_path / "agent_workspace",
                           archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda script: None
    return ctx


def test_the_first_start_opens_one_session_where_delfin_was_launched(home):
    tmp, _saved = home
    build = _Build()
    _widget, refs = AS.create_tab(_ctx(tmp), build=build)
    assert len(refs["sessions"]()) == 1
    assert build.contexts[0].agent_workspace == str(tmp / "launch")


def test_each_session_works_where_it_was_started(home):
    tmp, _saved = home
    build, ctx = _Build(), _ctx(tmp)
    _widget, refs = AS.create_tab(ctx, build=build)
    refs["open"](str(tmp / "calc"))

    first, second = refs["sessions"]()
    assert [c.agent_workspace for c in build.contexts] == [
        str(tmp / "launch"), str(tmp / "calc")]
    assert (first["tab"].layout.display, second["tab"].layout.display) == ("none", "")
    assert ctx.agent_state is second["refs"]["state"]


def test_the_session_on_screen_is_the_one_the_dashboard_sees(home):
    tmp, _saved = home
    build, ctx = _Build(), _ctx(tmp)
    _widget, refs = AS.create_tab(ctx, build=build)
    refs["open"](str(tmp / "calc"))
    first, second = refs["sessions"]()

    refs["activate"](first["key"])
    assert (first["tab"].layout.display, second["tab"].layout.display) == ("", "none")
    assert ctx.agent_state is first["refs"]["state"]

    second["ctx"].agent_status_html.value = "working in calc"
    assert ctx.agent_status_html.value != "working in calc"
    first["ctx"].agent_status_html.value = "working in launch"
    assert ctx.agent_status_html.value == "working in launch"


def test_the_sessions_that_were_open_come_back(home):
    tmp, saved = home
    build = _Build()
    _widget, refs = AS.create_tab(_ctx(tmp), build=build)
    refs["open"](str(tmp / "calc"))
    for rec, sid in zip(refs["sessions"](), ("s-launch", "s-calc")):
        rec["refs"]["state"]["active_session_id"] = sid
        saved[sid] = {"session_id": sid}
    refs["refresh"]()

    again = _Build()
    _widget, refs2 = AS.create_tab(_ctx(tmp), build=again)
    assert [c.initial_session_id for c in again.contexts] == ["s-launch", "s-calc"]
    assert [c.agent_workspace for c in again.contexts] == [
        str(tmp / "launch"), str(tmp / "calc")]


def test_a_closed_session_stops_and_is_not_reopened(home):
    tmp, saved = home
    build = _Build()
    _widget, refs = AS.create_tab(_ctx(tmp), build=build)
    second = refs["open"](str(tmp / "calc"))
    second["refs"]["state"]["active_session_id"] = "s-calc"
    saved["s-calc"] = {"session_id": "s-calc"}
    refs["refresh"]()

    refs["close"](second["key"])
    assert build.closed == [second["ctx"]]
    assert AS.load_open_sessions() == []

    (last,) = refs["sessions"]()
    refs["close"](last["key"])
    assert len(refs["sessions"]()) == 1, "the list is never left empty"


def test_one_conversation_is_open_in_one_session(home):
    tmp, _saved = home
    build = _Build()
    _widget, refs = AS.create_tab(_ctx(tmp), build=build)
    first = refs["sessions"]()[0]
    first["refs"]["state"]["active_session_id"] = "s1"
    second = refs["open"](str(tmp / "calc"))

    assert second["ctx"].session_open_elsewhere("s1") is True
    assert first["ctx"].session_open_elsewhere("s1") is False
    assert refs["open"]("", "s1") is first
    assert len(build.contexts) == 2


def test_only_the_first_session_adds_the_page_scripts(home):
    tmp, _saved = home
    build, ctx = _Build(), _ctx(tmp)
    _widget, refs = AS.create_tab(ctx, build=build)
    refs["open"](str(tmp / "calc"))
    build.contexts[0].add_init_js("first();")
    build.contexts[1].add_init_js("second();")
    assert "first();" in ctx.init_js_parts
    assert "second();" not in ctx.init_js_parts


def test_two_real_agent_tabs_live_side_by_side(home):
    tmp, _saved = home
    ctx = _ctx(tmp)
    _widget, refs = AS.create_tab(ctx)
    refs["open"](str(tmp / "calc"))
    first, second = refs["sessions"]()

    assert first["refs"]["workspace"]() == str(tmp / "launch")
    assert second["refs"]["workspace"]() == str(tmp / "calc")
    assert first["refs"]["state"] is not second["refs"]["state"]
    for rec in refs["sessions"]():
        refs["close"](rec["key"])
