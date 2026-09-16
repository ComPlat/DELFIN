"""The emergency stop is one button on the session list, and it names what
is open elsewhere.

`delfin-agent stop-all` needs a terminal; the person who has just lost track
of a session on another login node is looking at a dashboard. Two clicks,
so a stray one ends nothing.
"""

from __future__ import annotations

import ast
import inspect

import pytest

pytest.importorskip("ipywidgets")

import ipywidgets as widgets                                # noqa: E402

from delfin.agent import scheduler as S                     # noqa: E402
from delfin.dashboard import agent_sessions as AS           # noqa: E402
from delfin.dashboard.context import DashboardContext       # noqa: E402


def _build(ctx):
    state = {"active_session_id": "", "chat_messages": [], "streaming": False,
             "engine": None}
    ctx.agent_state = state
    return widgets.VBox(), {"state": state, "shutdown": lambda: None,
                            "workspace": lambda: ctx.agent_workspace}


@pytest.fixture
def sidebar(tmp_path, monkeypatch):
    monkeypatch.setattr(AS, "_OPEN_SESSIONS_PATH", tmp_path / "open.json")
    monkeypatch.setattr(S, "_GLOBAL", S.Scheduler(path=tmp_path / "cron.json"))
    monkeypatch.delenv("DELFIN_RESUME_SESSION", raising=False)
    for name in ("launch", "calc", "agent_workspace", "archive", "office"):
        (tmp_path / name).mkdir()
    monkeypatch.setenv("DELFIN_LAUNCH_CWD", str(tmp_path / "launch"))
    monkeypatch.setattr(AS, "_saved_session", lambda sid: None)
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           agent_dir=tmp_path / "agent_workspace",
                           archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda script: None
    widget, _refs = AS.create_tab(ctx, build=_build)
    return widget.children[0]


def _with_class(box, name):
    return next(w for w in box.children if name in w._dom_classes)


def test_one_click_arms_the_stop_and_the_second_gives_it(sidebar, monkeypatch):
    given = []
    monkeypatch.setattr(AS, "_give_emergency_stop", lambda: given.append(1))
    button = _with_class(sidebar, "delfin-session-stopall")
    notice = _with_class(sidebar, "delfin-session-notice")

    button.click()
    assert given == [], "a single click ended everything"
    assert "delfin-armed" in button._dom_classes
    assert "again" in button.description

    button.click()
    assert given == [1]
    assert button.disabled
    assert "Emergency stop given" in notice.value


def test_the_stop_runs_in_a_process_of_its_own():
    src = inspect.getsource(AS._give_emergency_stop)
    assert '"stop-all", "--yes"' in src and "start_new_session=True" in src


def test_sessions_open_on_other_nodes_are_named():
    tree = ast.parse(inspect.getsource(AS))
    tick = next(n for n in ast.walk(tree)
                if isinstance(n, ast.FunctionDef) and n.name == "_tick")
    assert "_refresh_elsewhere" in {
        c.func.id for c in ast.walk(tick)
        if isinstance(c, ast.Call) and isinstance(c.func, ast.Name)}
