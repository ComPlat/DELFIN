"""What the agent has running in the background is on screen.

Report 20260915-132613: a background test suite ran for 40 minutes, a CI run
was watched and a wake-up was scheduled, and the dashboard showed none of
it -- while the subagent panel kept two explore agents that had finished
long before. One Background panel lists what is still out, and only that.
"""

from __future__ import annotations

import ast
import inspect
import json
import pathlib
import time
from types import SimpleNamespace

import pytest

from delfin.agent import background_view as BV
from delfin.agent import bash_jobs as BJ
from delfin.agent import job_monitor as JM
from delfin.agent import scheduler as SCH
from delfin.agent import subagents as SA

NOW = 1_789_480_000.0


@pytest.fixture
def ws(tmp_path, monkeypatch):
    monkeypatch.setattr(JM, "_AGENT_WATCH_INDEX_PATH", tmp_path / "index.json")
    monkeypatch.setattr(BJ, "live_jobs", lambda: [{"job_id": "4f064b45"}])
    monkeypatch.setattr(BJ, "_load_registry_file", lambda _ws: {"jobs": {
        "4f064b45": {"description": "Fast suite rerun", "started_at": NOW - 754},
        "6f2c39c4": {"description": "Killed first run", "started_at": NOW - 3000},
    }})
    monkeypatch.setattr(SA, "read_running", lambda **_: {
        "57a3ec39": {"type": "explore", "description": "Scratch-sink test",
                     "started_at": NOW - 80, "last_action": "read_file api_client.py"},
        "otherses": {"type": "plan", "description": "someone else's",
                     "started_at": NOW - 10},
    })
    monkeypatch.setattr(SA, "_entry_owned_by_us",
                        lambda entry: entry.get("type") != "plan")
    monkeypatch.setattr(SCH, "get_scheduler", lambda: SimpleNamespace(
        list_entries=lambda: [SimpleNamespace(
            id="fc8024e7cf", kind="once", reason="Wait for CI of 8feb1765",
            prompt="check CI", next_fire_at=NOW + 1500, disabled=False)]))
    JM.register_ci_watch(tmp_path, "ComPlat/DELFIN", "8feb1765", "main")
    return tmp_path


def test_every_kind_of_background_work_is_listed(ws):
    html = BV.render_html(BV.collect(ws, now=NOW), now=NOW)
    assert "Background · 4 running" in html
    assert "Fast suite rerun" in html and "12m34s" in html
    assert "Watch · CI" in html and "CI for 8feb1765 on main" in html
    assert "explore · Scratch-sink test" in html and "read_file api_client.py" in html
    assert "Wake-up" in html and "in 25m00s" in html


def test_what_has_ended_or_is_not_ours_is_not_listed(ws):
    html = BV.render_html(BV.collect(ws, now=NOW), now=NOW)
    assert "Killed first run" not in html
    assert "someone else" not in html


def test_nothing_out_is_an_empty_panel(tmp_path, monkeypatch):
    monkeypatch.setattr(JM, "_AGENT_WATCH_INDEX_PATH", tmp_path / "index.json")
    monkeypatch.setattr(BJ, "live_jobs", lambda: [])
    monkeypatch.setattr(SA, "read_running", lambda **_: {})
    monkeypatch.setattr(SCH, "get_scheduler",
                        lambda: SimpleNamespace(list_entries=lambda: []))
    assert BV.render_html(BV.collect(tmp_path, now=NOW), now=NOW) == ""


def test_labels_are_escaped():
    view = {"shells": [{"id": "x", "label": "<b>boom</b>", "since": NOW}],
            "watches": [], "agents": [], "wakeups": []}
    html = BV.render_html(view, now=NOW)
    assert "&lt;b&gt;boom&lt;/b&gt;" in html and "<b>boom</b>" not in html


def test_a_finished_background_agent_wakes_once_and_is_not_consumed(tmp_path, monkeypatch):
    monkeypatch.setattr(SA, "_PENDING_DIR", tmp_path)
    monkeypatch.setattr(SA, "_entry_owned_by_us", lambda entry: True)
    running = {"b2": {}}
    monkeypatch.setattr(SA, "read_running", lambda **_: running)
    for sa_id in ("a1", "b2"):
        (tmp_path / f"{sa_id}.json").write_text(json.dumps(
            {"sa_id": sa_id, "type": "explore", "description": f"run {sa_id}"}))
    seen: set = set()
    woke = BV.finished_background_agents(seen)
    assert [w["job_id"] for w in woke] == ["a1"]
    assert BV.finished_background_agents(seen) == []
    assert (tmp_path / "a1.json").exists(), "the turn drains the report, not the wake"


_TAB = pathlib.Path(inspect.getfile(
    __import__("delfin.dashboard.tab_agent", fromlist=["x"]))).read_text(encoding="utf-8")


def _fn(name: str) -> str:
    for node in ast.walk(ast.parse(_TAB)):
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return ast.unparse(node)
    raise AssertionError(f"{name} not found")


def test_the_panel_shows_what_is_out_not_what_has_finished():
    src = _fn("_refresh_subagent_panel")
    assert "background_view" in src and "render_html" in src
    assert "read_telemetry" not in src


def test_a_finished_background_agent_wakes_the_agent():
    assert "finished_background_agents" in _fn("_job_wake_tick")


def test_the_panel_refreshes_on_its_own():
    assert "_background_timer" in _fn("_wire_phase5_callbacks")
