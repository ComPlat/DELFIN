"""A watched job that finishes wakes the agent.

"Watch the CI and tell me" meant nothing until somebody typed: a finished
watched job was put into the NEXT turn's prompt, and no next turn came on
its own. The dashboard now looks every two minutes while the agent is idle
and the input box is empty, and sends the result the way a scheduled
wake-up is sent. It peeks under its own marker, so neither the daemon's
peek nor the turn's own report of the job is taken away.
"""

from __future__ import annotations

import ast
import inspect
import pathlib

import pytest

from delfin.agent import job_monitor as jm
from delfin.dashboard import tab_agent as T

_SRC = pathlib.Path(inspect.getfile(T)).read_text(encoding="utf-8")
_SHA = "d4ecf365"


def _fn(name: str) -> ast.FunctionDef:
    for node in ast.walk(ast.parse(_SRC)):
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return node
    raise AssertionError(f"{name} not found")


def _green(url):
    return {"workflow_runs": [{
        "name": "CI", "status": "completed", "conclusion": "success",
        "head_sha": _SHA + "0" * 32, "html_url": "https://github.com/r/runs/1",
        "jobs_url": "https://api.github.com/jobs/1"}]}


@pytest.fixture
def ws(tmp_path, monkeypatch):
    monkeypatch.setattr(jm, "_AGENT_WATCH_INDEX_PATH", tmp_path / "index.json")
    jm.register_ci_watch(tmp_path, "ComPlat/DELFIN", _SHA, "main")
    return tmp_path


def _look(ws, **kw):
    path = ws / ".delfin" / "agent_watched_jobs.json"
    data = jm.load_watched(path)
    for entry in data["jobs"].values():
        entry["last_checked"] = 0
    jm.save_watched(data, path)
    return jm.check_agent_jobs(ws, fetch_fn=_green, **kw)


def test_the_wake_the_daemon_and_the_turn_each_hear_of_it_once(ws):
    assert len(_look(ws, consume=False, marker="wake_notified")) == 1
    assert _look(ws, consume=False, marker="wake_notified") == []
    assert len(_look(ws, consume=False)) == 1, "the daemon's peek is its own"
    assert len(_look(ws, consume=True)) == 1, "the turn still gets the result"
    assert jm.load_watched(ws / ".delfin" / "agent_watched_jobs.json")["jobs"] == {}


def test_the_wake_message_names_what_finished_and_links_it():
    text = T._job_wake_prompt([{
        "kind": "ci", "job_id": f"ci:ComPlat/DELFIN@{_SHA}", "state": "FAILURE",
        "ok": False, "description": f"CI for {_SHA} on main",
        "signatures": ["CI › tests (py3.11) › Run fast test suite"],
        "url": "https://github.com/ComPlat/DELFIN/actions/runs/7"}])
    # "[watch" rather than "[watch]": the marker now carries who is
    # speaking. A turn nobody typed arrives through the same input box as
    # everything else, and unsaid it was read as the user asking.
    assert text.startswith("[watch")
    assert "FAILURE" in text and "Run fast test suite" in text
    assert text.count("runs/7") == 1
    assert T._job_wake_prompt([]) == ""


def test_a_wake_message_is_not_a_request_to_push():
    from delfin.agent.api_client import KitToolPermissions, _grant_push_from
    perms = KitToolPermissions(workspace=pathlib.Path("."))
    _grant_push_from(perms, T._job_wake_prompt([{
        "kind": "ci", "job_id": f"ci:ComPlat/DELFIN@{_SHA}", "state": "SUCCESS",
        "description": f"CI for {_SHA} on main"}]), new_request=True)
    assert not perms.push_grants.get("push")


def test_only_an_idle_agent_with_an_empty_box_is_woken():
    src = ast.unparse(_fn("_job_wake_tick"))
    assert "state.get('streaming')" in src
    assert "input_textarea.value" in src
    assert "marker='wake_notified'" in src
    # It sends through the one door for turns nobody typed, which holds it
    # after an emergency stop (test_an_emergency_stop_reaches_every_machine).
    assert "_send_on_its_own(_prompt)" in src
    assert ".daemon = True" in src


def test_nothing_is_asked_of_the_scheduler_when_nothing_is_watched():
    src = ast.unparse(_fn("_job_wake_tick"))
    assert src.index("load_watched") < src.index("check_agent_jobs")


def test_one_chain_however_many_engines_are_built():
    src = ast.unparse(_fn("_wire_phase5_callbacks"))
    assert "if state.get('_job_wake_timer') is None:" in src
