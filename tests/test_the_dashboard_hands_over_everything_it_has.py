"""Three hand-written field lists in the dashboard, each one lossy.

The headless save was fixed by forwarding the engine's exported state
WHOLESALE: listing the exporter's keys by hand is how a resumed session
came back thinner than the one that was saved. The dashboard kept its own
copy of that mistake on both sides of the round trip.

- ``_auto_save_session`` named eleven exporter keys and neither the run
  clock nor the outcome-cost baseline, so a resumed dashboard session was
  handed its whole wall-clock budget again and booked the entire
  session's spend as its first turn.
- the restore next to it named seven and dropped ``schema_version`` with
  them, which made every dashboard resume look like a v1 file -- so the
  v1 migration reset the run clock to zero even on files that had
  recorded it correctly.
- the handoff brief for a LIVE session carried no workspace, and the
  brief falls back to the real task store only when it knows where that
  store is. A fresh agent was briefed that nothing was outstanding while
  work was still open.
"""

from __future__ import annotations

import ast
import inspect
import pathlib
import time
from pathlib import Path

import pytest

from delfin.agent.engine import AgentEngine
from delfin.dashboard import tab_agent as T

_SRC = pathlib.Path(inspect.getfile(T)).read_text(encoding="utf-8")


@pytest.fixture
def store(monkeypatch, tmp_path):
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    import importlib

    from delfin.agent import session_store
    return importlib.reload(session_store)


def _engine(**attrs) -> AgentEngine:
    eng = AgentEngine.__new__(AgentEngine)
    eng.mode = "solo"
    eng.messages = [{"role": "user", "content": "hallo"}]
    eng.session_id = "s1"
    for key, value in attrs.items():
        setattr(eng, key, value)
    return eng


def _call_of(func_name: str) -> ast.Call:
    """The one call to ``func_name`` in the dashboard module."""
    found = [n for n in ast.walk(ast.parse(_SRC))
             if isinstance(n, ast.Call)
             and (getattr(n.func, "id", "") == func_name
                  or getattr(n.func, "attr", "") == func_name)]
    assert len(found) == 1, f"{func_name}: {len(found)} call sites"
    return found[0]


def _exporter_keys() -> set[str]:
    return set(_engine().export_state())


# ---------------------------------------------------------------------------
# the save
# ---------------------------------------------------------------------------

def test_the_run_clock_and_the_cost_baseline_survive_the_dashboard_shape(store):
    """The exact call shape ``_auto_save_session`` uses now: the exported
    state forwarded wholesale, plus the fields only the UI knows."""
    eng = _engine(_run_started_at=time.time() - 600.0,
                  _last_outcome_cost=1.25)
    estate = dict(eng.export_state())
    estate["session_id"] = "s1"
    store.save_session(
        chat_messages=[{"role": "user", "content": "hallo"}],
        workspace="/tmp/ws",
        **estate,
    )
    data = store.load_session("s1")
    assert data["run_elapsed_s"] >= 600.0
    assert data["last_outcome_cost"] == 1.25


def test_the_saved_clock_is_still_running_after_a_restore(store):
    """A resume that reset the clock was a way to launder a time budget."""
    eng = _engine(_run_started_at=time.time() - 600.0,
                  _last_outcome_cost=1.25)
    estate = dict(eng.export_state())
    estate["session_id"] = "s1"
    store.save_session(chat_messages=[], workspace="/tmp/ws", **estate)
    data = store.load_session("s1")

    resumed = _engine()
    resumed.restore_state({**data, "mode": "solo", "session_id": "s1"})
    assert time.time() - resumed._run_started_at >= 600.0
    assert resumed._last_outcome_cost == 1.25


def test_the_save_forwards_the_exported_state_instead_of_listing_it():
    call = _call_of("save_session")
    assert any(kw.arg is None for kw in call.keywords), (
        "the dashboard save still lists the exporter's keys by hand")


def test_the_save_names_nothing_the_exporter_already_produces():
    """A name on both sides is a TypeError inside a best-effort save --
    the session would silently stop being written at all."""
    call = _call_of("save_session")
    named = {kw.arg for kw in call.keywords if kw.arg}
    assert not (named & _exporter_keys()), (
        "explicit keyword also produced by export_state: "
        f"{sorted(named & _exporter_keys())}")


# ---------------------------------------------------------------------------
# the restore
# ---------------------------------------------------------------------------

def test_the_restore_is_handed_the_saved_file_itself():
    call = _call_of("restore_state")
    assert len(call.args) == 1 and isinstance(call.args[0], ast.Dict)
    assert None in call.args[0].keys, (
        "the dashboard restore still re-lists the fields it will read")


def test_the_restore_overrides_only_what_the_ui_owns():
    call = _call_of("restore_state")
    named = {k.value for k in call.args[0].keys
             if isinstance(k, ast.Constant)}
    assert named == {"mode", "session_id"}


def test_a_session_from_a_newer_schema_is_reported_not_raised(store):
    """restore_state refuses a file it would have to read lossily. Now
    that the version reaches it, the refusal must be shown to the user."""
    import delfin.agent.session_store as _ss

    eng = _engine()
    with pytest.raises(_ss.SessionSchemaError):
        eng.restore_state({"schema_version": _ss.SESSION_SCHEMA_VERSION + 1})
    body = _SRC.split("engine.restore_state(", 1)[1][:600]
    assert "except Exception" in body and "Session not restored" in body


# ---------------------------------------------------------------------------
# the handoff
# ---------------------------------------------------------------------------

def test_a_live_handoff_reads_the_real_task_list(tmp_path):
    """The brief falls back to the task store when the payload is empty --
    but only when it knows which workspace to look in."""
    from delfin.agent.agent_tasks import get_store
    from delfin.agent.session_store import build_handoff_brief

    get_store(tmp_path).create("finish the migration", session_id="live")
    live = {
        "session_id": "live",
        "chat_messages": [{"role": "user", "content": "los"}],
        "todo_payload": [],
    }
    assert "finish the migration" not in build_handoff_brief(live)
    assert "finish the migration" in build_handoff_brief(
        {**live, "workspace": str(tmp_path)})


def test_the_live_handoff_dictionary_carries_the_workspace():
    """The saved-session path gets it from the file; the live one has to
    put it there."""
    literals = []
    for node in ast.walk(ast.parse(_SRC)):
        if not isinstance(node, ast.Dict):
            continue
        keys = {k.value for k in node.keys if isinstance(k, ast.Constant)}
        if "pending_plan_body" in keys and "active_gate" in keys:
            literals.append(keys)
    assert literals, "the live-session handoff dictionary is gone"
    for keys in literals:
        assert "workspace" in keys
