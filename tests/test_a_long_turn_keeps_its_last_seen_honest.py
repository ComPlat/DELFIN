"""The stored session's "last seen" moves while a turn is still running.

Input:  an engine driven through ``_record_tool_trace`` -- the per-tool-call
        path the real turn takes -- with a redirected session store.
Output: whether ``heartbeat_at`` of the STORED session moves forward, how
        often it is rewritten, and what happens with no session id.
Semantics: ``delfin-agent sessions`` reads ``heartbeat_at`` as "last seen"
        (``cli_resume._load_row``, and the roster at :218). Only a save
        moved it, so a session working through one long turn was listed as
        last seen when the turn began -- while ``session_presence.touch``
        was refreshing the presence record beside it on every tool call.
        Two records of the same fact and only one of them refreshed;
        ``session_store.touch_heartbeat`` existed for exactly this and had
        no caller in the product.

Nothing here touches a real store: the sessions directory is redirected to
tmp_path and the session id is a stand-in.

Run: ~/micromamba/envs/delfin-ci/bin/python -m pytest \
     tests/test_a_long_turn_keeps_its_last_seen_honest.py -q
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import engine as engine_mod
from delfin.agent import session_store

SESSION = "s-long-turn-0001"


@pytest.fixture
def stored(tmp_path, monkeypatch):
    """A stored session whose heartbeat is old, in a throwaway store."""
    monkeypatch.setattr(session_store, "_SESSIONS_DIR", tmp_path)
    path = tmp_path / f"{SESSION}.json"
    path.write_text(json.dumps({"session_id": SESSION, "heartbeat_at": 1.0}))
    return path


def _engine(session_id: str = SESSION):
    """The engine as the existing presence tests build it: the real class,
    with the attributes the per-tool-call path reads."""
    eng = engine_mod.AgentEngine.__new__(engine_mod.AgentEngine)
    eng._trace_pending = []
    # kit_permissions is a read-only property over the CLIENT's
    # permissions, so the presence key is supplied where it looks.
    eng.client = type("_C", (), {"_permissions": type(
        "_P", (), {"presence_key": "pk-long-turn"})()})()
    eng.trace_session = lambda: session_id
    eng.session_id = session_id
    eng._last_session_heartbeat = 0.0
    return eng


def test_a_tool_call_moves_the_stored_heartbeat(stored):
    eng = _engine()
    eng._record_tool_trace("bash", "out")
    after = json.loads(stored.read_text())
    assert after["heartbeat_at"] > 1.0, after
    assert after["session_id"] == SESSION, "nothing else is rewritten"


def test_the_rewrite_is_rate_limited(stored):
    """The other direction, and the reason this is safe per tool call:
    touch_heartbeat rewrites the whole record, so a turn with a hundred
    tool calls must not mean a hundred rewrites."""
    eng = _engine()
    eng._record_tool_trace("bash", "out")
    first = stored.stat().st_mtime_ns
    for _ in range(20):
        eng._record_tool_trace("bash", "out")
    assert stored.stat().st_mtime_ns == first, (
        "the limiter is what makes this safe to call on every tool call")


def test_a_session_without_an_id_writes_nothing(tmp_path, monkeypatch):
    """A backend whose id arrives later from the stream starts with "".
    No file may be created for it: a nameless record on the roster is
    worse than a stale one."""
    monkeypatch.setattr(session_store, "_SESSIONS_DIR", tmp_path)
    eng = _engine(session_id="")
    eng._record_tool_trace("bash", "out")
    # The session RECORDS only: the redirected directory also holds the
    # stand-in home the presence side creates, and asserting on the whole
    # listing measured that instead.
    assert list(tmp_path.glob("*.json")) == []


def test_an_unstored_session_is_not_a_crash(tmp_path, monkeypatch):
    """The id exists but nothing is stored yet (before the first save).
    touch_heartbeat returns False and the turn continues."""
    monkeypatch.setattr(session_store, "_SESSIONS_DIR", tmp_path)
    eng = _engine()
    eng._record_tool_trace("bash", "out")
    assert list(tmp_path.glob("*.json")) == []


def test_the_product_reaches_the_writer(stored):
    """The reachability check, which is what was wrong: the capability
    existed, was tested, and nothing in delfin/ called it."""
    from pathlib import Path

    root = Path(engine_mod.__file__).resolve().parent
    callers = [p.name for p in root.rglob("*.py")
               if p.name != "session_store.py"
               and "touch_heartbeat" in p.read_text()]
    assert "engine.py" in callers, callers
    baseline = (root.parent.parent / "tests" / "fixtures"
                / "unmounted_capabilities.txt")
    if baseline.is_file():
        assert "touch_heartbeat" not in baseline.read_text()
