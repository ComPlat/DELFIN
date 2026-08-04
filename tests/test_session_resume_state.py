"""What a resumed session gets back — for every mode, not just one.

The engine exports its full state and the store keeps what it is handed.
Two fields were exported and then dropped on the way out, so a resumed
session silently started without them: the directory it had been pinned
to, and the context estimator's floor. The second one is visible — the
window bar reads far emptier than the conversation is, until the next
turn re-establishes the count from the provider — and the first one is
not visible at all, which is worse.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import pytest


@pytest.fixture
def store(monkeypatch, tmp_path):
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    import importlib

    from delfin.agent import session_store
    return importlib.reload(session_store)


def test_the_conversation_itself_round_trips(store):
    store.save_session("s1", mode="solo",
                       engine_messages=[{"role": "user", "content": "hallo"},
                                        {"role": "assistant", "content": "ja"}],
                       chat_messages=[{"role": "user", "content": "hallo"}])
    data = store.load_session("s1")
    assert len(data["engine_messages"]) == 2
    assert data["engine_messages"][0]["content"] == "hallo"


def test_the_estimator_floor_survives(store):
    """Without it the context bar under-reports a resumed session."""
    store.save_session("s1", mode="office", last_input_tokens=48_000)
    assert store.load_session("s1")["last_input_tokens"] == 48_000


def test_the_project_pin_survives(store):
    """restore_state reads it; it was never written."""
    store.save_session("s1", mode="solo", project_dir="/pfad/zum/projekt")
    assert store.load_session("s1")["project_dir"] == "/pfad/zum/projekt"


def test_everything_the_engine_exports_can_be_saved():
    """The exporter and the store are two halves of one round trip. When
    a field is added to one and not the other it is lost in silence."""
    import inspect

    from delfin.agent import session_store
    from delfin.agent.engine import AgentEngine

    exported = set()
    source = inspect.getsource(AgentEngine.export_state)
    for line in source.splitlines():
        line = line.strip()
        if line.startswith('"') and '":' in line:
            exported.add(line.split('"')[1])

    accepted = set(inspect.signature(session_store.save_session).parameters)
    # route/role_outputs and friends are named identically on both sides;
    # session_id is the argument the record is keyed by.
    for field in exported:
        assert field in accepted or field in {"session_id"}, (
            f"engine exports {field!r} and save_session cannot take it")


def test_the_restore_reads_what_the_store_keeps():
    """The dashboard hands restore_state a dict it builds by hand; a field
    kept by the store but not passed there is still lost."""
    import inspect

    from delfin.dashboard import tab_agent

    source = inspect.getsource(tab_agent)
    block = source.split("engine.restore_state({", 1)[1][:900]
    for field in ("engine_messages", "token_usage", "cost_usd",
                  "project_dir", "last_input_tokens"):
        assert field in block, field


def test_resume_state_is_not_mode_specific():
    """Nothing in the save or restore path branches on the mode."""
    import inspect

    from delfin.agent import session_store
    from delfin.agent.engine import AgentEngine

    for func in (session_store.save_session, AgentEngine.restore_state,
                 AgentEngine.export_state):
        source = inspect.getsource(func)
        assert '"office"' not in source and '"solo"' not in source
