"""The real chat -r path carries the refusals, not just the engine methods.

The engine-level round trip (test_what_a_continued_session_forgets) proves
export_state/restore_state move the denials. But Welle 4 taught the other
lesson: a module can be green while its CALLER drops what it produced.
So this file goes through the path `delfin-agent chat -r <id>` actually
takes:

  save_session()          -- the real store, writing a real session file
  load_session(sid)       -- reading it back, key by key
  cli._resume_or_create() -- the function the chat CLI calls, which
                            forwards **data (the whole stored dict) into
                            engine.restore_state()

If a layer in that chain drops the "refusals" blob, these tests go red.
"""

from __future__ import annotations

import types

from pathlib import Path

from delfin.agent import cli as agent_cli
from delfin.agent import session_store as S
from delfin.agent.engine import AgentEngine


def _engine_with_denial(tmp_path: Path) -> AgentEngine:
    from delfin.agent.api_client import KitToolPermissions
    from delfin.agent.refusal_memory import Refusal, RefusalMemory

    eng = AgentEngine.__new__(AgentEngine)
    eng.mode = "solo"
    eng.route = ["solo_agent"]
    eng.current_role_index = 0
    eng.role_outputs = {}
    eng.compaction_summaries = {}
    eng.messages = []
    eng.token_usage = {"input": 0, "output": 0}
    eng.cost_usd = 0.0
    eng.session_id = "refusal-rt-1"
    for spec in AgentEngine._SESSION_FIELDS:
        setattr(eng, spec.attr, spec.reset())

    class _Client:
        pass

    perms = KitToolPermissions(workspace=tmp_path)
    perms.denied_paths = {str(tmp_path / "secret.env")}
    perms.denied_actions = {"bash:rm -rf /": "rm -rf /"}
    perms.refusal_memory = RefusalMemory()
    perms.refusal_memory.record(Refusal(
        tool="read_file", target=str(tmp_path / "secret.env"),
        reason="operator said no", time="2026-09-26T00:00:00Z"))
    eng.client = _Client()
    eng.client._permissions = perms
    return eng


def _fresh_engine(tmp_path: Path) -> AgentEngine:
    """The resumed process: client rebuilt, permissions empty."""
    eng = _engine_with_denial(tmp_path)
    perms = eng.kit_permissions
    perms.denied_paths = set()
    perms.denied_actions = {}
    perms.refusal_memory = None
    return eng


def test_the_store_round_trip_carries_the_refusals(tmp_path,
                                                   monkeypatch):
    """cli._save_session is the real saver (it forwards export_state()
    wholesale); load_session reads the file back key by key. If the store
    or its loader drops the blob, the resumed gate starts empty."""
    monkeypatch.setattr(S, "_SESSIONS_DIR", tmp_path / "sessions")
    eng = _engine_with_denial(tmp_path)
    sid = agent_cli._save_session(eng, tmp_path)

    loaded = S.load_session(sid)

    assert loaded is not None
    blob = loaded.get("refusals", {})
    assert str(tmp_path / "secret.env") in blob.get("denied_paths", [])
    assert blob.get("denied_actions") == {"bash:rm -rf /": "rm -rf /"}
    assert any("operator said no" in str(e)
               for e in blob.get("refusal_memory", {}).get("entries", []))


def test_the_chat_resume_path_restores_onto_the_gate(tmp_path,
                                                     monkeypatch):
    """cli._resume_or_create is the function `chat -r` calls. It must hand
    the loaded dict (refusals included) to restore_state, which writes
    them back onto the permission object of THIS process."""
    monkeypatch.setattr(S, "_SESSIONS_DIR", tmp_path / "sessions")
    eng = _engine_with_denial(tmp_path)
    sid = agent_cli._save_session(eng, tmp_path)

    resumed = _fresh_engine(tmp_path)
    args = types.SimpleNamespace(
        session=sid, mode="solo", verbose=False, cwd=str(tmp_path))
    monkeypatch.setattr(AgentEngine, "_sync_task_session",
                         lambda self: None)

    agent_cli._resume_or_create(resumed, args)

    perms = resumed.kit_permissions
    assert str(tmp_path / "secret.env") in perms.denied_paths
    assert perms.denied_actions == {"bash:rm -rf /": "rm -rf /"}
    assert perms.refusal_memory is not None
    assert any("operator said no" in str(e)
               for e in perms.refusal_memory.to_dict().get("entries", []))
