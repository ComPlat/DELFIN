"""The agent `remember` tool — proactively save durable facts to project memory
(mirrors Claude Code's memory behaviour: the agent itself, not just the user,
can persist a typed memory mid-conversation)."""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import patch

import pytest

from delfin.agent import api_client as A
from delfin.agent.api_client import KitToolPermissions


def _perms(ws: Path) -> KitToolPermissions:
    return KitToolPermissions(workspace=ws, mode="bypassPermissions")


def test_remember_saves_typed_memory(tmp_path):
    home = tmp_path / "home"
    home.mkdir()
    ws = tmp_path / "ws"
    ws.mkdir()
    with patch.object(Path, "home", lambda: home):
        out = json.loads(A._doc_executor._execute_remember(
            {"text": "The user prefers German.", "type": "user"}, _perms(ws)))
        assert out["status"] == "ok" and out["type"] == "user"
        from delfin.agent.memory_store import list_typed_memories
        mems = list_typed_memories(ws)
    assert mems and len(mems) == 1
    # MEMORY.md index was created under the project's memory dir
    mem_md = home / ".delfin" / "projects"
    assert mem_md.exists()


def test_remember_requires_text(tmp_path):
    out = json.loads(A._doc_executor._execute_remember({"type": "user"},
                                                       _perms(tmp_path)))
    assert "error" in out


def test_remember_type_defaults_and_parses_prefix(tmp_path):
    home = tmp_path / "h"
    home.mkdir()
    ws = tmp_path / "w"
    ws.mkdir()
    with patch.object(Path, "home", lambda: home):
        # no explicit type → parsed/defaulted (a 'feedback:' prefix is honoured)
        out = json.loads(A._doc_executor._execute_remember(
            {"text": "feedback: always live-test changes. Why: mocks miss bugs."},
            _perms(ws)))
    assert out["status"] == "ok" and out["type"] == "feedback"


def test_remember_wired_end_to_end():
    src = (Path(__file__).resolve().parent.parent / "delfin" / "agent"
           / "api_client.py").read_text(encoding="utf-8")
    assert '"name": "remember"' in src              # advertised tool
    assert 'name == "remember"' in src              # dispatch route
    assert '"remember",' in src                      # dashboard allow-list
    addendum = (Path(__file__).resolve().parent.parent / "delfin" / "agent"
                / "pack" / "shared" / "memory_addendum.md")
    assert addendum.is_file()                        # prompt guidance exists
    pl = (Path(__file__).resolve().parent.parent / "delfin" / "agent"
          / "prompt_loader.py").read_text(encoding="utf-8")
    assert "memory_addendum" in pl                   # injected into the prompt
