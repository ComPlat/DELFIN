"""A session's name is never a slash command, and renaming from the phone works.

Seen live (2026-11-12): a session whose first line was the dashboard command
``/sync`` got the title ``/sync``, and ``name Test`` in its SimpleX group did
not rename it.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent import session_store
from delfin.agent.repl_commands import BUILTINS, looks_like_command


@pytest.fixture
def sessions(tmp_path, monkeypatch):
    monkeypatch.setattr(session_store, "_SESSIONS_DIR", Path(tmp_path))
    return tmp_path


def _saved(sessions, sid: str) -> dict:
    return json.loads((Path(sessions) / f"{sid}.json").read_text())


def _save(sessions, sid: str, first: str) -> None:
    session_store.save_session(
        session_id=sid, chat_messages=[{"role": "user", "content": first},
                                       {"role": "assistant", "content": "ok"}],
        workspace=str(sessions))


def test_a_builtin_command_is_command_shaped():
    assert looks_like_command("/sync")
    assert "/sync" in BUILTINS


@pytest.mark.parametrize("first, kept", [
    ("/sync", False),                       # the dashboard command as typed
    ("/sync name Butan", False),            # the command with its argument
    ("/status", False),
    ("/mode plan", False),
    ("/home/user/notes.md", True),          # a path someone pasted: a name
    ("Fix the parser", True),               # a real task
])
def test_the_first_line_names_a_session_only_when_it_is_not_a_command(
        sessions, first, kept):
    # An untitled session gets its stand-in name from its first user message;
    # a builtin command as the first line must never become the title -- the
    # next real prompt names it (and until one comes the session is untitled).
    _save(sessions, "cmd-test", first)
    assert _saved(sessions, "cmd-test").get("title") == (first if kept else "")
