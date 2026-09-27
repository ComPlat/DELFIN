"""The session_search tool answers over the real tool path.

Phase 3 of the searchable-sessions package: the tool is registered in
api_client (definition, dispatch, executor) exactly like s14's
skill_propose_patch pattern, and these tests drive it through
``_DocToolExecutor.execute`` -- the way the engine calls it.

Hits are marked as untrusted data (a snippet is text from a past
session's transcript), the tool only ever reads, and no test touches the
real ``~/.delfin``.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import session_store


@pytest.fixture
def tool_env(tmp_path, monkeypatch):
    """Redirected index + sources, and a bare executor."""
    from delfin.agent import api_client as A
    from delfin.agent import session_index

    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setattr(session_store, "_SESSIONS_DIR", home / "agent_sessions")
    monkeypatch.setattr(session_store, "_transcript_archive_path",
                        lambda: home / "transcript_archive")
    monkeypatch.setattr(session_index, "_index_path",
                        lambda: home / "session_index.sqlite")
    d = home / "agent_sessions"
    d.mkdir(parents=True)
    (d / "sess-tool.json").write_text(json.dumps({
        "session_id": "sess-tool", "title": "Fe complex spin state",
        "created_at": 1760000000,
        "chat_messages": [
            {"role": "user", "content": "how did we solve the Fe complex "
                                        "spin state problem?"},
            {"role": "assistant", "content": "broken symmetry PBE0 worked"}],
    }), encoding="utf-8")
    session_index.index_session("sess-tool")
    ex = A._DocToolExecutor.__new__(A._DocToolExecutor)
    # A permissions object with a workspace, the way a real client has
    # one. Without it the sandbox pre-gate rightly refuses the tool --
    # `session_search` is not on the no-sandbox allow-list and does not
    # ask to be (that list is a documented security decision).
    perms = A.KitToolPermissions(workspace=tmp_path, mode="default")
    return ex, perms


class TestSessionSearchTool:
    def test_the_tool_finds_an_indexed_session(self, tool_env):
        ex, perms = tool_env
        out = ex.execute("session_search",
                         {"query": "Fe complex"}, perms)
        assert '"sess-tool"' in out
        assert "Fe complex spin state" in out

    def test_hits_are_marked_as_untrusted_data(self, tool_env):
        ex, perms = tool_env
        out = ex.execute("session_search",
                         {"query": "Fe complex"}, perms)
        assert "UNTRUSTED" in out, (
            "snippets are text from past sessions: data, not instructions")

    def test_limit_is_passed_and_capped(self, tool_env):
        ex, perms = tool_env
        out = ex.execute("session_search",
                         {"query": "Fe complex", "limit": 99}, perms)
        payload = json.loads(out.strip().splitlines()[1]
                             if out.startswith("[UNTRUSTED") else out)
        assert isinstance(payload, list)
        assert len(payload) <= 20  # MAX_LIMIT is the hard ceiling

    def test_empty_query_is_an_error_not_a_crash(self, tool_env):
        ex, perms = tool_env
        out = ex.execute("session_search", {"query": ""}, perms)
        assert "error" in out

    def test_a_secret_is_not_retrievable(self, tool_env, monkeypatch):
        import time as _time

        from delfin.agent import session_index
        d = session_store._SESSIONS_DIR
        token = "AKIA" + "WXYZ" * 4
        (d / "sess-sec.json").write_text(json.dumps({
            "session_id": "sess-sec", "title": "credentials",
            "created_at": _time.time(),
            "chat_messages": [
                {"role": "user", "content": f"use {token} for zinc"}],
        }), encoding="utf-8")
        session_index.index_session("sess-sec")
        ex, perms = tool_env
        out = ex.execute("session_search", {"query": "zinc"}, perms)
        assert token not in out

    def test_the_tool_definition_is_registered(self):
        from delfin.agent import api_client as A
        names = [t.get("function", {}).get("name", "")
                 for t in A._DOC_TOOLS_OPENAI]
        assert "session_search" in names
