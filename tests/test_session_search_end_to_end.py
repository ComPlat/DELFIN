"""End-to-end: a finished session is findable through the real tool path.

Paket 5's integration test over the public call route, not module
internals: a session is created the way DELFIN creates one
(session_store.save_session), ended the way the session-end step will
call it (session_end_index.reindex_finished_session), and found through
the tool the model would call (_DocToolExecutor.execute("session_search",
...)). A secret that was in the session must not be retrievable.
"""

from __future__ import annotations

import pytest

from delfin.agent import session_store

_UNIQUE_TERM = "zirkonium-oxychloride-benchmark-2026"


@pytest.fixture
def real_home(tmp_path, monkeypatch):
    from delfin.agent import session_index

    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setattr(session_store, "_SESSIONS_DIR", home / "agent_sessions")
    monkeypatch.setattr(session_store, "_transcript_archive_path",
                        lambda: home / "transcript_archive")
    monkeypatch.setattr(session_index, "_index_path",
                        lambda: home / "session_index.sqlite")
    return home


def test_end_to_end_over_the_public_route(real_home, tmp_path):
    from delfin.agent import api_client as A
    from delfin.agent import session_end_index, session_index

    # 1. A session exists, created the way DELFIN creates sessions.
    session_store.save_session(
        "e2e-sess", mode="quick",
        title=f"the {_UNIQUE_TERM} episode",
        chat_messages=[
            {"role": "user",
             "content": f"how did we tune the {_UNIQUE_TERM} run?"},
            {"role": "assistant", "content": "lower lambda, more salt"},
        ],
        token_usage={"input": 10, "output": 5}, cost_usd=0.01)

    # Not indexed yet -> not findable.
    assert not session_index.search(_UNIQUE_TERM)

    # 2. Session ends: the session-end step runs.
    assert session_end_index.reindex_finished_session("e2e-sess") is True

    # 3. The model's tool finds it, over the real executor path.
    ex = A._DocToolExecutor.__new__(A._DocToolExecutor)
    perms = A.KitToolPermissions(workspace=tmp_path, mode="default")
    out = ex.execute("session_search", {"query": _UNIQUE_TERM}, perms)
    assert "e2e-sess" in out
    assert "zirkonium" in out
    assert "UNTRUSTED" in out  # hits are data, not instructions


def test_end_to_end_a_secret_stays_unfindable(real_home, tmp_path):
    from delfin.agent import api_client as A
    from delfin.agent import session_end_index

    token = "AKIA" + "E2E2" * 4
    session_store.save_session(
        "e2e-secret", mode="quick", title="credentials session",
        chat_messages=[
            {"role": "user", "content": f"rhodium work, key {token}"}],
        token_usage={"input": 1, "output": 1}, cost_usd=0.0)
    session_end_index.reindex_finished_session("e2e-secret")

    ex = A._DocToolExecutor.__new__(A._DocToolExecutor)
    perms = A.KitToolPermissions(workspace=tmp_path, mode="default")
    # The session is findable by its innocent content...
    out = ex.execute("session_search", {"query": "rhodium"}, perms)
    assert "e2e-secret" in out
    # ...but the credential is not retrievable, by content or by search.
    assert token not in out
    out2 = ex.execute("session_search", {"query": token}, perms)
    assert token not in out2
