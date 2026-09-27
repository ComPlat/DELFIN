"""Every session-end path also indexes the ended session (Paket 5).

Wiring test: the three public session-end doors (CLI cmd_chat's finally,
CLI cmd_run, dashboard tab_agent's distill thread) each call
index_at_session_end beside the skill-learning call, so a finished
session is searchable without anyone remembering to run a command. The
stage itself never raises (proven in test_session_end_reindex.py);
these tests pin that it is CALLED on every path.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import session_store


@pytest.fixture
def index_home(tmp_path, monkeypatch):
    from delfin.agent import session_index

    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setattr(session_store, "_SESSIONS_DIR", home / "agent_sessions")
    monkeypatch.setattr(session_store, "_transcript_archive_path",
                        lambda: home / "transcript_archive")
    monkeypatch.setattr(session_index, "_index_path",
                        lambda: home / "session_index.sqlite")
    return home


def _session_file(home, session_id, text):
    d = home / "agent_sessions"
    d.mkdir(parents=True, exist_ok=True)
    (d / f"{session_id}.json").write_text(json.dumps({
        "session_id": session_id, "title": "wiring",
        "created_at": 1760000000,
        "chat_messages": [{"role": "user", "content": text}],
    }), encoding="utf-8")


def test_the_hook_indexes_the_ended_session(index_home):
    from delfin.agent import session_end, session_index

    _session_file(index_home, "w-sess", "the europium antenna story")
    assert session_end.index_at_session_end("w-sess") is True
    assert session_index.search("europium")[0].session_id == "w-sess"


def test_all_three_session_end_doors_call_the_index_stage():
    import inspect

    from delfin.agent import cli
    from delfin.dashboard import tab_agent

    src_cli = inspect.getsource(cli)
    n_calls = src_cli.count("index_at_session_end(")
    # cmd_chat's finally AND cmd_run: two CLI doors.
    assert n_calls >= 2, (
        f"cli.py should call index_at_session_end on both doors, found {n_calls}")
    src_tab = inspect.getsource(tab_agent)
    assert "index_at_session_end(" in src_tab, (
        "the dashboard's session-end thread must index the ended session")


def test_the_stage_is_own_try_except_like_the_learning_stage():
    import inspect

    from delfin.agent import session_end
    src = inspect.getsource(session_end.index_at_session_end)
    assert "try:" in src and "except Exception" in src
