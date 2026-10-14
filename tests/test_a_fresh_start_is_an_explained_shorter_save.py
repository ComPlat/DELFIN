"""An answered fresh start is a shorter save the history guard accepts.

Input: a session whose history exceeds the model's window (a switch to a
smaller model), one turn through ``prepare_fresh_turn`` and the real
engine, answered, then ``save_session`` and ``load_session``.

Output: the saved engine history is the fresh one (the fused prompt and its
answer) and the transcript archive holds the history the fresh start cut.

Semantics: ``save_session`` keeps the stored history when a save would make
it shorter without an explanation. A fresh start is an explanation: it
archived what it cut, and its first message carries the fresh-context
marker. Restoring the archived history into the engine file would put the
session over the window again at every resume, and the next turn would cut
it again.
"""
from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pytest

from delfin.agent import session_store as ss
from delfin.agent.engine import AgentEngine
from delfin.agent.repl import _FRESH_CONTEXT_MARKER, prepare_fresh_turn

from tests.test_t2_engine_fresh_start import _pack_tree

SID = "fresh-start-save-0001"


@pytest.fixture
def home(tmp_path, monkeypatch):
    root = tmp_path / "home"
    root.mkdir()
    monkeypatch.setattr(Path, "home", lambda: root)
    sessions = root / ".delfin" / "agent_sessions"
    sessions.mkdir(parents=True)
    with patch("delfin.agent.session_store._SESSIONS_DIR", sessions):
        yield root


class _Answer:
    def stream_message(self, *, messages, system="", max_tokens=0, **kw):
        self.sent = list(messages or [])
        yield SimpleNamespace(type="text_delta", text="continuing the task")
        yield SimpleNamespace(type="message_delta", input_tokens=40,
                              output_tokens=3, cost_usd=0.0,
                              stop_reason="end_turn", text="")


def _engine(tmp_path):
    pack = _pack_tree(tmp_path / "pack_root")
    client = _Answer()
    with patch("delfin.agent.engine.create_client", return_value=client):
        engine = AgentEngine(repo_dir=pack, backend="cli", mode="quick",
                             pack_dir=pack)
    engine.client = client
    engine.session_id = SID
    return engine


def _history(n_pairs):
    out = []
    for i in range(n_pairs):
        out.append({"role": "user", "content": f"question {i} " + "x" * 200})
        out.append({"role": "assistant", "content": f"answer {i} " + "y" * 200})
    return out


def _chat(n):
    return [{"role": "user" if i % 2 == 0 else "assistant",
             "content": f"line {i}"} for i in range(n)]


def test_an_answered_fresh_start_is_saved_as_the_fresh_history(home, tmp_path):
    engine = _engine(tmp_path)
    before = _history(60)
    engine.messages = [dict(m) for m in before]
    ss.save_session(SID, mode="solo", chat_messages=_chat(120),
                    engine_messages=list(engine.messages))

    msg = prepare_fresh_turn(engine, "go on", context_budget=500)
    assert msg.startswith(_FRESH_CONTEXT_MARKER), "the turn started fresh"
    engine.stream_response(user_message=msg)
    assert [m["role"] for m in engine.messages] == ["user", "assistant"]

    ss.save_session(SID, mode="solo", chat_messages=_chat(122),
                    engine_messages=list(engine.messages))
    loaded = ss.load_session(SID)["engine_messages"]

    assert len(loaded) == 2, len(loaded)
    assert str(loaded[0]["content"]).startswith(_FRESH_CONTEXT_MARKER)
    assert loaded[1]["content"] == "continuing the task"
    assert not list((home / ".delfin" / "agent_sessions").glob(
        f"{SID}.json.history-kept-*")), "nothing was kept aside"

    records = ss.load_transcript_archive(SID)
    fresh = [r for r in records
             if (r.get("info") or {}).get("kind") == "fresh_restart"]
    assert fresh, records
    assert [m["content"] for m in fresh[-1]["messages"]] == [
        m["content"] for m in before]


def test_a_marker_without_an_archive_record_is_not_an_explanation(home):
    """The marker alone could be any text; the archive record is what says
    the cut history is recoverable. Without it the stored history stays."""
    ss.save_session(SID, mode="solo", chat_messages=_chat(120),
                    engine_messages=_history(60))
    forged = [{"role": "user", "content": _FRESH_CONTEXT_MARKER + "\n\ngo on"},
              {"role": "assistant", "content": "ok"}]
    ss.save_session(SID, mode="solo", chat_messages=_chat(122),
                    engine_messages=forged)
    assert len(ss.load_session(SID)["engine_messages"]) == 122


def test_the_marker_the_guard_reads_is_the_one_the_fresh_start_writes():
    assert ss._FRESH_PREFIX == _FRESH_CONTEXT_MARKER
