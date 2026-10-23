"""A compaction tried before a fresh start, through save, load and Stop.

Input: a session whose history exceeds the fresh-start budget, one turn
through ``prepare_fresh_turn`` and the real engine, with a scripted client
that answers the summariser's request and then either answers the turn or
stops before its first token. The archive and session store are real, under
a temporary home.

Output:

(a) the compaction fits, the turn answers, the session is saved and loaded:
    the saved history begins with the compaction summary, the save guard
    accepts the shorter history as explained (no ``history-kept-*`` file,
    no restored archive rows), and the archive holds the compacted messages.
(b) the compaction fits, Stop before the first answer: the engine holds the
    compacted history -- the turn's own row is taken back, the compaction is
    not undone and the engine is not emptied; summary-side archive plus the
    kept recent rows together are the whole pre-turn history.
(c) the compaction fails (raises) or does not fit, the turn starts fresh,
    Stop before the first answer: the fresh start's take-back restores the
    history it cut -- the original one, or the compacted one when the
    compaction ran and did not fit.
"""
from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pytest

from delfin.agent import repl as R
from delfin.agent import session_store as ss
from delfin.agent.engine import AgentEngine

from tests.test_t2_engine_fresh_start import _pack_tree

SID = "compaction-before-fresh-0001"
WINDOW = 4_000


@pytest.fixture
def home(tmp_path, monkeypatch):
    root = tmp_path / "home"
    root.mkdir()
    monkeypatch.setattr(Path, "home", lambda: root)
    sessions = root / ".delfin" / "agent_sessions"
    sessions.mkdir(parents=True)
    with patch("delfin.agent.session_store._SESSIONS_DIR", sessions):
        yield root


class _Client:
    """Answers the summariser; the turn itself answers or stops."""

    def __init__(self, how="answer"):
        self.how = how
        self.engine = None
        self.summary_calls = 0
        self.sent = None

    def stream_message(self, *, messages, system="", max_tokens=0, **kw):
        if "summarisation assistant" in (system or ""):
            self.summary_calls += 1
            yield SimpleNamespace(type="text_delta",
                                  text="## Goal\nRefine the cobalt complex.")
            yield SimpleNamespace(type="message_delta", input_tokens=0,
                                  output_tokens=0, cost_usd=0.0,
                                  stop_reason="end_turn", text="")
            return
        self.sent = list(messages or [])
        if self.how == "stop":
            self.engine.request_stop()
            return
        yield SimpleNamespace(type="text_delta", text="continuing the task")
        yield SimpleNamespace(type="message_delta", input_tokens=40,
                              output_tokens=3, cost_usd=0.0,
                              stop_reason="end_turn", text="")


def _engine(tmp_path, how="answer"):
    pack = _pack_tree(tmp_path / "pack_root")
    client = _Client(how)
    with patch("delfin.agent.engine.create_client", return_value=client):
        engine = AgentEngine(repo_dir=pack, backend="cli", mode="quick",
                             pack_dir=pack)
    engine.client = client
    client.engine = engine
    # An API-style backend: the summariser is a stateless side call there.
    engine.backend = "api"
    engine.context_window_tokens = WINDOW
    engine.session_id = SID
    return engine, client


def _history(n_pairs, size=400):
    out = []
    for i in range(n_pairs):
        out.append({"role": "user", "content": f"step {i} " + "x" * size})
        out.append({"role": "assistant", "content": f"done {i} " + "y" * size})
    return out


def _chat(n):
    return [{"role": "user" if i % 2 == 0 else "assistant",
             "content": f"line {i}"} for i in range(n)]


def _contents(messages):
    return [m["content"] for m in messages]


def _records(kind=None):
    records = ss.load_transcript_archive(SID)
    if kind is None:
        return [r for r in records
                if (r.get("info") or {}).get("kind") != "fresh_restart"]
    return [r for r in records if (r.get("info") or {}).get("kind") == kind]


def _forced_compaction_fails(monkeypatch):
    original = AgentEngine._compact_history

    def _compact(self, *, force=False):
        if force:
            raise RuntimeError("summariser unavailable")
        return original(self, force=force)
    monkeypatch.setattr(AgentEngine, "_compact_history", _compact)


# (a) ----------------------------------------------------------------------

@pytest.mark.parametrize("pass_info", [True, False],
                         ids=["with-compaction-record", "summary-alone"])
def test_an_answered_compaction_is_saved_as_the_compacted_history(
        home, tmp_path, pass_info):
    """With the compaction record (as the dashboard saves) and without it
    (the summary at the head of the history is the explanation on its own)."""
    engine, client = _engine(tmp_path)
    before = _history(30)
    engine.messages = [dict(m) for m in before]
    ss.save_session(SID, mode="solo", chat_messages=_chat(60),
                    engine_messages=list(engine.messages))

    msg = R.prepare_fresh_turn(engine, "go on")
    assert msg == "go on", "the compaction fit; no fresh start"
    assert client.summary_calls == 1
    engine.stream_response(user_message=msg)
    assert engine.messages[-1]["content"] == "continuing the task"

    info = engine.last_compaction_info if pass_info else None
    ss.save_session(SID, mode="solo", chat_messages=_chat(62),
                    engine_messages=list(engine.messages),
                    last_compaction_info=info)
    loaded = ss.load_session(SID)["engine_messages"]

    assert len(loaded) < len(before)
    assert str(loaded[0]["content"]).startswith(ss._SUMMARY_MARKER), (
        str(loaded[0]["content"])[:120])
    assert "cobalt" in loaded[0]["content"]
    assert loaded[-1]["content"] == "continuing the task"
    assert _contents(loaded) == _contents(engine.messages)
    assert not list((home / ".delfin" / "agent_sessions").glob(
        f"{SID}.json.history-kept-*")), "the guard kept the old history aside"
    archived = set(_contents(before[:-engine._KEEP_RECENT]))
    assert not archived & set(_contents(loaded)), (
        "an archived row came back into the saved history")

    compactions = _records()
    assert compactions, "the compaction archived what it summarised"
    assert _contents(compactions[-1]["messages"]) == _contents(
        before[:-engine._KEEP_RECENT])
    assert not _records("fresh_restart"), "no fresh start ran"


# (b) ----------------------------------------------------------------------

def test_a_stop_after_a_compaction_keeps_the_compacted_history(home, tmp_path):
    engine, client = _engine(tmp_path, how="stop")
    before = _history(30)
    engine.messages = [dict(m) for m in before]

    msg = R.prepare_fresh_turn(engine, "go on")
    assert msg == "go on", "the compaction fit; no fresh start"
    compacted = [dict(m) for m in engine.messages]
    assert len(compacted) < len(before)

    engine.stream_response(user_message=msg)
    assert client.sent is not None, "the turn reached the transport"

    assert engine.messages, "the stop emptied the engine"
    assert _contents(engine.messages) == _contents(compacted), (
        "the stop took back more (or less) than the turn's own row")
    assert str(engine.messages[0]["content"]).startswith(ss._SUMMARY_MARKER)
    assert not getattr(engine, "start_fresh", False)
    # Nothing lost: the summarised rows are in the archive, the recent
    # ones are still in the engine.
    recent = before[-engine._KEEP_RECENT:]
    assert _contents(engine.messages[-len(recent):]) == _contents(recent)
    archived = _records()[-1]["messages"]
    assert _contents(archived) + _contents(recent) == _contents(before)
    assert not _records("fresh_restart")


# (c) ----------------------------------------------------------------------

def test_a_stop_after_a_failed_compaction_restores_the_fresh_cut(
        home, tmp_path, monkeypatch):
    _forced_compaction_fails(monkeypatch)
    engine, client = _engine(tmp_path, how="stop")
    before = _history(30)
    engine.messages = [dict(m) for m in before]

    msg = R.prepare_fresh_turn(engine, "go on")
    assert msg.startswith(R._FRESH_CONTEXT_MARKER), "the turn started fresh"
    engine.stream_response(user_message=msg)
    assert client.sent is not None and len(client.sent) == 1

    assert _contents(engine.messages) == _contents(before)
    assert engine.start_fresh is False
    assert _records("fresh_restart"), "the fresh start archived its cut"


def test_a_stop_after_a_compaction_that_did_not_fit_restores_the_compacted(
        home, tmp_path):
    """The compaction ran and is kept; the fresh start then cut the
    compacted history, and the take-back gives that one back."""
    engine, client = _engine(tmp_path, how="stop")
    # The four kept recent messages alone exceed the window.
    before = _history(10) + _history(2, size=10_000)
    engine.messages = [dict(m) for m in before]

    msg = R.prepare_fresh_turn(engine, "go on")
    assert msg.startswith(R._FRESH_CONTEXT_MARKER), "the turn started fresh"
    assert client.summary_calls == 1
    engine.stream_response(user_message=msg)

    assert engine.messages, "the stop emptied the engine"
    assert str(engine.messages[0]["content"]).startswith(ss._SUMMARY_MARKER)
    recent = before[-engine._KEEP_RECENT:]
    assert _contents(engine.messages[-len(recent):]) == _contents(recent)
    cut = _records("fresh_restart")[-1]["messages"]
    assert _contents(engine.messages) == _contents(cut)
    assert _contents(_records()[-1]["messages"]) + _contents(recent) == \
        _contents(before)
    assert engine.start_fresh is False
