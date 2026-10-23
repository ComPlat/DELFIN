"""A save never silently drops the model history a resume needs.

A session is saved with what its engine holds at that moment. When the
engine had been rebuilt without its history (a dashboard session after a
reload, a model switch whose hand-over was lost), that was one message, and
the save replaced the session's whole history with it. Every later resume
restored that one message while the transcript still showed everything:
seen 2026-10-14 on a dashboard session, 902 transcript entries and 1 engine
message, and the model no longer knew its task.

Shorter is legitimate after a compaction, after a new cycle (the transcript
is cleared too) and for a sliding-window trim. Anything else keeps the
stored history, appends what is new, and keeps the file as it was beside it.
Through ``save_session`` and ``load_session`` on a redirected folder.
"""
from __future__ import annotations

from unittest.mock import patch

import pytest

from delfin.agent import session_store as ss

SID = "resume-keeps-history-0001"


@pytest.fixture
def store(tmp_path):
    d = tmp_path / "agent_sessions"
    d.mkdir()
    with patch("delfin.agent.session_store._SESSIONS_DIR", d):
        yield d


def _history(n_pairs):
    out = []
    for i in range(n_pairs):
        out.append({"role": "user", "content": f"question {i}"})
        out.append({"role": "assistant", "content": f"answer {i}"})
    return out


def _chat(n):
    return [{"role": "user" if i % 2 == 0 else "assistant",
             "content": f"line {i}"} for i in range(n)]


def _save(engine_messages, chat_messages, **kw):
    ss.save_session(SID, mode="solo", chat_messages=chat_messages,
                    engine_messages=engine_messages, **kw)
    return ss.load_session(SID)["engine_messages"]


def test_a_fresh_engine_does_not_replace_the_stored_history(store):
    _save(_history(100), _chat(400))
    new_turn = [{"role": "user", "content": "go on with the task"},
                {"role": "assistant", "content": "continuing"}]
    loaded = _save(new_turn, _chat(402))
    assert len(loaded) == 202
    assert loaded[:200] == _history(100)
    assert loaded[-2:] == new_turn
    assert list(store.glob(f"{SID}.json.history-kept-*")), "the file as it was"


def test_what_the_session_already_had_is_not_added_twice(store):
    _save(_history(100), _chat(400))
    loaded = _save(_history(100)[-2:] + [{"role": "user", "content": "next"}],
                   _chat(401))
    assert len(loaded) == 201 and loaded[-1]["content"] == "next"


def test_a_compacted_history_is_saved_as_it_is(store):
    _save(_history(100), _chat(400))
    compacted = [{"role": "user", "content": "[Conversation summary — older "
                  "messages compacted]\nwhat happened"},
                 {"role": "assistant", "content": "Understood."}]
    assert _save(compacted, _chat(402)) == compacted


def test_a_changed_compaction_record_explains_a_shorter_history(store):
    _save(_history(100), _chat(400))
    short = [{"role": "user", "content": "after compaction"}]
    assert _save(short, _chat(401),
                 last_compaction_info={"at": 1, "removed": 199}) == short


def test_a_new_cycle_starts_empty(store):
    _save(_history(100), _chat(400))
    fresh = [{"role": "user", "content": "a new task"}]
    assert _save(fresh, _chat(1)) == fresh


def test_a_sliding_window_trim_is_saved_as_it_is(store):
    _save(_history(100), _chat(400))
    trimmed = _history(100)[20:]
    assert _save(trimmed, _chat(402)) == trimmed


def test_a_first_save_is_written_as_it_is(store):
    assert _save([{"role": "user", "content": "hello"}], _chat(1)) == [
        {"role": "user", "content": "hello"}]
