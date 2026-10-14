"""A compacted history is an explained shorter save, working-state block or not.

Input: a session over the fresh-start budget, one terminal turn
(``repl.run_turn``) whose forced compaction fits and writes a working-state
block ahead of the summary, then a save.

Output: the save keeps the compacted history (no ``history-kept-*`` file, no
archived row restored), whether or not the save passes the compaction
record. The headless saver (``cli._save_session``: terminal, ``-p``/run,
scheduler) writes the record, and ``restore_state`` reads it back, so the
record a resumed session compares against is the one it saved.

Semantics: ``save_session`` keeps the stored history when a save shortens
it without an explanation. A compaction is one. It was recognised by the
first message starting with ``"[Conversation summary"``, but the compaction
composes the working-state block AHEAD of the summary, so a save without the
compaction record kept the long history and appended the compacted one --
the next resume restored rows the compaction had archived. The summary is
now accepted at the head of the first message or directly after a
working-state block, by the exact markers the compaction writes.
"""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from delfin.agent import repl as R
from delfin.agent import session_store as ss
from delfin.agent import working_state as ws
from delfin.agent.engine import AgentEngine

from tests.test_a_compaction_before_a_fresh_start_saves_and_stops import (
    SID, _chat, _contents, _engine, _history)

_BLOCK = ws._HEADER + "\nOpen tasks (task tool):\n  refine the cobalt complex\n"


@pytest.fixture
def home(tmp_path, monkeypatch):
    root = tmp_path / "home"
    root.mkdir()
    monkeypatch.setattr(Path, "home", lambda: root)
    sessions = root / ".delfin" / "agent_sessions"
    sessions.mkdir(parents=True)
    monkeypatch.setattr(ws, "build_working_state_block",
                        lambda *a, **k: _BLOCK)
    with patch("delfin.agent.session_store._SESSIONS_DIR", sessions):
        yield root


def _terminal_compaction(tmp_path):
    engine, client = _engine(tmp_path)
    before = _history(30)
    engine.messages = [dict(m) for m in before]
    return engine, client, before


def _run(engine, client):
    R.run_turn(engine, "go on", sink=lambda item: None)
    assert client.summary_calls == 1, "the compaction ran"
    assert R._FRESH_CONTEXT_MARKER not in str(client.sent), "no fresh start"
    first = engine.messages[0]["content"]
    assert first.startswith(ws._HEADER), first[:80]
    assert "\n[Conversation summary — older messages compacted]\n" in first


def test_a_save_without_the_record_keeps_the_compacted_history(home, tmp_path):
    engine, client, before = _terminal_compaction(tmp_path)
    ss.save_session(SID, mode="solo", chat_messages=_chat(60),
                    engine_messages=list(engine.messages))
    _run(engine, client)

    # The display transcript keeps growing; no compaction record passed.
    ss.save_session(SID, mode="solo", chat_messages=_chat(62),
                    engine_messages=list(engine.messages))
    loaded = ss.load_session(SID)["engine_messages"]

    assert _contents(loaded) == _contents(engine.messages)
    assert not list((home / ".delfin" / "agent_sessions").glob(
        f"{SID}.json.history-kept-*"))
    archived = set(_contents(before[:-engine._KEEP_RECENT]))
    assert not archived & set(_contents(loaded)), (
        "an archived row came back into the saved history")


def test_the_headless_saver_writes_the_compaction_record(home, tmp_path):
    from delfin.agent import cli
    engine, client, _before = _terminal_compaction(tmp_path)
    _run(engine, client)
    cli._save_session(engine, tmp_path)
    saved = ss.load_session(SID)
    assert saved["last_compaction_info"] == engine.last_compaction_info
    assert saved["last_compaction_info"]["kind"] == "summary"


def test_a_resume_reads_the_compaction_record_back(home, tmp_path):
    """Saved but not restored, the record would differ from the resumed
    engine's empty one at the first save, and that difference alone counts
    as an explanation -- any shorter save would then pass the guard."""
    from delfin.agent import cli
    engine, client, _before = _terminal_compaction(tmp_path)
    _run(engine, client)
    cli._save_session(engine, tmp_path)

    resumed, _client = _engine(tmp_path / "second")
    assert resumed.last_compaction_info is None
    resumed.restore_state(ss.load_session(SID))
    assert resumed.last_compaction_info == engine.last_compaction_info


@pytest.mark.parametrize("first, explained", [
    # What the compaction writes: the summary alone, or behind the block.
    ("[Conversation summary — older messages compacted]\nGoal", True),
    (_BLOCK + "\n[Conversation summary — older messages compacted]\nGoal",
     True),
    # The in-turn compaction's own row, accepted at the head as before.
    ("[Conversation summary — earlier rounds of this turn compacted]\nx",
     True),
    # A block with no summary behind it, the marker quoted mid-sentence,
    # or the marker far beyond any block the compaction can write.
    (_BLOCK + "\nno summary here", False),
    ("see [Conversation summary — older messages compacted]\nGoal", False),
    (_BLOCK + "y" * 5000
     + "\n[Conversation summary — older messages compacted]\nGoal", False),
    ("[Conversation summary of my own]\nGoal", False),
], ids=["summary", "block-then-summary", "in-turn-summary", "block-only",
        "quoted-marker", "marker-beyond-a-block", "other-bracket-text"])
def test_the_markers_the_guard_accepts(first, explained):
    assert ss._begins_with_compaction_summary(first) is explained


def test_the_markers_the_guard_reads_are_the_ones_the_compaction_writes():
    import inspect
    from delfin.agent import in_turn_compaction as itc
    assert ss._WORKING_STATE_HEADER == ws._HEADER
    assert ss._IN_TURN_SUMMARY_MARKER == itc._SUMMARY_HEADER
    assert ss._SUMMARY_MARKER in inspect.getsource(
        AgentEngine._compact_history)
