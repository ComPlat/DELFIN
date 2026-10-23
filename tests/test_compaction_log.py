"""The compaction log: one durable record per compaction event.

Control test (red on the commit that introduced compaction_log.py without
the engine call): ``_compact_history`` must write a record through the
PUBLIC write path — the same call real operation makes — for every one of
the shapes ``last_compaction_info`` can take.
"""

from __future__ import annotations

from unittest.mock import MagicMock

from delfin.agent import compaction_log
from delfin.agent.engine import AgentEngine


SESSION = "ctl-compaction-log"


def _bare_engine(tmp_path, monkeypatch):
    eng = AgentEngine.__new__(AgentEngine)
    eng.messages = []
    eng.role_outputs = {}
    eng.compaction_summaries = {}
    eng.token_usage = {"input": 0, "output": 0}
    eng.cost_usd = 0.0
    eng.context_window_tokens = 100_000
    eng.auto_compact_pct = 0.80
    eng.last_compaction_info = None
    eng.session_id = SESSION
    eng.backend = "api"
    eng.client = MagicMock()
    eng.current_role_index = 0
    eng.route = ["solo_agent"]
    eng.repo_dir = tmp_path
    return eng


def _big_history():
    big = "y" * 25_000
    return [
        {"role": "user", "content": f"msg-{i} " + big}
        for i in range(20)
    ]


class TestEngineWritesTheLog:
    def test_full_compaction_leaves_one_record(self, tmp_path, monkeypatch):
        eng = _bare_engine(tmp_path, monkeypatch)
        eng.messages = _big_history()
        assert eng._should_auto_compact()
        eng._compact_history()
        records = compaction_log.read_compactions(SESSION)
        assert len(records) == 1
        rec = records[0]
        assert rec["kind"] in ("summary", "deterministic_digest")
        assert rec["messages_compacted"] >= 1
        assert rec["tokens_before"] > rec["tokens_after"]
        assert rec["session"] == SESSION

    def test_record_names_the_working_state_sections(self, tmp_path,
                                                     monkeypatch):
        eng = _bare_engine(tmp_path, monkeypatch)
        # A machine turn with a test outcome so the working-state block
        # has at least one section to carry.
        hist = _big_history()
        hist[0] = {"role": "user", "content":
                   "[Command results]\n1 passed"}
        eng.messages = hist
        eng._compact_history()
        rec = compaction_log.read_compactions(SESSION)[0]
        assert isinstance(rec["state_block_sections"], dict)

    def test_forced_compaction_is_flagged(self, tmp_path, monkeypatch):
        eng = _bare_engine(tmp_path, monkeypatch)
        eng.messages = _big_history()
        eng._compact_history(force=True)
        rec = compaction_log.read_compactions(SESSION)[0]
        assert rec["forced"] is True

    def test_record_carries_the_loss_comparison(self, tmp_path, monkeypatch):
        eng = _bare_engine(tmp_path, monkeypatch)
        hist = _big_history()
        # A machine turn whose facts the block may or may not carry —
        # the record must NAME what is lost either way.
        hist[0] = {"role": "user", "content":
                   "[Command results]\ngate tests/test_zz.py -> 2 passed"}
        eng.messages = hist
        eng._compact_history()
        rec = compaction_log.read_compactions(SESSION)[0]
        assert "lost" in rec
        assert set(rec["lost"]) == {
            "denials", "tests", "files", "instructions"}


class TestNoLogBeforeCompaction:
    def test_short_history_writes_nothing(self, tmp_path, monkeypatch):
        eng = _bare_engine(tmp_path, monkeypatch)
        eng.messages = [
            {"role": "user", "content": "small"},
            {"role": "assistant", "content": "ok"},
        ]
        eng._compact_history()
        assert compaction_log.read_compactions(SESSION) == []


class TestRecordShape:
    def test_sections_are_counts_not_contents(self):
        block = (
            "[Working state]\n"
            "Last test outcomes:\n"
            "  2 passed\n"
            "\n"
            "Recently worked on:\n"
            "  delfin/agent/engine.py\n"
        )
        secs = compaction_log._state_block_sections(block)
        assert secs["Last test outcomes:"] == 1
        assert secs["Recently worked on:"] == 1
        # Headings only — the log must never carry the block's contents.
        assert "2 passed" not in str(secs)

    def test_read_is_empty_for_unknown_session(self):
        assert compaction_log.read_compactions("no-such-session-xyz") == []
