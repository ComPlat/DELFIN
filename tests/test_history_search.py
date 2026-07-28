"""Tests for delfin.agent.history_search — model-callable recall over
the live message list + archived pre-compaction transcripts.

Archives are built through the real writer
(session_store.archive_pre_compaction_transcript) against a fake home,
so the reader is exercised against the exact on-disk format.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import history_search as hs
from delfin.agent import session_store as ss


@pytest.fixture
def fake_home(monkeypatch, tmp_path):
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    monkeypatch.setattr(ss, "_SESSIONS_DIR", tmp_path / ".delfin" / "agent_sessions")
    return tmp_path


def _archive_file(sid: str) -> Path:
    return Path.home() / ".delfin" / "transcript_archive" / f"{sid}.jsonl"


# ---------------------------------------------------------------------------
# history_search
# ---------------------------------------------------------------------------

def test_live_and_archive_hits_merged_and_ranked(fake_home):
    ss.archive_pre_compaction_transcript(
        "sid-merge",
        [
            {"role": "user", "content": "please compute the singlet fission "
                                        "rate for tetracene with ORCA"},
            {"role": "assistant", "content": "totally unrelated small talk"},
        ],
        info={"messages_compacted": 2},
    )
    live = [
        {"role": "user", "content": "back to the singlet fission topic"},
        {"role": "assistant", "content": "nothing relevant here"},
    ]
    hits = hs.history_search("sid-merge", "singlet fission", messages=live)
    real = [h for h in hits if "score" in h]
    assert real, "expected hits from both sources"
    sources = {h["source"] for h in real}
    assert sources == {"live", "archive"}
    # sorted by score descending
    scores = [h["score"] for h in real]
    assert scores == sorted(scores, reverse=True)
    # every hit carries a resolvable ref + snippet around the match
    for h in real:
        assert h["ref"].startswith(("live:", "archive:"))
        assert "fission" in h["snippet"].lower()
    # the unrelated messages never match
    assert all("unrelated" not in h["snippet"] for h in real)
    assert all("nothing relevant" not in h["snippet"] for h in real)


def test_ranking_ties_break_newest_first(fake_home):
    # Two identical texts, one archived (older) and one live (newer):
    # deterministic ordering must put the live (newest) one first.
    ss.archive_pre_compaction_transcript(
        "sid-tie", [{"role": "user", "content": "benchmark the polariton model"}]
    )
    live = [{"role": "user", "content": "benchmark the polariton model"}]
    hits = hs.history_search("sid-tie", "polariton benchmark", messages=live)
    real = [h for h in hits if "score" in h]
    assert len(real) == 2
    assert real[0]["source"] == "live"
    assert real[1]["source"] == "archive"


def test_substring_fallback_on_no_query_tokens(fake_home):
    # "S1" tokenises to nothing (len < 3) -> substring matching must kick in.
    ss.archive_pre_compaction_transcript(
        "sid-sub", [{"role": "assistant", "content": "the S1 state lies at 2.1 eV"}]
    )
    hits = hs.history_search("sid-sub", "S1", messages=[
        {"role": "user", "content": "no such token here"},
    ])
    real = [h for h in hits if "score" in h]
    assert len(real) == 1
    assert real[0]["source"] == "archive"
    assert "S1 state" in real[0]["snippet"]
    assert real[0]["score"] == 1.0


def test_archive_only_session(fake_home):
    ss.archive_pre_compaction_transcript(
        "sid-arch", [{"role": "user", "content": "record the dielectric constant"}]
    )
    hits = hs.history_search("sid-arch", "dielectric constant")   # messages=None
    real = [h for h in hits if "score" in h]
    assert real and all(h["source"] == "archive" for h in real)
    assert real[0]["ts"] is not None   # compacted_at flows through


def test_live_only_session_without_archive(fake_home):
    hits = hs.history_search("sid-noarch", "gradient descent", messages=[
        {"role": "user", "content": "explain gradient descent please"},
    ])
    real = [h for h in hits if "score" in h]
    assert len(real) == 1
    assert real[0]["source"] == "live"
    assert real[0]["ref"] == "live:0"


def test_corrupt_jsonl_lines_skipped(fake_home):
    sid = "sid-corrupt"
    ss.archive_pre_compaction_transcript(
        sid, [{"role": "user", "content": "first block about ruthenium"}]
    )
    # Corrupt the file the way a crash mid-append would: torn JSON,
    # a non-object line, and a blank line.
    with _archive_file(sid).open("a", encoding="utf-8") as fh:
        fh.write('{"compacted_at": 1, "messages": [{"role"\n')
        fh.write("42\n")
        fh.write("\n")
    ss.archive_pre_compaction_transcript(
        sid, [{"role": "user", "content": "second block about ruthenium"}]
    )
    hits = hs.history_search(sid, "ruthenium")
    real = [h for h in hits if "score" in h]
    assert len(real) == 2   # both valid records found, garbage skipped, no raise
    # Refs use raw line numbers, so the post-garbage record resolves too.
    newest = [h for h in real if "second" in h["snippet"]][0]
    got = hs.history_get(sid, newest["ref"])
    assert got.get("text") == "second block about ruthenium"


def test_capped_scan_note_and_old_messages_skipped(fake_home, monkeypatch):
    monkeypatch.setattr(hs, "SCAN_CAP", 5)
    sid = "sid-cap"
    msgs = [{"role": "user", "content": "oldest zirconium mention"}]
    msgs += [{"role": "user", "content": f"filler number {i}"} for i in range(9)]
    ss.archive_pre_compaction_transcript(sid, msgs)
    hits = hs.history_search(sid, "zirconium")
    # The zirconium message fell outside the most-recent-5 window ...
    assert not [h for h in hits if "score" in h]
    # ... and the result says so.
    assert hits and hits[-1].get("capped") is True
    assert "capped" in hits[-1]["note"]
    # Within the window everything still works, note still appended.
    hits2 = hs.history_search(sid, "", messages=None)
    real2 = [h for h in hits2 if "score" in h]
    assert 0 < len(real2) <= 5
    assert hits2[-1].get("capped") is True


def test_search_never_raises_on_invalid_inputs(fake_home):
    assert hs.history_search("sid-x", "", messages=None) == []
    # None query + non-string content: still a list, no exception
    got = hs.history_search("sid-x", None, messages=[{"role": "u", "content": 5}])
    assert isinstance(got, list)
    # messages entries that are not dicts are skipped, not fatal
    hits = hs.history_search("sid-x", "abcdefu", messages=["not-a-dict", 7])
    assert hits == []


# ---------------------------------------------------------------------------
# history_get
# ---------------------------------------------------------------------------

def test_history_get_returns_full_text_and_truncates(fake_home):
    long_text = ("H" * 3000) + " MIDDLE-MARKER-XYZ " + ("T" * 3000)
    ss.archive_pre_compaction_transcript(
        "sid-long", [{"role": "assistant", "content": long_text}]
    )
    # Untruncated when it fits
    full = hs.history_get("sid-long", "archive:0:0", max_chars=10_000)
    assert full["truncated"] is False
    assert full["text"] == long_text
    assert full["total_chars"] == len(long_text)
    # Head+tail with marker when it does not
    cut = hs.history_get("sid-long", "archive:0:0", max_chars=1000)
    assert cut["truncated"] is True
    assert "[truncated:" in cut["text"]
    assert cut["text"].startswith("HHH")
    assert cut["text"].endswith("TTT")
    assert len(cut["text"]) < 1200   # ~max_chars + marker


def test_history_get_live_ref(fake_home):
    live = [
        {"role": "user", "content": "keep the basis set def2-TZVP"},
        {"role": "assistant", "content": "noted"},
    ]
    got = hs.history_get("sid-live", "live:0", messages=live)
    assert got["source"] == "live"
    assert got["role"] == "user"
    assert got["text"] == "keep the basis set def2-TZVP"
    # live ref without a live list must fail soft, not crash
    no_live = hs.history_get("sid-live", "live:0")
    assert "error" in no_live


def test_history_get_bad_refs(fake_home):
    ss.archive_pre_compaction_transcript(
        "sid-refs", [{"role": "user", "content": "hello"}]
    )
    assert "error" in hs.history_get("sid-refs", "bogus")
    assert "error" in hs.history_get("sid-refs", "archive:5:0")     # no such record
    assert "error" in hs.history_get("sid-refs", "archive:0:9")     # index range
    assert "error" in hs.history_get("sid-refs", "live:-1", messages=[])
    assert "error" in hs.history_get("sid-nothing", "archive:0:0")  # no archive


# ---------------------------------------------------------------------------
# path safety
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("bad_sid", [
    "../evil", "..", "a/../b", "a/b", "a\\b", ".hidden", "", "x" * 200,
])
def test_traversal_session_ids_rejected(fake_home, bad_sid):
    hits = hs.history_search(bad_sid, "anything", messages=[
        {"role": "user", "content": "anything"},
    ])
    assert hits and "error" in hits[0]
    got = hs.history_get(bad_sid, "archive:0:0")
    assert "error" in got
    # And nothing escaped the archive dir
    outside = fake_home.parent / "evil.jsonl"
    assert not outside.exists()
