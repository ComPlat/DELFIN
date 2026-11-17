"""F2 P2 red test — the long-horizon replay runner (deterministic, CI-safe).

Target: delfin/agent/longhorizon/replay.py (new, owned). ``replay`` plays a
scripted long-horizon transcript through the REAL record writers into a
redirected store and reduces it to metrics. The red control at TRANSCRIPT
level (mode a, brief item 3): two transcripts of the same task that differ
only in whether the planted fact survives to the LATE records must score
different — the recalling one higher on late recall, the forgetting one
lower. No model call anywhere; deterministic, runs in CI.
"""
from __future__ import annotations

import pytest


@pytest.fixture
def redirected_store(tmp_path, monkeypatch):
    """Redirect the real record writers into a per-test store so replay
    writes and reads the same scratch dir (never the real ~/.delfin)."""
    from delfin.agent import compaction_log, tool_trace, turn_metrics
    monkeypatch.setattr(turn_metrics, "_DIR", tmp_path / "turn_metrics")
    monkeypatch.setattr(tool_trace, "_DIR", tmp_path / "tool_traces")
    monkeypatch.setattr(compaction_log, "_LOG_DIR", tmp_path / "compaction_log")
    return tmp_path


def _long_transcript(forget: bool) -> list:
    """A ~120-record long-horizon transcript of the SAME task (a multi-file
    code change; planted fact: the value -37.460 given at turn 3, probed
    late). Several compactions land along the way. ``forget`` drops the fact
    from every record but the early one (turn 3), so nothing late carries it.

    Record count matters: a '100+ turn' task must survive several
    compactions, not one snapshot — so the transcript carries three
    compactions and ~120 records.
    """
    transcript = []
    for i in range(1, 121):
        transcript.append({"kind": "turn", "provider": "kit",
                           "model": "kit.deepseek-v4-flash",
                           "role": "assistant", "total_ms": 5000,
                           "output_chars": 60, "tool_calls": 1,
                           "input_tokens": 200, "output_tokens": 80,
                           "cached_tokens": 60})
        # the planted fact: an early record always carries it; the late
        # periodic records (turns 53, 103) carry it only when NOT forgetting.
        output = ""
        if i == 3 or (not forget and i % 50 == 3):
            output = "E = -37.460"
        transcript.append({"kind": "tool", "tool": "read_file",
                           "input": f"src/mod_{i%5}.py",
                           "output": output})
        if i in (45, 80, 110):
            transcript.append({"kind": "compaction", "compaction_kind": "auto",
                               "tokens_before": 450000,
                               "tokens_after": 90000,
                               "messages_compacted": 40})
    return transcript


def test_recalling_transcript_scores_higher_late_recall(redirected_store):
    """The red control: of two identical transcripts that differ only in
    whether the exact planted number survives to the late records, the
    recalling one must report a higher late recall."""
    from delfin.agent.longhorizon import replay

    facts = [{"fact_id": "mass", "value": "-37.460", "probe": 5}]

    r_recall = replay.replay(_long_transcript(forget=False),
                             session="s-rec", planted_facts=facts)
    r_forget = replay.replay(_long_transcript(forget=True),
                             session="s-for", planted_facts=facts)

    rec = r_recall["metrics"]["facts_recalled_late"] / max(
        r_recall["metrics"]["facts_total"], 1)
    forg = r_forget["metrics"]["facts_recalled_late"] / max(
        r_forget["metrics"]["facts_total"], 1)
    assert rec > forg, (
        f"recalling must score higher: recall={rec} forget={forg}")
    # the session really spanned a long horizon with several compactions
    assert r_recall["metrics"]["turns"] >= 100
    assert r_recall["metrics"]["compactions"] >= 3


def test_replay_records_are_real_writer_shapes(redirected_store):
    """The replay wrote through the REAL writers: the records the reducer
    read carry the real fields (ts, model, tokens), not my dict shape."""
    from delfin.agent.longhorizon import replay
    r = replay.replay([
        {"kind": "turn", "provider": "kit", "model": "m",
         "total_ms": 100, "output_chars": 5, "output_tokens": 1},
        {"kind": "tool", "tool": "read_file", "input": "a.txt",
         "output": "x"},
    ], session="s-shape")
    turns = r["records"]["turns"]
    tools = r["records"]["tools"]
    assert turns and "ts" in turns[0] and turns[0]["model"] == "m"
    assert tools and "ts" in tools[0] and tools[0]["tool"] == "read_file"
    assert r["metrics"]["turns"] == 1


def test_a_garbage_event_never_raises(redirected_store):
    """An unknown record kind or a malformed event is skipped, never raises
    (the writers' own never-raises contract)."""
    from delfin.agent.longhorizon import replay
    r = replay.replay([
        {"kind": "nonsense", "whatever": 1},
        {"kind": "turn", "provider": "kit", "model": "m",
         "total_ms": 1, "output_chars": 1, "output_tokens": 1},
    ], session="s-garbage")
    assert r["metrics"]["turns"] == 1
