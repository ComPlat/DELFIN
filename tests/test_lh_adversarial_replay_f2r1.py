"""F2 P2a reviewer adversarial tests — the replay runner (delfin/agent/longhorizon/replay.py).

The claim under attack (replay.py docstring + test_lh_replay.py): a
transcript event is "a record the matching real writer accepts", and a
"garbage event is skipped, never raises". Each test here plays a
transcript that a real user of replay would plausibly produce and checks
the promise the module makes about it. Written against 21a07d812.
"""
from __future__ import annotations

import pytest


@pytest.fixture
def redirected_store(tmp_path, monkeypatch):
    from delfin.agent import compaction_log, tool_trace, turn_metrics
    monkeypatch.setattr(turn_metrics, "_DIR", tmp_path / "turn_metrics")
    monkeypatch.setattr(tool_trace, "_DIR", tmp_path / "tool_traces")
    monkeypatch.setattr(compaction_log, "_LOG_DIR", tmp_path / "compaction_log")
    return tmp_path


def test_a_malformed_turn_event_is_skipped_never_raises(redirected_store):
    """An event carrying a key the real writer does not accept must be
    skipped, not crash the run.

    _emit passes ``**ev`` straight into the writer call (replay.py:43) —
    an unexpected keyword raises TypeError at CALL BINDING, before the
    writer's own never-raises try/except (turn_metrics.py:136) can run.
    test_lh_replay.py names this promise ("a malformed event is skipped,
    never raises") but tests only an unknown KIND; a malformed TURN event
    is the case a transcript authoring bug produces.
    """
    from delfin.agent.longhorizon import replay
    r = replay.replay([
        {"kind": "turn", "provider": "kit", "model": "m", "total_ms": 1,
         "bogus_key": 1},
        {"kind": "turn", "provider": "kit", "model": "m", "total_ms": 2,
         "output_chars": 1, "output_tokens": 1},
    ], session="s-malformed")
    assert r["metrics"]["turns"] == 1


def test_round_field_survives_the_replay(redirected_store):
    """tool_trace.record accepts ``round`` (tool_trace.py:82); the real
    callers record it. The transcript's contract says a tool event is a
    record the real writer accepts — replay must not silently drop a
    field the writer keeps (replay.py:47-55 passes explicit kwargs and
    omits it)."""
    from delfin.agent.longhorizon import replay
    r = replay.replay([
        {"kind": "tool", "tool": "bash", "input": "ls", "output": "x",
         "round": 7},
    ], session="s-round")
    recs = r["records"]["tools"]
    assert recs, "tool record lost entirely"
    assert recs[0].get("round") == 7


def test_a_case_typed_kind_is_not_silently_a_zero(redirected_store):
    """{"kind": "Turn"} (capitalised, a typo) is ignored by the dispatch
    (replay.py:64-65): the whole transcript reduces to turns=0 with no
    signal at all. A run whose only defect is a typo must be
    DIAGNOSABLE, not silently zero — the reducer reports None rather
    than a guessed 0, and replay must surface what it dropped."""
    from delfin.agent.longhorizon import replay
    r = replay.replay([
        {"kind": "Turn", "provider": "kit", "model": "m", "total_ms": 1,
         "output_chars": 1, "output_tokens": 1},
    ], session="s-typo")
    m = r["metrics"]
    assert m["turns"] == 0, "a typo'd kind must not be counted as a turn"
    # and the drop must be visible: replay reports what it did not
    # recognise instead of silently discarding it.
    assert r.get("unrecognized_events", 0) >= 1
