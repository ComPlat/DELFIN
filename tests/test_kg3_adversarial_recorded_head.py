"""Adversarial contract: the request head measured on REAL engine rows.

KG3 reviewer test, guardrail 4 of phase 3. A real AgentEngine
(throwaway pack, fake client) runs two consecutive turns through the
real repl turn path; an instance-level recorder around
engine.loader.compose_sections captures the section rows the engine
actually built. prompt_head.stable_prefix_len must report a real,
non-zero shared prefix on those recorded PromptSection rows — a module
that reads only synthetic dict/tuple shapes fails here (finding 3,
red on 267a5818, green since 135d653b).
"""

from unittest.mock import patch

from delfin.agent.engine import AgentEngine
from delfin.agent import prompt_head as ph

from tests.test_a_compaction_before_a_fresh_start_saves_and_stops import (
    _engine, _pack_tree, _Client)


def _recorded_rows(tmp_path):
    engine, client = _engine(tmp_path, how="answer")

    recorded = []
    real_compose = engine.loader.compose_sections

    def recording_compose(*args, **kwargs):
        rows = real_compose(*args, **kwargs)
        recorded.append(list(rows))
        return rows

    engine.loader.compose_sections = recording_compose

    with patch.object(type(client), "stream_message", client.stream_message):
        first = list(engine.stream_response("Turn one: measure the cache."))
        assert first, "turn one produced nothing"
        second = list(engine.stream_response("Turn two: measure the cache."))
        assert second, "turn two produced nothing"

    assert len(recorded) >= 2, (
        f"expected compose_sections for both turns, got {len(recorded)}")
    return recorded


def test_recorded_turns_share_a_stable_prefix(tmp_path):
    recorded = _recorded_rows(tmp_path)

    turn_a, turn_b = recorded[-2], recorded[-1]
    report = ph.stable_prefix_len(turn_a, turn_b)
    print(f"turn_a rows: {len(turn_a)}; turn_b rows: {len(turn_b)}")
    for r in turn_b:
        print(f"  row: {r!r}")
    assert report["bytes"] > 0, (
        f"prompt_head cannot read real PromptSection rows: {report}")
    # No STABLE (layer-0) row after the shared prefix: that would make
    # the tail uncacheable by design.
    assert report["rows_all_stable"], (
        f"a stable row sits after the shared prefix: {report}")
    # The safety ordering holds on the rows the engine really sent.
    problems = ph.check_ordering(turn_b)
    assert problems == [], f"ordering broken on recorded rows: {problems}"
