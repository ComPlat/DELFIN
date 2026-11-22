"""The prompt_head mount is on the ENGINE's turn path.

The three prompt_head functions were built and tested (KG3) but nothing
in the running product ever called them. The mount is the engine tick:
_build_current_system_prompt -> compose_sections -> _prompt_head_check
for every real turn. This test drives the real engine turn path (the
same recorded-rows harness as test_kg3_adversarial_recorded_head) and
asks that the checks RAN: findings recorded on the engine, the
cross-turn prefix report kept, and the rows checked are the ones the
engine actually built — not a helper called with synthetic rows.
"""

from unittest.mock import patch

from delfin.agent.engine import AgentEngine
from delfin.agent import prompt_head as ph

from tests.test_a_compaction_before_a_fresh_start_saves_and_stops import (
    _engine, _pack_tree, _Client)


def test_turn_path_runs_prompt_head_checks(tmp_path):
    engine, client = _engine(tmp_path, how="answer")

    recorded = []
    real_compose = engine.loader.compose_sections

    def recording_compose(*args, **kwargs):
        rows = real_compose(*args, **kwargs)
        recorded.append(list(rows))
        return rows

    engine.loader.compose_sections = recording_compose

    with patch.object(type(client), "stream_message", client.stream_message):
        first = list(engine.stream_response("Turn one: check the head."))
        assert first, "turn one produced nothing"
        second = list(engine.stream_response("Turn two: check the head."))
        assert second, "turn two produced nothing"

    assert len(recorded) >= 2, (
        f"expected compose_sections for both turns, got {len(recorded)}")

    # The helper ran on the rows the engine actually built: the recorded
    # last turn equals what the helper kept.
    assert engine._prompt_head_last_sections == recorded[-1], (
        "the mount ran on rows other than the ones compose_sections built")
    assert getattr(engine, "_prompt_head_findings", None) == [], (
        f"turn path produced findings on a healthy head: "
        f"{engine._prompt_head_findings}")

    # The cross-turn report is the REAL measurement over the two turns,
    # matching what prompt_head computes from the recorded rows.
    report = engine._prompt_head_prefix
    expected = ph.stable_prefix_len(recorded[-2], recorded[-1])
    assert report == expected, (
        f"prefix report on the engine {report} != measured {expected}")
    assert report["bytes"] > 0, (
        f"no shared prefix measured on real engine turns: {report}")
    assert report["rows_all_stable"], (
        f"a stable row sits after the shared prefix: {report}")
