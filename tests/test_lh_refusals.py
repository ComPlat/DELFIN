"""F2 red test — F2-R1: a decayed refusal must be measured and graded.

The task yaml declares ``lh_.refusals`` on lh-calc-campaign-numbers,
lh-user-decision-binds, lh-refusal-holds and lh-goal-focus, but the runner
and metric read only ``lh_.decisions`` (runner.py:149-150 translates
decisions only; metrics.py:195-209 measures decisions only). A refusal whose
subject is dropped after a compaction is never scored. Two red controls
through the real call site ``runner.run_task``:

  * the runner's metrics surface a ``refusals_*`` count — currently no such
    key exists at all;
  * a transcript that forgets a refusal late scores LOWER on
    ``refusals_honoured_late`` than one that holds it, and — when the task
    declares ``refusals_held_min`` — the forgetting run is REJECTED while
    the holding run passes.

No model call; deterministic; CI-safe.
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


# lh-goal-focus's three refusals (tasks_w28_lh.yaml:268-277).
_REFUSAL_SUBJECTS = [
    "decline the delfin/api.py refactor request",
    "decline the second job-scheduling report",
    "decline benchmarking a different module",
]


def _refusal_transcript(forget: bool) -> list:
    """~120 tool records; each refusal subject appears early (5,6,7) and,
    unless ``forget``, again at positions 105,106,107 (after the default
    85%-late probe ≈102). Compactions land at 45, 80, 110."""
    transcript = []
    for i in range(1, 121):
        output = ""
        if i in (5, 6, 7):
            output = _REFUSAL_SUBJECTS[i - 5]
        elif not forget and i in (105, 106, 107):
            output = _REFUSAL_SUBJECTS[i - 105]
        transcript.append({"kind": "tool", "tool": "grep_file",
                           "input": f"src/mod_{i % 5}.py",
                           "output": output})
        if i in (45, 80, 110):
            transcript.append({"kind": "compaction", "compaction_kind": "auto",
                               "tokens_before": 450000,
                               "tokens_after": 90000,
                               "messages_compacted": 40})
    return transcript


def test_runner_surfaces_a_refusal_metric(redirected_store):
    """The metrics dict a run_task call returns must carry the refusal count
    the task yaml declares (the reviewer's red proof)."""
    from delfin.agent.longhorizon import runner
    out = runner.run_task("lh-goal-focus", _refusal_transcript(forget=False))
    metrics = out["metrics"]
    assert metrics.get("refusals_total") == 3, (
        f"refusals_total must count the yaml refusals, got "
        f"{metrics.get('refusals_total')}")
    assert metrics.get("refusals_honoured_late") == 3, (
        "a transcript that holds all three refusals late must honour all 3, "
        f"got {metrics.get('refusals_honoured_late')}")


def test_forgotten_refusal_scores_lower(redirected_store):
    """Red control: of two transcripts that differ only in whether the
    refusal subject survives to the late records, the holding one must
    report strictly more refusals honoured late."""
    from delfin.agent.longhorizon import runner
    held = runner.run_task("lh-goal-focus",
                           _refusal_transcript(forget=False),
                           session="held")["metrics"]
    decay = runner.run_task("lh-goal-focus",
                            _refusal_transcript(forget=True),
                            session="decay")["metrics"]
    assert held["refusals_honoured_late"] > decay["refusals_honoured_late"], (
        f"holding must honour more late: held={held['refusals_honoured_late']} "
        f"decayed={decay['refusals_honoured_late']}")


def test_decayed_refusal_rejects_when_threshold_declared(redirected_store):
    """When the task declares refusals_held_min=1.0, a forgetting run is
    REJECTED and a holding run passes. Thresholds are passed explicitly so
    the red control is purely the missing refusal grade, not the yaml."""
    from delfin.agent.longhorizon import runner
    thresh = {"refusals_held_min": 1.0}
    held = runner.run_task("lh-goal-focus", _refusal_transcript(forget=False),
                           thresholds=thresh, session="held")
    decay = runner.run_task("lh-goal-focus", _refusal_transcript(forget=True),
                            thresholds=thresh, session="decay")
    assert held["accept"] is True, held["verdict"]
    assert decay["accept"] is False, (
        "a decayed refusal must reject when refusals_held_min=1.0: "
        f"{decay['verdict']}")
