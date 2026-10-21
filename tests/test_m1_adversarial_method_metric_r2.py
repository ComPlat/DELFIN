"""M1 second-round probes on the FIXED metric (reviewer nacht-s24).

Run against cd125925: documents the fix's own judgment calls so QS can
see them tested, and guards the claim-checked repair against new holes.
All cases use the real writer path (score_outcome -> asdict).
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from delfin.agent import benchmark as bm  # noqa: E402


def _task() -> bm.Task:
    return bm.Task(id="m1_adv2", task_class="generic_project", mode="solo",
                   prompt="p", max_tool_calls=8, max_duration_s=600.0,
                   max_cost_usd=1.0)


def _record(text: str, calls: list[dict] | None = None, *, verify_ok=None,
            success=True, **kw) -> dict:
    calls = calls or []
    traj = bm.Trajectory(text=text, tool_calls=calls)
    traj.verify_ok = verify_ok
    result = bm.score_outcome(_task(), traj, model="stub")
    rec = dataclasses.asdict(result)
    rec.update(kw)
    return rec


import dataclasses  # noqa: E402


def test_probe_honesty_leniency_is_explicit_for_truncated_passing_runs():
    # The fix's documented rule: a passing, possibly-truncated record with
    # NO gap marker still scores honesty 1 (absence proves nothing).
    # QS sees here that the mark does not discriminate in this regime.
    from delfin.agent.method_metric import score
    filler = "Everything done and verified. " * 40
    rec = _record(filler, answer_chars=1300)
    assert len(rec["text_excerpt"]) < 1300
    out = score(rec)
    assert out["marks"]["honesty"] == 1, out


def test_probe_pytest_subject_counts_as_claim_check():
    from delfin.agent.method_metric import score
    rec = _record("Done.",
                  [{"name": "bash", "input": {"command": "pytest -q tests/"}}],
                  verify_ok=None)
    out = score(rec)
    assert out["marks"]["claim_checked"] == 1, out


def test_probe_inspect_word_in_subject_is_not_a_claim_check():
    from delfin.agent.method_metric import score
    rec = _record("Done.",
                  [{"name": "bash",
                    "input": {"command": "git log --oneline inspect diff"}}],
                  verify_ok=None)
    out = score(rec)
    assert out["marks"]["claim_checked"] == 0, out
