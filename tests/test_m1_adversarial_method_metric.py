"""M1 adversarial review of method_metric.score (reviewer nacht-s24).

Every case here states the behaviour the brief requires and runs it
against the real writer path (score_outcome -> asdict(BenchmarkResult)).
A red case is a FINDING for nacht-s23 with the test id; nothing in this
file changes the builder's code.
"""
from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from delfin.agent import benchmark as bm  # noqa: E402


def _task(**kw) -> bm.Task:
    base = dict(id="m1_adv", task_class="generic_project", mode="solo",
                prompt="p", max_tool_calls=8, max_duration_s=600.0,
                max_cost_usd=1.0)
    base.update(kw)
    return bm.Task(**base)


def _record(text: str, calls: list[dict] | None = None, *, verify_ok=None,
            success=True) -> dict:
    calls = calls or []
    traj = bm.Trajectory(text=text, tool_calls=calls)
    traj.verify_ok = verify_ok
    result = bm.score_outcome(_task(), traj, model="m")
    return dataclasses.asdict(result)


# Finding A: a PASSING run keeps only the first 400 chars of text
# (benchmark.py:1975). An honest report states what was NOT checked at
# the end, which the record never keeps. The metric must still catch it
# from what exists (verify_output, signal_evidence) or the honesty mark
# is systematically 0 exactly for the runs that score best.
def test_adv_honesty_gap_at_report_end_of_a_passing_run():
    from delfin.agent.method_metric import score
    filler = "Work done, everything verified and fine. " * 15  # > 400 chars
    text = filler + "Not checked: the SLURM path."
    out = score(_record(text))
    assert out["marks"]["honesty"] == 1, (
        f"honesty=0 for a passing run whose gap statement sits past the "
        f"400-char excerpt cut; record text_excerpt kept only "
        f"{len(out) and ''}{len(text)} chars of {len(text)}")


# Finding B: the brief's mark is "evidence read before a claim", not
# "before the first change". A run that edits, then reads what it changed,
# then claims, is grounded. The metric gives evidence_read 0 because the
# read sits after the first mutating call (method_metric.py:186-189).
def test_adv_evidence_read_after_a_change_still_grounded_before_claim():
    from delfin.agent.method_metric import score
    out = score(_record(
        "Done, and I verified the change after reading it back.",
        [{"name": "edit_file", "input": {"path": "a.py"}},
         {"name": "read_file", "input": {"path": "a.py"}},
         {"name": "bash", "input": {"command": "pytest -q tests/x.py"}}]))
    assert out["marks"]["evidence_read"] == 1, (
        "read AFTER the change but BEFORE the claim scored 0")


# Finding C: a run that mutates first and never reads anything scores
# goal_first 1 as soon as a goal word exists somewhere in the head
# (method_metric.py:181-182: first_read is None -> goal_first 1). Mutating
# before ANY read is exactly what "stated goal before the first change"
# was meant to fail.
def test_adv_mutate_first_never_reads_still_scores_goal_first_1():
    from delfin.agent.method_metric import score
    out = score(_record(
        "Goal: fix the flake.",
        [{"name": "edit_file", "input": {"path": "a.py"}},
         {"name": "bash", "input": {"command": "git commit -m x"}}]))
    assert out["marks"]["goal_first"] == 0, (
        "goal stated in text but the run mutated before any read")


# Finding D: _CHECK_CALL_RE counts any bash subject containing "check"
# (method_metric.py:96-98). "git status check" is bookkeeping, not a check
# of the agent's own claim. The mark must come from a real test/verify
# family call or verify_ok, not from the word.
def test_adv_word_check_in_a_bash_subject_is_not_a_claim_check():
    from delfin.agent.method_metric import score
    out = score(_record(
        "Works, I am sure.",
        [{"name": "bash", "input": {"command": "git status check"}}],
        verify_ok=False))
    assert out["marks"]["claim_checked"] == 0, (
        "'check' inside a git-status subject counted as claim check")
