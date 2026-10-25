"""SP1 review round 1 (nacht-s15 findings F1-F3), red controls.

Each case reproduces a finding on the code BEFORE the fix:
- F1: run_pilot scores a zero-citation answer as citations_ok.
- F2: a nesting call through child_runner reaches the real runner.
- F3: refusal paths return a 4-attr SimpleNamespace, not SubagentResult.

Fakes only: the real runner is replaced at module level (the same
binding SA1's own docstring says tests monkeypatch, explore_pilot.py:43-45).
"""
import types

from delfin.agent import explore_pilot as ep


# --------------------------------------------------------------------------
# F1 — zero citations must not pass the bar
# --------------------------------------------------------------------------

def test_f1_zero_citation_text_is_not_fully_cited():
    rep = ep.verify_citations("I read some files and here is my answer.")
    row = ep.TaskRow(task_id="t")
    row.arms["plain"] = ep.ArmResult(
        arm="plain", success=True, tool_calls=5, duration_s=10.0,
        final_text="prose", citations_ok=True)
    # run_pilot itself must set citations_ok False for prose-only text.
    task = ep.PilotTask(task_id="t", prompt="p")
    seen = {}

    def fake_runner(task, arm="plain"):
        seen[arm] = True
        return {"success": True, "final_text": "no citations at all",
                "tool_calls": 3, "duration_s": 2.0}
    rows = ep.run_pilot([task], runner=fake_runner)
    assert rows[0].arms["plain"].citations_ok is False


def test_f1_uncited_win_is_not_a_win():
    """pilot_verdict: an explore arm that wins on speed with ZERO
    citations while plain cited correctly must not count as a win at
    equal correctness."""
    row = ep.TaskRow(task_id="t")
    row.arms["plain"] = ep.ArmResult(
        arm="plain", success=True, tool_calls=5, duration_s=10.0,
        citations_ok=True)
    row.arms["explore"] = ep.ArmResult(
        arm="explore", success=True, tool_calls=1, duration_s=1.0,
        citations_ok=False)
    v = ep.pilot_verdict([row])
    assert v.equal_correctness is False


def test_f1_cited_answer_still_passes():
    """The positive case: a real path:line claim against this tree."""
    rep = ep.verify_citations(
        "see delfin/agent/explore_pilot.py:36 for the type",
        root=".")
    assert not rep.bad
    assert rep.claims == 1


# --------------------------------------------------------------------------
# F2 — nesting must be refused on the real path
# --------------------------------------------------------------------------

def test_f2_child_runner_refuses_nesting_before_the_real_runner():
    calls = []

    def fake_real(**kw):
        calls.append(kw)
        return types.SimpleNamespace(final_text="ran", error="")

    saved = ep._real_run_subagent
    ep._real_run_subagent = fake_real
    try:
        parent = ep.pilot_slot()
        child = ep.child_runner(parent_perms=None, slot=parent)
        result = child(subagent_type="explore", description="d",
                       prompt="try to nest")
        # The child slot is depth 1: the guard must refuse BEFORE the
        # real runner is reached.
        assert calls == [], "nesting call reached the real runner"
        assert getattr(result, "error", ""), "refusal must carry a reason"
    finally:
        ep._real_run_subagent = saved


def test_f2_parent_slot_depth0_still_allowed():
    calls = []

    def fake_real(**kw):
        calls.append(kw)
        return types.SimpleNamespace(final_text="ran", error="")

    saved = ep._real_run_subagent
    ep._real_run_subagent = fake_real
    try:
        runner = ep.child_runner(parent_perms=None, slot=ep.pilot_slot(),
                                 runner=lambda **kw: None)
        # runner= given: the caller owns the child path entirely.
        assert callable(runner)
        top = ep.pilot_slot()
        r = ep.run_pilot_subagent(
            subagent_type="explore", description="d", prompt="ok",
            parent_client=None, parent_perms=None, _pilot_slot=top)
        assert calls, "depth-0 explore delegation must be allowed"
        assert not getattr(r, "error", "")
    finally:
        ep._real_run_subagent = saved


# --------------------------------------------------------------------------
# F3 — refusal shape equals the real SubagentResult
# --------------------------------------------------------------------------

def test_f3_refusal_returns_real_subagentresult():
    from delfin.agent.subagents import SubagentResult
    r = ep.run_pilot_subagent(
        subagent_type="plan", description="d", prompt="p",
        parent_client=None, parent_perms=None)
    assert isinstance(r, SubagentResult)
    for attr in ("final_text", "error", "tool_calls", "input_tokens",
                 "output_tokens", "elapsed_s", "sa_id"):
        assert hasattr(r, attr), f"SubagentResult needs {attr}"
