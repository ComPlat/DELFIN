"""SP1 adversarial tests (reviewer nacht-s15) — attack the pilot harness.

No model call anywhere: every runner is a fake. The live measurement is
the operator's run, never ours.
"""
from types import SimpleNamespace

from delfin.agent import explore_pilot as ep


def _fake_runner(text_by_arm=None, **over):
    calls = []

    def runner(task, arm="plain"):
        calls.append((task.task_id, arm))
        out = {"success": True, "input_tokens": 10, "output_tokens": 20,
               "tool_calls": 3, "duration_s": 1.0,
               "final_text": (text_by_arm or {}).get(arm, "")}
        out.update(over)
        return out

    runner.calls = calls
    return runner


# --- contract pins (green: what the adapter built in this package must meet)


def test_run_pilot_runs_both_arms_per_task(tmp_path):
    runner = _fake_runner({"plain": "see delfin/agent/cli.py:332",
                           "explore": "see delfin/agent/cli.py:520"})
    rows = ep.run_pilot([ep.PilotTask(task_id="t1", prompt="p")],
                        runner=runner, root=tmp_path)
    assert runner.calls == [("t1", "plain"), ("t1", "explore")]
    row = rows[0]
    assert set(row.arms) == {"plain", "explore"}
    for arm in ("plain", "explore"):
        r = row.arms[arm]
        assert r.success is True
        assert r.input_tokens == 10 and r.output_tokens == 20
        assert r.tool_calls == 3 and r.duration_s == 1.0
        assert r.final_text.startswith("see delfin/agent/cli.py:")
        assert r.citation_report.claims == 1


def test_runner_returning_none_is_tolerated(tmp_path):
    def runner(task, arm="plain"):
        return None
    rows = ep.run_pilot([ep.PilotTask(task_id="t", prompt="p")],
                        runner=runner, root=tmp_path)
    assert rows[0].arms["plain"].success is False
    assert rows[0].arms["plain"].final_text == ""


def test_pilot_verdict_mixed_evidence_is_none():
    def mk(tc, dc, ok=True):
        r = ep.ArmResult(arm="x", success=True, tool_calls=tc,
                         duration_s=dc, citations_ok=ok)
        return r
    rows = [ep.TaskRow(task_id="a", arms={
        "plain": mk(5, 5.0), "explore": mk(3, 4.0)}),
        ep.TaskRow(task_id="b", arms={
        "plain": mk(5, 5.0), "explore": mk(7, 6.0)})]
    v = ep.pilot_verdict(rows)
    assert v.explore_wins is None
    assert "1" in v.detail and "losses" in v.detail


# --- attacks on the citation bar


def test_zero_citation_text_is_not_a_fully_cited_answer(tmp_path):
    # The bar: answers MUST carry path:line citations. A text with ZERO
    # claims must not score citations_ok=True — that treats "never cited"
    # as "all cited" and lets a lazy arm through the same bar untouched.
    report = ep.verify_citations("I read some files and found the answer.",
                                 root=tmp_path)
    assert report.claims == 0
    assert not report.bad
    # citations_ok in run_pilot: an uncited answer must be marked NOT ok.
    runner = _fake_runner({"plain": "I read some files and found the answer."})
    rows = ep.run_pilot([ep.PilotTask(task_id="t", prompt="p")],
                        runner=runner, root=tmp_path)
    assert rows[0].arms["plain"].citations_ok is False


def test_prose_abbreviation_is_not_a_path_claim(tmp_path):
    # "e.g.:3" in ordinary prose matches the citation regex (file "e.g",
    # line 3) and is reported as a missing file — a false citation
    # failure on perfectly ordinary English.
    report = ep.verify_citations("Use the guard, e.g.:3, or the slot.",
                                 root=tmp_path)
    assert report.claims == 0
    assert report.bad == []


def test_real_citation_against_real_file(tmp_path):
    f = tmp_path / "mod.py"
    f.write_text("x = 1\n# the guard lives here\n", encoding="utf-8")
    good = ep.verify_citations("mod.py:2", root=tmp_path)
    assert good.claims == 1 and not good.bad
    bad_line = ep.verify_citations("mod.py:99", root=tmp_path)
    assert len(bad_line.bad) == 1 and "out of range" in bad_line.bad[0]
    bad_file = ep.verify_citations("nope.py:1", root=tmp_path)
    assert len(bad_file.bad) == 1 and "missing" in bad_file.bad[0]


# --- attacks on the guard


def test_check_delegate_refuses_everything_but_explore_depth0_plain():
    for typ in ("plan", "general-purpose", "verifier", "Explore", ""):
        ok, why = ep.check_delegate(typ)
        assert not ok and "explore" in why
    ok, why = ep.check_delegate("explore", depth=1)
    assert not ok and "nesting" in why
    for iso in ("worktree", ""):
        pass
    ok, why = ep.check_delegate("explore", isolation="worktree")
    assert not ok and "isolation" in why


def test_check_delegate_refuses_a_second_running_delegate():
    slot = ep.pilot_slot()
    ep.claim_slot(slot)
    ok, why = ep.check_delegate("explore", slot=slot)
    assert not ok and "one" in why
    ep.reap_slot(slot)
    ok, why = ep.check_delegate("explore", slot=slot)
    assert ok and why == ""


def test_run_pilot_subagent_never_raises_and_never_reaches_real_runner(
        monkeypatch):
    seen = []

    def fake_real(**kw):
        seen.append(kw)
        return SimpleNamespace(subagent_type=kw["subagent_type"],
                               final_text="ok", error="")

    monkeypatch.setattr(ep, "_real_run_subagent", fake_real)
    r = ep.run_pilot_subagent(subagent_type="plan", description="d",
                              prompt="p", parent_client=None,
                              parent_perms=None)
    assert not seen and r.error and "plan" in r.error
    r2 = ep.run_pilot_subagent(subagent_type="explore", description="d",
                               prompt="p", parent_client=None,
                               parent_perms=None)
    assert len(seen) == 1 and r2.final_text == "ok"
    monkeypatch.setattr(ep, "_real_run_subagent",
                        lambda **kw: (_ for _ in ()).throw(RuntimeError("x")))
    r3 = ep.run_pilot_subagent(subagent_type="explore", description="d",
                               prompt="p", parent_client=None,
                               parent_perms=None)
    assert r3.final_text == "" and "delegate failed" in r3.error


def test_child_runner_refuses_nesting_on_the_real_path(monkeypatch):
    seen = []

    def fake_real(**kw):
        seen.append(kw)
        return SimpleNamespace(final_text="ok", error="")

    monkeypatch.setattr(ep, "_real_run_subagent", fake_real)
    parent = ep.pilot_slot()
    child = ep.child_runner(parent_perms=None, slot=parent)
    r = child(subagent_type="explore", description="d", prompt="p")
    assert not seen and r.error and "nesting" in r.error


def test_refusal_result_carries_the_full_subagentresult_shape(monkeypatch):
    # explore_pilot.py:130-131 claims the refusal result is "the same shape
    # the real runner returns". The real runner returns a SubagentResult
    # (subagents.py:1530-1553) whose consumers read tool_calls,
    # input_tokens, output_tokens, elapsed_s and sa_id. A SimpleNamespace
    # with four attributes makes every such reader crash on the exact
    # path the wrapper exists for.
    from delfin.agent.subagents import SubagentResult

    r = ep.run_pilot_subagent(subagent_type="plan", description="d",
                              prompt="p", parent_client=None,
                              parent_perms=None)
    assert isinstance(r, SubagentResult)
