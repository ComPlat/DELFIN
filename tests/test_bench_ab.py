"""A before/after comparison must be able to say "I cannot tell".

Input: fake runners with fixed outcomes, and a real git repository of two
commits that differ by one benchmark yaml. Output: the usual pass or fail.

SEMANTICS. Each test here fails when the check it names is deleted from
``bench_ab``. The four numbers a single-run probe cannot produce are the ones
under test:

* a clear improvement (0/3 -> 3/3) must be called improved;
* a clear regression (3/3 -> 0/3) must be called regressed;
* the SAME outcome on both sides (3/3 -> 3/3) must say
  ``no measurable difference`` -- not "improved", not "0.0 improvement";
* 3/3 -> 2/3 stays unresolved at N=3, which is the case a rate-only report
  gets wrong and an interval gets right.

WHY THE INTERVAL IS CROSS-CHECKED AGAINST A PUBLISHED EXAMPLE. A statistical
routine tested only against the numbers it produces is tested against itself.
Newcombe (1998, Stat Med 17:873) tabulates his method 10 interval for
9/10 vs 3/10 as 0.1705 to 0.8090; that number is independent of this
repository and of me, so it is the oracle.

NO MODEL IS CALLED. ``runner`` is a fake in every orchestration test. The one
test that touches the real machinery drives the real ``safe_git`` and the real
worktree factory on a throwaway repository, because a fake worktree would prove
only that the fake works.

HOST. Nothing here skips on what this machine has. git is a hard dependency of
the module under test and of this repository's own tooling, so its absence is
an error the test states, not a condition it steps over (see
tests/test_a_test_is_bound_to_delfin_not_to_a_machine.py).
"""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import pytest


class _LazyModule:
    """``delfin.agent.bench_ab``, imported on FIRST ATTRIBUTE ACCESS.

    Constructing this proxy imports nothing, so the module-level ``AB`` below
    does not import ``bench_ab``. That matters because ``bench_ab`` does not
    exist on ``origin/main``: a module-level import would turn every test in
    this file -- including the control that must pass on both sides -- into a
    collection error on the base, and a file that is wholly red there proves
    only that a file was added. With the proxy, the control passes on the base
    and on this branch, and each test that names new behaviour fails on the base
    with ``No module named 'delfin.agent.bench_ab'`` -- the reason the test
    names, rather than a side effect of collection.
    """

    def __getattr__(self, name: str):
        from delfin.agent import bench_ab as _module

        return getattr(_module, name)


AB = _LazyModule()

from delfin.agent.benchmark import Task  # noqa: E402


# ---------------------------------------------------------------------------
# Fixtures and fakes
# ---------------------------------------------------------------------------


def _task(tid: str = "t-alpha", *, behavior: str = "verify",
          task_class: str = "generalist", prompt: str = "do the thing") -> Task:
    return Task(id=tid, task_class=task_class, mode="solo", prompt=prompt,
                expected_signals=(), behavior=behavior)


class _FakeAccounting:
    """The nightly Accounting protocol, recording what it was asked."""

    def __init__(self, estimate: float = 0.1) -> None:
        self._estimate = estimate
        self.charged = 0.0
        self.estimates = 0

    def estimate(self, task) -> float:
        self.estimates += 1
        return self._estimate

    def charge(self, core_hours: float) -> None:
        self.charged += float(core_hours)


class _FakeWorktree:
    def __init__(self, path: Path) -> None:
        self.path = path
        self.branch = "ab/fake"
        self.base_ref = ""


def _factory():
    """A worktree factory that makes a real directory per side, nothing else.

    Real directories matter: ``run`` builds paths from ``side.path``, and a
    nonexistent one would make the fake runner's failure look like a worktree
    failure.

    THE CONTAINMENT, and why it is written this way. An earlier version of this
    helper took the caller's root and registered ``shutil.rmtree(root)`` at
    exit. Two call sites passed the checkout, so the run created
    ``tests/wt-*`` inside the repository (the suite's leak guard reported it,
    correctly) and the exit hook then removed the whole working tree:
    ``rmtree(".")`` deleted this checkout, its ``.git`` included. It was
    recovered from the harness's undo store, which is luck, not a design.

    So: the root is ALWAYS a fresh ``mkdtemp`` -- no caller can aim it at
    anything -- and ``cleanup`` removes exactly that one path and nothing a
    caller hands it. A fake for measuring an orchestrator has no business
    deleting a tree.

    WHY CLEANUP TAKES A ``target`` IT MAY ONLY REFUSE. It is called with the
    root it made (or with nothing); anything else is there so a test can
    demonstrate the refusal. The old guard also required the root's PARENT to
    be named ``tmp*``, which the suite's own ``tempfile.tempdir`` fixture makes
    false -- so it could never delete, and nothing called it either, and every
    ``_run()`` left its directories behind.
    """
    import tempfile

    root = Path(tempfile.mkdtemp(prefix="ab-fake-wt-"))
    made: list[_FakeWorktree] = []

    def make(repo, commit, *, parent=None, prefix="ab"):
        p = root / f"wt-{len(made)}-{commit[:7]}"
        p.mkdir(parents=True, exist_ok=True)
        info = _FakeWorktree(p)
        info.base_ref = commit
        made.append(info)
        return info

    def cleanup(target: Path | None = None) -> None:
        """Remove THIS helper's root, and refuse every other path."""
        if target is not None and Path(target).resolve() != root.resolve():
            return
        shutil.rmtree(root, ignore_errors=True)

    make.made = made
    make.root = root
    make.cleanup = cleanup
    return make


def _runner(outcomes: dict[str, tuple[list[bool], list[bool]]]):
    """A runner returning fixed per-side outcomes.

    ``outcomes`` maps task id -> (base passes, head passes). The side is known
    from ``spec.side.name``, so one table serves both arms and a test cannot
    accidentally assert a result the table did not specify.
    """
    calls: list[tuple[str, str]] = []

    def run(task, spec):
        tid = str(task.id)
        if tid not in outcomes:
            raise AssertionError(f"the runner was asked for {tid!r}, which this "
                                 f"test did not script")
        base, head = outcomes[tid]
        passes = list(base if spec.side.name == "base" else head)
        calls.append((tid, spec.side.name))
        return AB.SideResult(passes=passes, model=spec.model,
                          task_digest=AB.task_digest(task))

    run.calls = calls
    return run


def _run(tmp_path: Path | None = None, outcomes=None, *, tasks=None, **kw):
    """One comparison with fakes, against a fixed task population.

    ``tasks`` defaults to a single fake task and is always passed through to
    ``run(tasks=...)``: without it ``run`` would load the packaged suite (265
    tasks at the time of writing, changing every wave) and the fake runner
    would be asked for ids the test never scripted.
    """
    tasks = list(tasks) if tasks is not None else [_task()]
    runner = kw.pop("runner", None) or _runner(outcomes or {})
    factory = kw.pop("worktree_factory", None) or _factory()
    # Both refs default to the real checkout's HEAD, which is a commit in
    # *this* repository -- resolve_commit needs no fixture for that, and the
    # worktrees are fakes, so nothing is written into the tree.
    return AB.run(
        base=kw.pop("base", "HEAD"), head=kw.pop("head", "HEAD"),
        model=kw.pop("model", "test-model"),
        repeat=kw.pop("repeat", 3),
        repo=kw.pop("repo", Path(__file__).resolve().parents[1]),
        accounting=kw.pop("accounting", _FakeAccounting()),
        tasks=tasks,
        runner=runner, worktree_factory=factory,
        log=lambda *_a, **_k: None,
        **kw,
    )


# ---------------------------------------------------------------------------
# The interval: the published oracle, then the four cases
# ---------------------------------------------------------------------------


def test_the_interval_reproduces_newcombes_published_example():
    """The oracle is the paper, not this code.

    Newcombe (1998) method 10, 9/10 vs 3/10: 0.1705 to 0.8090. The argument
    order here is (base, head), so the published case is the HEAD arm passing
    9 of 10 against a BASE of 3 of 10 -- quoted with its direction rather than
    silently reversed.
    """
    lo, hi = AB.newcombe_interval(3, 10, 9, 10)
    assert lo == pytest.approx(0.1705, abs=5e-4)
    assert hi == pytest.approx(0.8090, abs=5e-4)


def test_the_interval_is_exactly_antisymmetric_in_the_two_sides():
    """Swapping the arms mirrors the interval about zero, exactly.

    A difference of proportions negates when the sides swap, so its interval
    must too. This is the check that caught a real defect while this module was
    written: with the two Wilson spreads paired to the wrong sides the interval
    for 0/3 -> 3/3 came out [1.000, 1.000] -- total certainty from three
    samples, which is impossible and would have been read as a certain win.
    """
    for a, b in [((0, 3), (3, 3)), ((3, 3), (2, 3)), ((1, 3), (2, 3)),
                 ((3, 10), (9, 10))]:
        lo, hi = AB.newcombe_interval(*a, *b)
        rlo, rhi = AB.newcombe_interval(*b, *a)
        assert lo == pytest.approx(-rhi, abs=1e-9)
        assert hi == pytest.approx(-rlo, abs=1e-9)


def test_three_of_three_against_none_of_three_is_not_certain():
    """A three-sample sweep is evidence, not proof.

    The lower bound has to stay clear of 1.0. An interval of [1.0, 1.0] would
    tell a reader that 3 vs 0 cannot be anything but a total improvement, and
    a reader is entitled to act on that number.
    """
    lo, hi = AB.newcombe_interval(0, 3, 3, 3)
    assert hi == 1.0
    assert lo < 0.5, f"a 3-vs-0 sweep claimed a lower bound of {lo}, which " \
                     f"three samples cannot support"


def test_an_empty_side_yields_the_whole_range_and_so_covers_zero():
    lo, hi = AB.newcombe_interval(0, 0, 3, 3)
    assert (lo, hi) == (-1.0, 1.0)


def test_a_clear_improvement_is_called_improved():
    rep = _run(Path(__file__).parent,
               {"t-alpha": ([False] * 3, [True] * 3)})
    row = rep.tasks[0].as_dict()
    assert (row["pass_base"], row["pass_head"]) == (0, 3)
    assert row["diff"] == 1.0
    assert row["contains_zero"] is False
    assert row["verdict"] == AB.VERDICT_IMPROVED


def test_a_clear_regression_is_called_regressed():
    rep = _run(Path("."), {"t-alpha": ([True] * 3, [False] * 3)})
    row = rep.tasks[0].as_dict()
    assert row["diff"] == -1.0
    assert row["contains_zero"] is False
    assert row["verdict"] == AB.VERDICT_REGRESSED


def test_the_same_outcome_on_both_sides_says_no_measurable_difference():
    """The case that matters most: identical sides must NOT report a change."""
    rep = _run(Path("."), {"t-alpha": ([True] * 3, [True] * 3)})
    row = rep.tasks[0].as_dict()
    assert row["diff"] == 0.0
    assert row["contains_zero"] is True
    assert row["verdict"] == AB.VERDICT_NO_DIFFERENCE
    text = AB.format_report(rep)
    assert AB.NO_DIFFERENCE in text
    # The summary line counts the verdicts; "0 improved" is the point, so the
    # assertion reads the count rather than the presence of the word.
    summary = text.split("task(s):")[-1]
    assert "0 improved" in summary and "0 regressed" in summary, (
        f"a comparison that found nothing reported: {summary.strip()[:200]}")
    # And the task ROW must not carry the word as its verdict.
    rows = [ln for ln in text.splitlines() if ln.strip().startswith("t-alpha")]
    assert rows, "the report has no row for the measured task"
    assert AB.VERDICT_IMPROVED not in rows[0]
    assert AB.NO_DIFFERENCE in rows[0]
    assert rep.n_no_difference == 1 and rep.n_improved == 0


def test_a_one_sample_gain_stays_unresolved_at_three_repeats():
    """3/3 -> 2/3 is not a regression at N=3, and saying so is the point.

    This is the case a rate-only report gets wrong: the head is 33 points
    lower and the honest answer is still "not resolvable". If this test ever
    starts failing with AB.VERDICT_REGRESSED, the interval has been loosened.
    """
    rep = _run(Path("."), {"t-alpha": ([True] * 3, [True, True, False])})
    row = rep.tasks[0].as_dict()
    assert row["diff"] == pytest.approx(-1 / 3)
    assert row["verdict"] == AB.VERDICT_NO_DIFFERENCE


def test_an_unmeasured_side_is_not_counted_as_a_failure():
    """A side that measured nothing has no rate -- and no verdict either."""
    rep = _run(Path("."), {"t-alpha": ([], [True] * 3)})
    row = rep.tasks[0].as_dict()
    assert row["rate_base"] is None
    assert row["n_base"] == 0
    assert row["verdict"] == AB.VERDICT_UNMEASURED
    assert row["verdict"] != AB.VERDICT_REGRESSED
    assert "not measured" in AB.format_report(rep)


# ---------------------------------------------------------------------------
# The refusals: a pair that is not a pair
# ---------------------------------------------------------------------------


def test_a_different_model_on_the_two_sides_is_refused():
    """Two models measured against each other is not an A/B of the code."""

    def runner(task, spec):
        return AB.SideResult(passes=[True],
                          model="model-a" if spec.side.name == "base"
                          else "model-b",
                          task_digest=AB.task_digest(task))

    with pytest.raises(AB.AbRefused) as exc:
        _run(Path("."), {}, runner=runner, repeat=1)
    assert "same model" in str(exc.value)


def test_the_requested_model_must_be_the_one_that_ran():
    """A caller passing --model X and a side resolving Y is a mis-filed run."""

    def runner(task, spec):
        return AB.SideResult(passes=[True], model="something-else",
                          task_digest=AB.task_digest(task))

    with pytest.raises(AB.AbRefused) as exc:
        _run(Path("."), {}, runner=runner, model="the-asked-model", repeat=1)
    assert "not the requested" in str(exc.value)


def test_a_task_whose_text_changed_between_the_refs_is_refused():
    """Same id, different question: the difference is not the change's effect."""
    task = _task()

    def runner(task_arg, spec):
        digest = "aaaa" if spec.side.name == "base" else "bbbb"
        return AB.SideResult(passes=[True], model=spec.model, task_digest=digest)

    with pytest.raises(AB.AbRefused) as exc:
        _run(Path("."), {}, tasks=[task], runner=runner, repeat=1)
    assert "DIFFERENT question" in str(exc.value)


def test_the_same_task_text_on_both_sides_is_not_refused():
    """The control for the test above: equal digests must pass.

    Without this pair, a refusal that fires on everything would look like a
    working guard.
    """
    rep = _run(Path("."), {"t-alpha": ([True], [True])}, repeat=1)
    assert rep.tasks[0].verdict == AB.VERDICT_NO_DIFFERENCE


def test_a_repeats_of_zero_is_refused():
    with pytest.raises(AB.AbRefused) as exc:
        _run(Path("."), {}, repeat=0)
    assert "at least 1" in str(exc.value)


def test_no_model_is_refused():
    """Without a model each side may resolve its own; that is not a pair."""
    with pytest.raises(AB.AbRefused) as exc:
        _run(Path("."), {}, model="")
    assert "--model" in str(exc.value)


def test_an_unresolvable_ref_is_refused_and_not_defaulted():
    """A silent fallback to HEAD measures the wrong baseline."""
    with pytest.raises(AB.AbRefused) as exc:
        _run(Path("."), {}, base="no/such/ref-exists-here")
    assert "does not name a commit" in str(exc.value)


def test_a_tree_that_is_not_a_repository_says_so(tmp_path):
    """A copied tree has no ref to fetch, so "Fetch it" is the wrong advice.

    Measured: the same refusal served a directory that is not a repository at
    all, and both real-run sides were told to fetch a ref from a tree with no
    ``.git``. The control is the sibling case directly below, where a real
    repository genuinely lacks the ref and "Fetch it" IS the right advice: that
    one keeps its wording, so this is a branch on the cause, not a rewording of
    one message.
    """
    not_a_repo = tmp_path / "copied-tree"
    not_a_repo.mkdir()
    with pytest.raises(AB.AbRefused) as exc:
        _run(tmp_path, {}, base="origin/main", repo=not_a_repo)
    text = str(exc.value)
    assert "not a git repository" in text, text
    assert "Fetch it" not in text, text
    # The ref is still named: a refusal that hides what failed cannot be acted on.
    assert "origin/main" in text, text


def test_a_selection_matching_nothing_is_refused():
    """An empty population reports a difference of zero and means nothing."""
    with pytest.raises(AB.AbRefused) as exc:
        _run(Path("."), {}, task_ids=["no-such-task-id-anywhere"])
    assert "matched no task" in str(exc.value)


# ---------------------------------------------------------------------------
# The budget gate: the nightly round's own rule
# ---------------------------------------------------------------------------


def test_an_attempt_over_the_daily_budget_is_refused_and_named():
    accounting = _FakeAccounting(estimate=10.0)
    rep = _run(Path("."), {"t-alpha": ([True], [True])}, repeat=1,
               budget=5.0, accounting=accounting)
    row = rep.tasks[0].as_dict()
    # 10 + 0 > 5 on the very first attempt: nothing ran, and the report says
    # so rather than reporting a rate for work that never happened.
    assert accounting.charged == 0.0
    assert row["verdict"] == AB.VERDICT_UNMEASURED
    assert rep.budget_refused, "the refusal was not recorded"
    assert "exceeds" in rep.budget_refused[0]


def test_the_budget_gate_lets_work_through_below_the_limit():
    """The control for the refusal: under the limit, attempts run.

    Without this, a gate stuck closed would look like a working gate.
    """
    accounting = _FakeAccounting(estimate=0.1)
    rep = _run(Path("."), {"t-alpha": ([True], [True])}, repeat=1,
               budget=5.0, accounting=accounting)
    assert not rep.budget_refused
    assert accounting.charged == pytest.approx(0.2)   # two sides, one attempt
    assert rep.tasks[0].verdict == AB.VERDICT_NO_DIFFERENCE


def test_a_negative_estimate_never_reduces_the_days_spend():
    """A negative estimate would let later work fit for free."""
    accounting = _FakeAccounting(estimate=-3.0)
    rep = _run(Path("."), {"t-alpha": ([True], [True])}, repeat=1,
               budget=5.0, accounting=accounting)
    assert accounting.charged == 0.0
    assert any("negative" in r for r in rep.budget_refused)


# ---------------------------------------------------------------------------
# Selection: by tag, and by the files a branch changed
# ---------------------------------------------------------------------------


def test_tag_selection_matches_behavior_and_task_class_as_whole_words():
    tasks = [_task("a", behavior="verify", task_class="science_analysis"),
             _task("b", behavior="plan", task_class="office"),
             _task("c", behavior="", task_class="generalist")]
    got, missing = AB.select_by_tags(tasks, ["verify"])
    assert [t.id for t in got] == ["a"] and missing == []
    got, missing = AB.select_by_tags(tasks, ["generalist"])
    assert [t.id for t in got] == ["c"]
    got, missing = AB.select_by_tags(tasks, ["plan", "nosuchtag"])
    assert [t.id for t in got] == ["b"]
    assert missing == ["nosuchtag"], (
        "an unmatched tag must be reported: a typo that silently ran the rest "
        "reports a number for a question nobody asked")


def test_a_tag_is_not_matched_as_a_substring():
    """`plan` must not select a task class that merely starts with it."""
    tasks = [_task("a", behavior="", task_class="planning_notes")]
    got, missing = AB.select_by_tags(tasks, ["plan"])
    assert got == [] and missing == ["plan"]


def test_a_changed_benchmark_yaml_selects_exactly_the_tasks_it_defines():
    tasks = [_task("only-in-that-file"), _task("somewhere-else")]
    meta = {"only-in-that-file": {"file": "tasks_w99_zz9.yaml", "skill": ""},
            "somewhere-else": {"file": "tasks.yaml", "skill": ""}}
    got, unmatched = AB.select_by_changes(
        tasks, ["delfin/agent/pack/benchmark/tasks_w99_zz9.yaml"],
        metadata=meta)
    assert [t.id for t, _rule in got] == ["only-in-that-file"]
    assert "tasks_w99_zz9.yaml" in got[0][1]
    assert unmatched == []


def test_a_changed_skill_selects_the_tasks_that_declare_it():
    tasks = [_task("uses-skill"), _task("does-not")]
    meta = {"uses-skill": {"file": "tasks.yaml", "skill": "data-figures"},
            "does-not": {"file": "tasks.yaml", "skill": ""}}
    got, _ = AB.select_by_changes(
        tasks, ["delfin/agent/pack/skills/data-figures.md"], metadata=meta)
    assert [t.id for t, _r in got] == ["uses-skill"]
    assert got[0][1].startswith("skill:")


def test_a_change_to_the_turn_path_selects_the_whole_suite():
    """Every task goes through the turn path, so every task is evidence."""
    tasks = [_task("a"), _task("b"), _task("c")]
    got, _ = AB.select_by_changes(tasks, ["delfin/agent/api_client.py"])
    assert sorted(t.id for t, _r in got) == ["a", "b", "c"]


def test_a_specialist_change_selects_the_tasks_that_delegate():
    tasks = [_task("delegating", prompt="Delegate to it (subagent tool) now"),
             _task("plain", prompt="read the file")]
    got, _ = AB.select_by_changes(tasks, ["delfin/agent/subagents.py"])
    assert [t.id for t, _r in got] == ["delegating"]


def test_a_changed_file_no_rule_covers_is_named_and_widens_to_the_suite():
    """The conservative reading, stated rather than silent."""
    tasks = [_task("a"), _task("b")]
    got, unmatched = AB.select_by_changes(tasks, ["delfin/agent/brand_new.py"])
    assert unmatched == ["delfin/agent/brand_new.py"]
    assert sorted(t.id for t, _r in got) == ["a", "b"]
    assert all(r.startswith("unmapped") for _t, r in got)


def test_task_metadata_reads_skill_names_the_task_dataclass_drops():
    """`skill:` is a selection key, and _coerce_task deliberately drops it.

    So the selection route reads it from the suite files; this pins that the
    two do not drift -- if the pack ever renames the key, the routing here is
    silently empty and this fails.
    """
    meta = AB.task_metadata()
    assert meta, "no packaged task carried readable metadata"
    with_skill = {k: v for k, v in meta.items() if v.get("skill")}
    assert with_skill, ("not one packaged task declared a skill, so the "
                        "skill route in this module selects nothing")
    sample_id = sorted(with_skill)[0]
    assert with_skill[sample_id]["file"].startswith("tasks")


# ---------------------------------------------------------------------------
# The report says how it was built, before it says any number
# ---------------------------------------------------------------------------


def test_the_report_states_the_interval_method_and_why():
    rep = _run(Path("."), {"t-alpha": ([True], [True])}, repeat=1)
    text = AB.format_report(rep)
    assert AB.INTERVAL_METHOD in text
    assert AB.INTERVAL_WHY in text
    assert text.index(AB.INTERVAL_METHOD) < text.index("t-alpha"), (
        "the method was printed after the numbers it explains")


def test_a_same_commit_comparison_says_it_measures_noise_alone():
    """origin/main against itself cannot show a change, and must say so."""
    rep = _run(Path("."), {"t-alpha": ([True], [True])}, repeat=1,
               base="HEAD", head="HEAD")
    assert rep.base.commit == rep.head.commit
    text = AB.format_report(rep)
    assert "SAME COMMIT" in text
    assert "false positive" in text


def test_the_report_names_the_rule_that_selected_each_task():
    rep = _run(Path("."), {"t-alpha": ([True], [True])}, repeat=1,
               tags=["verify"])
    assert "selected by:" in AB.format_report(rep)


# ---------------------------------------------------------------------------
# Through the product entry point: cli.main
# ---------------------------------------------------------------------------


def _cli(argv, **patch):
    """Run ``cli.main`` with ``bench_ab.run`` replaced by a scripted fake."""
    from delfin.agent import cli as _cli_mod
    from delfin.agent import bench_ab as _ab_mod

    class _Patch:
        def __init__(self):
            self.seen = None
            self._orig = _ab_mod.run

        def __enter__(self):
            def fake_run(**kw):
                self.seen = kw
                return AB.ABReport(
                    base=AB.Side("base", kw["base"], "a" * 40),
                    head=AB.Side("head", kw["head"], "b" * 40),
                    model=kw["model"], repeat=kw["repeat"],
                    tasks=[AB.TaskAB("t-alpha",
                                  AB.SideResult(passes=[True, True, True]),
                                  AB.SideResult(passes=[True, True, True]))])
            _ab_mod.run = fake_run
            return self

        def __exit__(self, *exc):
            _ab_mod.run = self._orig
            return False

    with _Patch() as p:
        rc = _cli_mod.main(argv)
    return rc, p.seen


def test_bench_ab_reaches_the_comparison_through_cli_main(capsys):
    """A test that calls run() directly only proves run() exists.

    This drives the registered verb with the argv a person would type, so it
    fails if the subparser, the dispatch branch or the option names break --
    which is how a verb gets added and stays unreachable.
    """
    rc, seen = _cli(["bench", "ab", "--model", "m1", "--base", "origin/main",
                     "--head", "HEAD", "--repeat", "3",
                     "--task", "t-alpha", "--tags", "verify"])
    out = capsys.readouterr().out
    assert rc == 0
    assert seen is not None, "cli.main never reached bench_ab.run"
    assert seen["model"] == "m1"
    assert seen["repeat"] == 3
    assert seen["task_ids"] == ["t-alpha"]
    assert seen["tags"] == ["verify"]
    assert "bench ab" in out
    assert AB.INTERVAL_METHOD in out


def test_bench_ab_comma_separated_tasks_and_tags_reach_the_module(capsys):
    rc, seen = _cli(["bench", "ab", "--model", "m",
                     "--task", "t-a,t-b", "--tags", "verify,plan"])
    assert rc == 0
    assert seen["task_ids"] == ["t-a", "t-b"]
    assert seen["tags"] == ["verify", "plan"]


def test_a_comparison_that_measured_nothing_exits_one(capsys):
    """Exit 0 means a comparison took place. An outage is not a comparison.

    Measured 2026-11-25 through this same entry point with no provider key in
    the environment: every attempt died at engine init, the table read
    "0/2 vs 0/2, no measurable difference", and the exit code was 0 -- so a
    wave reading only the exit code would take an endpoint outage for a result
    about the change. What would make this fail: returning 0 when every task's
    verdict is "not measured", or returning 1 when a task DID measure (the
    control below).
    """
    from delfin.agent import cli as _cli_mod
    from delfin.agent import bench_ab as _ab_mod

    orig = _ab_mod.run

    def nothing_measured(**_kw):
        return AB.ABReport(
            base=AB.Side("base", "origin/main", "a" * 40),
            head=AB.Side("head", "HEAD", "b" * 40),
            model="m", repeat=2,
            tasks=[AB.TaskAB(
                "t-alpha",
                AB.SideResult(passes=[], unmeasured=2,
                              errors=["engine init failed: no key"]),
                AB.SideResult(passes=[], unmeasured=2,
                              errors=["engine init failed: no key"]))])

    _ab_mod.run = nothing_measured
    try:
        rc = _cli_mod.main(["bench", "ab", "--model", "m"])
    finally:
        _ab_mod.run = orig
    err = capsys.readouterr().err
    assert rc == 1, (
        "a comparison in which nothing was measured must not exit 0: "
        f"exit {rc} reads as a result about the change")
    assert "nothing was measured" in err, (
        f"the reason must say it, got: {err[-300:]!r}")


def test_a_refusal_exits_two_and_says_why(capsys):
    """A refusal is not a success with an empty table."""

    from delfin.agent import cli as _cli_mod
    from delfin.agent import bench_ab as _ab_mod

    orig = _ab_mod.run

    def refuse(**_kw):
        raise AB.AbRefused("the two sides did not use the same model")

    _ab_mod.run = refuse
    try:
        rc = _cli_mod.main(["bench", "ab", "--model", "m"])
    finally:
        _ab_mod.run = orig
    err = capsys.readouterr().err
    assert rc == 2
    assert "refused" in err and "same model" in err
    assert "Traceback" not in err, (
        "a refusal is a decision, not a crash; a traceback buries the reason")


def test_the_verb_refuses_to_start_without_a_model():
    """--model is required, and argparse says so before any work begins."""
    with pytest.raises(SystemExit) as exc:
        _cli(["bench", "ab", "--base", "origin/main"])
    assert exc.value.code == 2


# ---------------------------------------------------------------------------
# The real machinery, on a throwaway repository
# ---------------------------------------------------------------------------


def _git(repo: Path, *args: str) -> str:
    r = subprocess.run(["git", *args], cwd=str(repo), capture_output=True,
                       text=True)
    assert r.returncode == 0, f"git {' '.join(args)} failed: {r.stderr}"
    return r.stdout


def _repo_with_two_suites(root: Path) -> tuple[Path, str, str]:
    """A real repository: two commits, differing in one benchmark yaml.

    Returns (repo, base_commit, head_commit). The suite files live where
    ``benchmark.load_tasks`` looks, so ``--changes`` can be exercised against
    them for real rather than against a fabricated file list.
    """
    repo = root / "repo"
    pack = repo / "delfin" / "agent" / "pack" / "benchmark"
    pack.mkdir(parents=True)
    _git(repo, "init", "--quiet", "--initial-branch=main")
    _git(repo, "config", "user.email", "t@example.invalid")
    _git(repo, "config", "user.name", "t")
    _git(repo, "config", "commit.gpgsign", "false")
    (repo / "README.md").write_text("fixture\n")
    (pack / "tasks.yaml").write_text(
        "version: 1\ntasks:\n"
        "  - id: kept\n    task_class: generalist\n    mode: solo\n"
        "    prompt: do it\n")
    _git(repo, "add", ".")
    _git(repo, "commit", "--quiet", "-m", "base")
    base = _git(repo, "rev-parse", "HEAD").strip()
    (pack / "tasks.yaml").write_text(
        "version: 1\ntasks:\n"
        "  - id: kept\n    task_class: generalist\n    mode: solo\n"
        "    prompt: do it\n"
        "  - id: added\n    task_class: generalist\n    mode: solo\n"
        "    prompt: do the other thing\n")
    _git(repo, "add", ".")
    _git(repo, "commit", "--quiet", "-m", "head")
    head = _git(repo, "rev-parse", "HEAD").strip()
    return repo, base, head


def test_git_is_present_for_the_real_machinery_this_module_drives():
    """Stated, not skipped over.

    git is what resolves the two refs and makes the two checkouts. A test that
    skipped here would be a test that runs nowhere; this one names the
    dependency instead, the way the host-binding ratchet asks.
    """
    assert shutil.which("git"), (
        "git is required: bench_ab resolves refs with it and builds the two "
        "sides as worktrees")


def test_changed_files_is_the_merge_base_comparison_not_the_branch_tip(
        tmp_path):
    """``...`` and not ``..``: a branch behind main must not report main's moves.

    Measured need: with two dots, every file main advanced is listed as if the
    branch had touched it, and ``--changes`` then selects the wrong tasks.
    """
    repo, base, head = _repo_with_two_suites(tmp_path)
    got = AB.changed_files(repo, base, head)
    assert got == ["delfin/agent/pack/benchmark/tasks.yaml"]


def test_resolve_commit_pins_the_name_to_a_commit(tmp_path):
    repo, base, head = _repo_with_two_suites(tmp_path)
    assert AB.resolve_commit(repo, "main") == head
    assert AB.resolve_commit(repo, base) == base


def test_the_real_worktree_factory_makes_two_different_trees(tmp_path):
    """The default factory, not a fake: two sides must be two checkouts.

    A fake worktree proves the fake works. This drives
    ``worktree.enter_worktree`` at the two commits and asserts each tree
    carries its own commit's content -- which is the property the whole
    comparison rests on.
    """
    from delfin.agent import worktree as W

    repo, base, head = _repo_with_two_suites(tmp_path)
    infos = []
    try:
        for commit in (base, head):
            infos.append(AB.make_worktree(repo, commit,
                                          parent=tmp_path / "wts"))
        assert infos[0].path != infos[1].path
        base_yaml = (infos[0].path / "delfin" / "agent" / "pack" /
                     "benchmark" / "tasks.yaml").read_text()
        head_yaml = (infos[1].path / "delfin" / "agent" / "pack" /
                     "benchmark" / "tasks.yaml").read_text()
        assert "added" not in base_yaml, "the base tree carries head's content"
        assert "added" in head_yaml
        assert _git(infos[0].path, "rev-parse", "HEAD").strip() == base
        assert _git(infos[1].path, "rev-parse", "HEAD").strip() == head
    finally:
        for info in infos:
            W.exit_worktree(info, keep_if_changed=False)


def test_the_selection_follows_a_real_branch_change(tmp_path):
    """End to end: real git, real suite files, the rule that fired.

    The base and head suites differ by one ADDED task, so a selection by the
    changed file must pick the file's tasks and the report must name the file.
    """
    repo, base, head = _repo_with_two_suites(tmp_path)
    files = AB.changed_files(repo, base, head)
    assert files == ["delfin/agent/pack/benchmark/tasks.yaml"]
    # The selection reads the pack of the repository BEING MEASURED, so the
    # fixture's own suite file is what decides. Reading the installed pack
    # would describe a different tree than the one under test.
    got, unmatched = AB.select_by_changes(
        [_task("kept"), _task("added")], files,
        pack_dir=repo / "delfin" / "agent" / "pack" / "benchmark")
    assert unmatched == []
    assert sorted(t.id for t, _r in got) == ["added", "kept"], (
        "a changed task-pack file must select the tasks that file defines")
    assert all("tasks.yaml" in r for _t, r in got)


# ---------------------------------------------------------------------------
# Nothing was written where it does not belong
# ---------------------------------------------------------------------------


def test_the_comparison_leaves_no_worktree_behind(tmp_path):
    """The sides are temporary; a kept worktree would poison the next run."""
    from delfin.agent import worktree as W

    repo, base, head = _repo_with_two_suites(tmp_path)
    made: list = []

    def factory(repo_arg, commit, *, parent=None, prefix="ab"):
        info = W.enter_worktree(repo_arg, branch_prefix=prefix,
                                parent=tmp_path / "wts", base_ref=commit)
        made.append(info)
        return info

    _run(tmp_path, {"t-alpha": ([True], [True])}, repeat=1, repo=repo,
         base=base, head=head, worktree_factory=factory)
    assert len(made) == 2
    for info in made:
        assert not info.path.exists(), (
            f"worktree {info.path} outlived the comparison")


def test_a_runner_that_raises_does_not_take_the_comparison_down(tmp_path):
    """One side failing is recorded, not fatal: the other side still reports.

    A runner crash is a measurement that did not happen, which is a fact about
    the attempt -- not a reason to discard the arm that did run.
    """
    def runner(task, spec):
        if spec.side.name == "head":
            raise RuntimeError("endpoint said no")
        return AB.SideResult(passes=[True], model=spec.model,
                          task_digest=AB.task_digest(task))

    rep = _run(tmp_path, {}, runner=runner, repeat=1)
    row = rep.tasks[0].as_dict()
    assert row["n_base"] == 1 and row["n_head"] == 0
    assert row["verdict"] == AB.VERDICT_UNMEASURED
    assert any("endpoint said no" in e for e in row["errors"])


def test_the_json_payload_carries_the_interval_and_its_method():
    rep = _run(Path("."), {"t-alpha": ([True], [False])}, repeat=1)
    payload = json.loads(json.dumps(rep.as_dict(), default=str))
    assert payload["interval_method"] == AB.INTERVAL_METHOD
    row = payload["tasks"][0]
    for key in ("n_base", "n_head", "rate_base", "rate_head", "diff",
                "interval_low", "interval_high", "contains_zero", "verdict"):
        assert key in row, f"the machine-readable report lost {key}"
