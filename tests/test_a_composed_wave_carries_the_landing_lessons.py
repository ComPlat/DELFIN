"""A composed wave says, in the briefs themselves, how a landing is decided.

What a landing must not do was learned by measurement and then typed by hand
into one wave's tasks -- outside the product, so the next wave did not have it.
Measured on w42i/w48, the cost of each missing line:

* the bundle was built on a SECOND branch, so every reviewer's checkout had an
  empty diff and four reviews were written against nothing;
* a package whose test was red on the bundle was merged anyway;
* a package that reported "lands only with scratch/<file>.patch" was merged
  without the patch -- the code never reached the base;
* a package no entry point reaches landed as a capability;
* commit messages carried the product name and the ``Co-Authored-By:`` /
  ``Agent: ... via ...`` trailers, which are public text on a merged branch;
* an integrator branch 0 commits ahead of the base was bundled as if it carried
  work, and the reviewer's ``collect`` held only the lander's scratch/, so no
  reviewer ever read the integrators' reports.

These tests drive ``compose.compose()`` -- the function ``delfin-agent swarm
wave compose`` and the ``wave_compose`` MCP tool both call -- on a two-package
spec, and assert on the tasks it produces. ``test_the_arm_and_integrator_briefs_
keep_what_they_had`` is the control: it passes on the base and after, because
those briefs are not what this changes.
"""
from __future__ import annotations

import pytest

from delfin.agent.swarm import compose

# Two families: compose refuses reviewers of the lander's family, so a spec
# needs two real ids from land.MODEL_FAMILY.
_M = {"glm": {"provider": "blablador", "model": "alias-glm-5.3-flash"},
      "ds": {"provider": "kit", "model": "kit.deepseek-v4-flash"}}

_PK = [{"id": "p1", "title": "one", "body": "do one"},
       {"id": "p2", "title": "two", "body": "do two"}]


def _composed(run: str = "w1", **kw) -> dict:
    """One wave composed through the product's own entry, by task id."""
    spec = dict(arms=[_M["glm"], _M["ds"]], integrator=_M["ds"], lander=_M["ds"],
                reviewer=_M["glm"], fixer=_M["ds"], base="origin/main")
    spec.update(kw)
    return {t["id"]: t for t in compose.compose(run, list(_PK), **spec)}


# -- the lander's brief ----------------------------------------------------------

def test_the_lander_builds_on_its_own_branch_and_names_it():
    """Rule A. Fails on a brief that says "one branch" without naming it, or
    names a second one: the reviewers check out wave/<run>/zz-land."""
    brief = _composed()["zz-land"]["brief"]
    assert "wave/w1/zz-land" in brief, brief[:400]
    assert "second branch" in brief.lower(), "rule A is missing: the bundle goes on THIS branch"


def test_the_lander_carries_every_landing_rule_by_name():
    """Each rule has a heading distinctive enough that a brief without it is
    caught: this fails the moment a rule stops being composed in."""
    brief = _composed()["zz-land"]["brief"]
    for heading in ("ONE BRANCH", "NOTHING RED LANDS", "NEEDS A PATCH",
                    "NO ENTRY POINT", "PUBLIC TEXT", "0 COMMITS AHEAD"):
        assert heading in brief, f"{heading!r} missing from the lander's brief"


def test_the_lander_is_told_to_apply_or_decline_a_patch():
    """Rule C, in operational words: apply the patch as its own commit, or say
    why the package is out."""
    brief = _composed()["zz-land"]["brief"].lower()
    assert "scratch/" in brief and ".patch" in brief
    assert "as its own commit" in brief, "the patch must become a commit with its test"
    assert "leave the package out" in brief, "or the package is out, with a reason"


def test_the_lander_may_not_land_what_no_entry_point_reaches():
    """Rule D, asserted on the BUNDLING decision.

    Not on "entry point" / "import is not a caller": both phrases are already
    in DELIVERY_RULES, so a brief with no landing rule at all still contains
    them and that assertion is green on the base -- it measured nothing. The
    phrase asserted here exists only in the landing rule.
    """
    brief = _composed()["zz-land"]["brief"].lower()
    assert "does not land as a capability" in brief, brief[:600]


def test_the_lander_strips_the_agent_trailers_from_public_commits():
    """Rule E names both trailer forms the product writes mechanically, so a
    member that only hears 'no personal data' still leaves them in."""
    brief = _composed()["zz-land"]["brief"]
    assert "Co-Authored-By" in brief, "the co-author trailer is public text"
    assert "Agent:" in brief, "the `Agent: ... via ...` trailer is public text"
    assert "assistant product" in brief.lower()


def test_the_lander_checks_a_branch_actually_has_commits():
    """Rule F carries the command, not just the advice."""
    brief = _composed()["zz-land"]["brief"]
    assert "git log --oneline origin/main.." in brief
    assert "wc -l" in brief


# -- the reviewer's task ---------------------------------------------------------

def test_the_reviewer_collects_every_integrator_of_the_run():
    """The reviewers read the integrators' reports, not only the lander's.
    Fails on collect == [lander], which is what left four reviewers blind."""
    by = _composed()
    want = {"w1/zz-land", "w1/p1-int", "w1/p2-int"}
    for key in compose.REVIEW_FOCUS:
        got = set(by[f"zz-rev-{key}"]["collect"])
        assert got == want, (key, sorted(got))


def test_a_reviewer_of_a_continued_wave_still_sees_every_integrator():
    """``arm_refs`` names arms that already ran, so the integrator reads and
    collects THOSE; its own report still lands in this run. The reviewers must
    collect this run's integrator dirs either way -- the case that went blind
    was a reviewer holding only the lander's scratch/."""
    by = _composed(arm_refs={"p1": ["old/p1-x"]})
    assert set(by["zz-rev-safety"]["collect"]) == {
        "w1/zz-land", "w1/p1-int", "w1/p2-int"}
    assert by["p1-int"]["collect"] == ["old/p1-x"], "the integrator reads the old arms"


def test_the_reviewer_counts_the_bundle_before_anything_else():
    """An empty bundle is the first finding, and the count is how it is seen."""
    brief = _composed()["zz-rev-behaviour"]["brief"]
    assert "git log --oneline origin/main..HEAD | wc -l" in brief
    assert "first line" in brief.lower()


def test_the_reviewer_runs_the_suite_the_way_the_lander_recorded_it():
    """A missing pytest is the wrong interpreter, not a red suite."""
    brief = _composed()["zz-rev-tests"]["brief"].lower()
    assert "wrong interpreter" in brief, brief[:400]
    assert "the lander" in brief and "recorded" in brief


# -- the fixer's brief -----------------------------------------------------------

def test_the_fixer_carries_the_rules_that_bind_its_own_lines():
    """(B) nothing red and (E) public text, plus the branch it works on."""
    brief = _composed()["zz-fix"]["brief"]
    assert "NOTHING RED LANDS" in brief and "PUBLIC TEXT" in brief
    assert "wave/w1/zz-land" in brief, "the fixer works on the lander's branch"


def test_the_fixer_gets_a_subset_and_says_which():
    """The subset is a decision, not an accident: the fixer changes no merge,
    so rule A is not its rule. Adding it is a deliberate edit to FIXER_RULES."""
    brief = _composed()["zz-fix"]["brief"]
    assert "ONE BRANCH" not in brief
    assert set(compose.FIXER_RULES) <= set(compose.LANDING_RULES)


# -- the rules are data ----------------------------------------------------------

def test_the_rules_are_data_a_brief_can_pick_a_subset_from():
    """Six rules, keyed, each naming what its absence cost. Fails if a rule is
    welded into the lander's prose (then the fixer cannot carry it)."""
    assert list(compose.LANDING_RULES) == ["A", "B", "C", "D", "E", "F"]
    assert tuple(compose.LANDER_RULES) == ("A", "B", "C", "D", "E", "F")
    for key, text in compose.LANDING_RULES.items():
        assert len(text) > 80, key
        assert text == text.strip() and not text.endswith(" ")


# -- the product entry point that composes a wave --------------------------------

def test_the_composed_rules_reach_the_file_the_cli_writes(tmp_path, capsys):
    """`delfin-agent swarm wave compose RUN SPEC --out TASKS` is how a wave is
    actually built (cli_swarm._cmd_wave -> mcp_server.h_compose -> compose).
    Asserting here rather than only on compose.compose() means a brief that
    lost its rules fails on the artifact an operator starts from."""
    import json

    from delfin.agent import cli_swarm
    from delfin.agent.cli import build_parser

    spec = {"packages": _PK, "arms": [_M["glm"], _M["ds"]], "integrator": _M["ds"],
            "lander": _M["ds"], "reviewer": _M["glm"], "fixer": _M["ds"],
            "base": "origin/main"}
    spec_file = tmp_path / "spec.json"
    spec_file.write_text(json.dumps(spec), encoding="utf-8")
    out = tmp_path / "tasks.json"
    ns = build_parser().parse_args(
        ["swarm", "wave", "compose", "w9", str(spec_file), "--out", str(out)])
    assert ns.swarm_cmd == "wave" and ns.wave_cmd == "compose"
    assert cli_swarm.run_args(ns) == 0, capsys.readouterr().out

    tasks = {t["id"]: t for t in json.loads(out.read_text(encoding="utf-8"))}
    land = tasks["zz-land"]["brief"]
    assert "wave/w9/zz-land" in land and "NOTHING RED LANDS" in land
    assert "Co-Authored-By" in land and "Agent:" in land
    assert set(tasks["zz-rev-claims"]["collect"]) == {
        "w9/zz-land", "w9/p1-int", "w9/p2-int"}
    assert "NOTHING RED LANDS" in tasks["zz-fix"]["brief"]


def test_every_stage_that_commits_carries_the_public_text_rule():
    """An arm and an integrator commit too, on branches somebody else merges.

    Rule E is about their commit messages, and the product splices both
    trailer forms into every ``git commit`` a model bash call runs
    (api_client._commit_trailer_rewrite) regardless of what its brief says:
    a member that is never told so cannot strip them. Red on the base, whose
    arm and integrator briefs carry no landing rule at all -- the repo-wide
    personal-data scan cannot see this, because the trailer is written at
    commit time, not in the tree.
    """
    by = _composed()
    for tid in ("p1-a1", "p1-a2", "p1-int", "p2-int"):
        assert "Co-Authored-By" in by[tid]["brief"], tid
        assert "PUBLIC TEXT" in by[tid]["brief"], tid


# -- control: what this must NOT change -----------------------------------------

def test_the_arm_and_integrator_briefs_keep_what_they_had():
    """CONTROL -- green on the base and after. The delivery rules and the real
    runs stay where they were, and no landing rule leaks into an arm's brief
    (an arm integrates nothing: it would only dilute its own package text)."""
    by = _composed()
    real = compose.REAL_RUNS.strip()
    for tid in ("p1-a1", "p1-int", "zz-land"):
        assert compose.DELIVERY_RULES.strip() in by[tid]["brief"], tid
    assert real in by["p1-int"]["brief"] and real not in by["p1-a1"]["brief"]
    for tid in ("p1-a1", "p1-a2", "p1-int"):
        assert "NOTHING RED LANDS" not in by[tid]["brief"], tid


def test_the_composed_shape_is_still_the_shape_every_wave_had():
    """CONTROL -- the task list itself is unchanged: same ids, same order, same
    dependencies, same family rule."""
    ids = [t["id"] for t in compose.compose("w", list(_PK),
                                            arms=[_M["glm"], _M["ds"]], integrator=_M["ds"],
                                            lander=_M["ds"], reviewer=_M["glm"], fixer=_M["ds"])]
    assert ids == ["p1-a1", "p1-a2", "p1-int", "p2-a1", "p2-a2", "p2-int", "zz-land",
                   "zz-rev-behaviour", "zz-rev-tests", "zz-rev-safety", "zz-rev-claims",
                   "zz-fix"]
    with pytest.raises(Exception, match="not independent"):
        compose.compose("w", list(_PK), arms=[_M["ds"]], integrator=_M["ds"],
                        lander=_M["ds"], reviewer=_M["ds"], fixer=_M["ds"])
