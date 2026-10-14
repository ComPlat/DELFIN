"""Harness support for the ``generalist`` benchmark category.

Four mechanisms, each pinned here against the real loader, scorer, runner
and CLI parser (no model call):

* ``prompt_en``: a task carries an English prompt beside the German one;
  ``english_variant`` turns it into a task with its own id (``<id>.en``)
  so the two languages never share a history row or a baseline.
* ``select_tasks``: ``--task`` accepts shell-style globs
  (``generalist_*``); a pattern that matched nothing is returned so the
  CLI can say so.
* ``QUALITY: <0-100>`` printed by an acceptance script is the graded
  part of the verdict. For a task with no expected signals it replaces
  the flat 15/20 routing share, and it is kept on the result
  (``verify_quality``) and through replicate aggregation (median).
* ``workspace_for(task_class="generalist")`` is the generalist fixture,
  and that fixture is restored after every attempt.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import benchmark as bm
from delfin.agent import benchmark_runner as br

_ROOT = Path(__file__).resolve().parents[1]


def _task(**kw) -> bm.Task:
    base = dict(id="generalist_x", task_class="generalist", mode="solo",
                prompt="Mach X.", verify="x.py")
    base.update(kw)
    return bm.Task(**base)


# --- prompt_en ------------------------------------------------------------

def test_a_task_file_carries_an_english_prompt(tmp_path):
    p = tmp_path / "t.yaml"
    p.write_text(
        "tasks:\n"
        "  - id: g1\n"
        "    task_class: generalist\n"
        "    prompt: Mach X.\n"
        "    prompt_en: Do X.\n", encoding="utf-8")
    (t,) = bm.load_tasks(p)
    assert t.prompt == "Mach X."
    assert t.prompt_en == "Do X."


def test_the_english_variant_has_its_own_id_and_prompt():
    t = _task(prompt_en="Do X.")
    en = bm.english_variant(t)
    assert en is not None
    assert en.id == "generalist_x.en"
    assert en.prompt == "Do X."
    # Everything else is the same task: same fixture, same verdict.
    assert (en.setup, en.verify, en.task_class, en.max_tool_calls) == (
        t.setup, t.verify, t.task_class, t.max_tool_calls)


def test_a_task_without_an_english_prompt_has_no_variant():
    assert bm.english_variant(_task()) is None


# --- select_tasks ---------------------------------------------------------

def test_a_glob_selects_a_family_and_an_exact_id_still_works():
    tasks = [_task(id="generalist_a"), _task(id="generalist_b"),
             _task(id="gen_other")]
    got, unmatched = bm.select_tasks(tasks, ["generalist_*", "gen_other"])
    assert [t.id for t in got] == ["generalist_a", "generalist_b",
                                   "gen_other"]
    assert unmatched == []


def test_a_pattern_that_matched_nothing_is_reported():
    got, unmatched = bm.select_tasks([_task(id="generalist_a")],
                                     ["generalist_a", "nosuch_*", "typo"])
    assert [t.id for t in got] == ["generalist_a"]
    assert unmatched == ["nosuch_*", "typo"]


def test_no_pattern_selects_everything():
    tasks = [_task(id="a"), _task(id="b")]
    assert bm.select_tasks(tasks, []) == (tasks, [])


def test_comma_separated_patterns_are_split():
    tasks = [_task(id="generalist_a"), _task(id="x")]
    got, _ = bm.select_tasks(tasks, ["generalist_*,x"])
    assert [t.id for t in got] == ["generalist_a", "x"]


def test_the_cli_accepts_a_language_flag():
    from delfin.agent import cli
    parser = cli.build_parser()
    ns = parser.parse_args(["bench", "run", "--model", "m",
                            "--task", "generalist_*", "--lang", "en"])
    assert ns.lang == "en"
    ns = parser.parse_args(["bench", "run", "--model", "m"])
    assert ns.lang == "de"


# --- QUALITY from the acceptance script -----------------------------------

def _traj(ok, out):
    tr = bm.Trajectory(text="done", duration_s=1.0)
    tr.verify_ok = ok
    tr.verify_output = out
    return tr


def test_the_acceptance_quality_is_kept_and_grades_the_routing_share():
    r = bm.score_outcome(_task(), _traj(True, "PASS a\nQUALITY: 40\nOK"))
    assert r.verify_quality == 40
    assert r.components["routing_pts"] == 8          # 20 * 40 / 100
    assert r.success is True


def test_a_failed_acceptance_keeps_its_partial_quality():
    r = bm.score_outcome(_task(), _traj(False, "FAIL a\nQUALITY: 75\nFAIL"))
    assert r.success is False
    assert r.verify_quality == 75
    assert r.components["routing_pts"] == 15


def test_without_a_quality_line_nothing_changes():
    r = bm.score_outcome(_task(), _traj(True, "OK"))
    assert r.verify_quality is None
    assert r.components["routing_pts"] == 15


def test_a_task_with_signals_keeps_the_signal_share():
    t = _task(expected_signals=(bm.Signal(pattern="done", against="text"),))
    r = bm.score_outcome(t, _traj(True, "QUALITY: 0\nOK"))
    assert r.verify_quality == 0
    assert r.components["routing_pts"] == 15


@pytest.mark.parametrize("line", ["QUALITY: 101", "QUALITY: -3",
                                  "QUALITY: abc", "quality 50"])
def test_a_malformed_quality_line_is_ignored(line):
    r = bm.score_outcome(_task(), _traj(True, line))
    assert r.verify_quality is None


def test_replicates_keep_the_median_acceptance_quality():
    rs = [bm.score_outcome(_task(), _traj(ok, f"QUALITY: {q}"))
          for ok, q in ((True, 100), (False, 20), (True, 80))]
    agg = bm.aggregate_replicates(rs)
    assert agg.verify_quality == 80
    assert agg.per_run_verify_quality == [100, 20, 80]


# --- the fixture workspace ------------------------------------------------

def test_the_generalist_class_writes_into_its_own_fixture():
    ws = br.workspace_for(_ROOT, task_class="generalist")
    assert ws == _ROOT / "tests" / "fixtures" / "generalist_workspace"


def test_the_generalist_fixture_is_restored_after_every_attempt():
    rel = Path("tests") / "fixtures" / "generalist_workspace"
    assert rel in br._BEHAVIOR_WS_RELS
