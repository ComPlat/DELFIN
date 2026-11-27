"""The office benchmark grades the FILE, and the file check has teeth.

The office suite was graded on the answer: a sentence about the work. A
person doing the work is graded on the workbook, and the two are not the
same question -- a confident sentence over an untouched file reads exactly
like a finished job. Package O1 adds twelve everyday tasks whose verdict
comes from ``accept/o1_office_accept.py``, which opens what the attempt
wrote.

TWO DIRECTIONS, and both are necessary:

* **rejects an untouched workspace.** Otherwise every task passes for free
  and the suite measures nothing.
* **accepts a correct one.** Otherwise no task can ever pass, and the suite
  measures the checker.

A third direction is what makes the first two worth anything: for each case,
ONE plausible wrong answer -- the half-done result a real attempt produces,
not a corrupt file -- must be REJECTED, and the message must name the defect
rather than a crash. Each variant is the trap the case was built around:
rows dropped from a split, a sheet per raw value that never merges the
case/space variants, an inner join, only the second copy of a duplicate
marked, one date format understood out of two, an amount left as text, a
merge that stops after the first row, a case-sensitive duplicate test.

The reference solution lives in ``accept/o1_office_reference.py``: it is the
thing that "correct" means here, and it is committed beside the checker so
the two cannot drift.

The suite is run through the real acceptance entry point
(``benchmark_runner.run_acceptance``), not by importing a helper, because
that is what a benchmark run calls -- and it is the thing that has to keep
working.
"""
from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from delfin.agent.benchmark import load_tasks
from delfin.agent.benchmark_runner import acceptance_path, run_acceptance

_ROOT = Path(__file__).resolve().parents[1]
_PACK = _ROOT / "delfin" / "agent" / "pack" / "benchmark"
_SETUP = _PACK / "setup" / "o1_office_fixtures.py"
_REFERENCE = _PACK / "accept" / "o1_office_reference.py"
_ACCEPT = "o1_office_accept.py"

#: Every case, with the task id that grades it in office mode. The solo twin
#: asks the same question under another id; ``test_both_modes_ask_the_same_question``
#: pins that, so this list does not have to name both.
CASES = (
    "blaettern", "zusammenfuehren", "dubletten", "monatssummen",
    "serienbrief", "formular", "pdf_tabelle", "bereinigen", "pivot",
    "csv_deutsch", "regional_rollup", "rechnungs_dubletten",
)

#: case -> a plausible wrong answer this checker must reject.
DEFECTS = {
    "blaettern": "raw",
    "zusammenfuehren": "inner",
    "dubletten": "second-only",
    "monatssummen": "dotted-only",
    "serienbrief": "first-only",
    "formular": "one-only",
    "pdf_tabelle": "text-hours",
    "bereinigen": "raw",
    "pivot": "sums",
    "csv_deutsch": "text",
    "regional_rollup": "north-only",
    "rechnungs_dubletten": "strict",
}


def _env() -> dict:
    import os
    env = dict(os.environ)
    prior = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = os.pathsep.join(
        [str(_ROOT)] + [p for p in prior.split(os.pathsep) if p])
    return env


def _build(workspace: Path) -> None:
    """Build the twelve fixtures the way the runner does: the setup script."""
    proc = subprocess.run([sys.executable, str(_SETUP), str(workspace)],
                          capture_output=True, text=True, timeout=300,
                          env=_env(), cwd=str(workspace))
    assert proc.returncode == 0, proc.stderr[-1500:]


def _solve(workspace: Path, defect: str = "") -> None:
    argv = [sys.executable, str(_REFERENCE), str(workspace)]
    if defect:
        argv += ["--defect", defect]
    proc = subprocess.run(argv, capture_output=True, text=True, timeout=600,
                          env=_env(), cwd=str(_ROOT))
    assert proc.returncode == 0, proc.stderr[-1500:]


@pytest.fixture(scope="module")
def built(tmp_path_factory) -> Path:
    """One untouched fixture workspace, shared: building it takes ~0.5 s."""
    ws = tmp_path_factory.mktemp("o1-untouched")
    _build(ws)
    return ws


# ---------------------------------------------------------------------------
# The suite is registered in both modes
# ---------------------------------------------------------------------------


def test_every_case_has_an_office_task_and_a_solo_twin():
    tasks = {t.id: t for t in load_tasks()}
    missing = []
    for case in CASES:
        office = f"o1_office_{case}"
        solo = f"o1_solo_{case}"
        if office not in tasks:
            missing.append(office)
        if solo not in tasks:
            missing.append(solo)
    assert not missing, f"not registered: {missing}"


def test_both_modes_ask_the_same_question():
    """A solo twin is the same question in another mode, not another task.

    If the prompts drifted, a comparison between the two modes would be
    comparing two questions and could not say which mode did better.
    """
    tasks = {t.id: t for t in load_tasks()}
    drifted = []
    for case in CASES:
        office, solo = tasks[f"o1_office_{case}"], tasks[f"o1_solo_{case}"]
        if office.mode != "office":
            drifted.append(f"{office.id}: mode is {office.mode!r}")
        if solo.mode != "solo":
            drifted.append(f"{solo.id}: mode is {solo.mode!r}")
        if office.prompt != solo.prompt:
            drifted.append(f"{case}: the two prompts differ")
        if office.verify != solo.verify:
            drifted.append(f"{case}: the two graders differ")
        if office.setup != solo.setup:
            drifted.append(f"{case}: the two fixtures differ")
    assert not drifted, "\n  ".join(drifted)


def test_the_office_tasks_run_in_office_mode_where_the_workspace_is():
    """The existing rule, restated for this suite: office mode means the
    folder. A solo twin runs in the same fixture under its own class."""
    from delfin.agent.benchmark_runner import workspace_for

    tasks = {t.id: t for t in load_tasks()}
    for case in CASES:
        office = tasks[f"o1_office_{case}"]
        assert workspace_for(_ROOT, mode=office.mode,
                             task_class=office.task_class) is not None, case


# ---------------------------------------------------------------------------
# The fixtures describe the column that was measured
# ---------------------------------------------------------------------------


def test_the_collaborator_column_has_the_measured_shape(built):
    """The messiness is the question. A tidy fixture measures something else.

    Counted here, not trusted from a comment: the numbers are the ones the
    request failed on (737 rows, 192 distinct values, 49 empty cells, 47
    with a character a sheet name may not carry, 18 over 31 characters, 4
    differing only in case or spacing, 7 naming several people).
    """
    import importlib.util

    spec = importlib.util.spec_from_file_location("o1_setup_check", _SETUP)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    got = mod.self_check()
    assert got["data_rows"] == 737, got
    assert got["raw_values"] == 192, got
    assert got["empty_cells"] == 49, got
    assert got["with_forbidden_char"] == 47, got
    assert got["longer_than_31"] == 18, got
    assert got["case_or_space_variants"] == 4, got
    assert got["naming_several_people"] == 7, got
    # And the shape has a consequence, as an IDENTITY rather than a guess:
    # the four variants merge into four sheets that already exist (-4), each
    # of the seven two-person cells disappears as a sheet of its own and its
    # two people appear (+14 - 7), and the empty cells get one sheet (+1).
    # 192 - 4 + 14 - 7 + 1 = 196. Written out because the first version of
    # this assertion guessed "fewer sheets than values" -- which is the
    # wrong shape for this fixture, and would have hidden a fixture that
    # lost its two-person cells.
    expect_sheets = (got["raw_values"] - got["case_or_space_variants"]
                     + 2 * got["naming_several_people"]
                     - got["naming_several_people"] + 1)
    assert got["sheets"] == expect_sheets, got
    assert got["sheets"] != got["raw_values"], (
        "the merging and the splitting cancelled out, so the fixture no "
        "longer exercises either")


def test_a_sheet_title_never_carries_a_forbidden_character(built):
    """The expected map is what the checker compares against, so it has to
    satisfy the rule the fixture states -- otherwise the task is impossible.
    """
    import importlib.util
    import openpyxl

    spec = importlib.util.spec_from_file_location("o1_setup_titles", _SETUP)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    titles = list(mod.expected_sheets())
    assert len(titles) == len(set(titles)), "two sheets with one title"
    bad = [t for t in titles
           if any(c in t for c in mod.FORBIDDEN_IN_TITLE)
           or len(t) > mod.MAX_TITLE or not t]
    assert not bad, f"expected sheet titles Excel would refuse: {bad[:5]}"
    # And openpyxl agrees they are writable, which is the check that would
    # fail if the rule and Excel's rule ever diverged.
    wb = openpyxl.Workbook()
    wb.remove(wb.active)
    for title in titles[:20]:
        wb.create_sheet(title=title)
    wb.close()


# ---------------------------------------------------------------------------
# Untouched: rejected
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("case", CASES)
def test_an_untouched_workspace_is_rejected(built, case):
    ok, output = run_acceptance(_ACCEPT, built, root=_ROOT,
                               task_id=f"o1_office_{case}")
    assert ok is False, (
        f"{case}: an untouched workspace was accepted ({output[:300]})")


# ---------------------------------------------------------------------------
# Reference: accepted
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def solved(tmp_path_factory) -> Path:
    ws = tmp_path_factory.mktemp("o1-solved")
    _build(ws)
    _solve(ws)
    return ws


@pytest.mark.parametrize("case", CASES)
def test_the_reference_solution_is_accepted(solved, case):
    ok, output = run_acceptance(_ACCEPT, solved, root=_ROOT,
                               task_id=f"o1_office_{case}")
    assert ok is True, f"{case}: {output[:600]}"


def test_the_solo_twin_grades_the_same_file_the_same_way(solved):
    """Same workspace, same verdict under the solo id: the mode is the only
    difference, which is what makes the comparison between them meaningful.
    """
    for case in CASES:
        ok, output = run_acceptance(_ACCEPT, solved, root=_ROOT,
                                   task_id=f"o1_solo_{case}")
        assert ok is True, f"solo {case}: {output[:400]}"


# ---------------------------------------------------------------------------
# Half-done: rejected, for the reason the case is about
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("case", CASES)
def test_a_plausible_wrong_answer_is_rejected(tmp_path, case):
    """The other half of a real check: a result that looks finished.

    Each variant is what a script written without reading HINWEISE.txt
    produces -- an inner join, the second copy only, one date format -- and
    every one of them is the kind of answer that survives a review. If this
    passes for a case, that case's check measures nothing.
    """
    ws = tmp_path / "defect"
    ws.mkdir()
    _build(ws)
    _solve(ws, defect=f"{case}:{DEFECTS[case]}")
    ok, output = run_acceptance(_ACCEPT, ws, root=_ROOT,
                               task_id=f"o1_office_{case}")
    assert ok is False, (
        f"{case}: the half-done answer ({DEFECTS[case]}) was ACCEPTED, so the "
        f"check cannot tell it from a correct one ({output[:400]})")
    # Rejected for the case's own reason, not by a crash: an exception in the
    # checker also yields False, and would leave every case unpassable.
    assert "cannot be graded" not in output, (
        f"{case}: rejected by an exception rather than a finding: {output[:400]}")


def test_a_missing_result_file_is_not_graded_as_a_wrong_answer(built):
    """The distinction the scorer reads: False is a finding about the file;
    None is "not run". A missing script must not look like a failed attempt.
    """
    ok, output = run_acceptance("o1_no_such_accept.py", built, root=_ROOT)
    assert ok is None, f"a missing acceptance script read as {ok!r}"
    assert "not found" in output


def test_the_setup_refuses_to_build_over_its_own_output(built, tmp_path):
    """A workspace that was not cleaned is reported as unmeasured rather than
    silently half-rebuilt -- the rule every setup in this pack follows.
    """
    proc = subprocess.run([sys.executable, str(_SETUP), str(built)],
                          capture_output=True, text=True, timeout=120,
                          env=_env(), cwd=str(built))
    assert proc.returncode != 0
    assert "refusing to overwrite" in proc.stderr


def test_the_setup_says_what_it_wants_when_given_nothing():
    proc = subprocess.run([sys.executable, str(_SETUP)],
                          capture_output=True, text=True, timeout=60,
                          env=_env(), cwd=str(_ROOT))
    assert proc.returncode != 0
    assert "usage" in (proc.stdout + proc.stderr).lower()
