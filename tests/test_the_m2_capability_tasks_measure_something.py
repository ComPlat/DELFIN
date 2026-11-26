"""Each new m2 capability task is graded right from wrong.

Rule 5: an assertion that cannot fail is not a test. For every task in
``tasks_w49_m2.yaml`` this file records TWO answers -- one a correct
agent would write, one a broken capability would produce -- and drives
both through the REAL ``benchmark_runner.run_acceptance``, which launches
``python <script> <workspace>`` exactly as a benchmark run does.

The failing answer is not random text. It is the specific wrong answer
the capability's failure mode yields:

  grep_file      every file containing the name, including the prose note
                 and the unrelated local closure -- what a plain text
                 search that does not read the code reports
  bash           the pooled mean over all 13 rows, 116.00 s, which is no
                 run's mean -- what an agent that cannot run code guesses
  read_file      the number from the first and only sheet it opened, or
                 the spare pump's -- a single read, no comparison
  read_document  150 C, the SUPERSEDED limit in the .txt beside the
                 .docx -- the answer an agent that cannot open the .docx
                 gives, and the reason this task exists
  search_calcs   the methanol folder ranked into the ethanol comparison --
                 a ranking across methods and molecules

A verifier that passed both answers would be measuring nothing, so each
test asserts the pass AND the fail.

Nothing here depends on what this host has installed: the fixtures are
built by the pack's own setup script with the standard library, and a
``.docx`` is written as a zip of XML parts rather than through a document
package.
"""

from __future__ import annotations

import ast
import json
from pathlib import Path

import pytest

from delfin.agent.benchmark import load_tasks
from delfin.agent.benchmark_runner import run_acceptance, run_setup

ROOT = Path(__file__).resolve().parents[1]

# Correct answers, as an agent that used the capability would write them.
# The prose is deliberately short: the verifiers grade FACTS, and a long
# answer must not be able to pass by mentioning everything.
_GOOD = {
    "m2_a1_grep_file_where_is_it_defined": (
        "g1/answer.md",
        "# Where `_settle` lives\n\n"
        "Defined once, in `src/core.py`, as `def _settle(grid, dt)`.\n\n"
        "Called once for real: `step()` in `src/core.py` calls\n"
        "`_settle(grid, 0.01)` inside its loop. That is the only call "
        "site of the helper.\n\n"
        "`src/run.py` imports the name, but the `_settle` defined there "
        "is a local closure (`def _settle(): return grid`), a separate "
        "and unrelated definition -- not a call of the core helper.\n\n"
        "`notes/README.md` mentions it in prose only; it is documentation, "
        "not a definition or a call.\n"),
    "m2_a2_bash_mean_cycle_time": (
        "b2/answer.md",
        "# Mean cycle times\n\n"
        "| run | rows | mean cycle time (s) |\n"
        "|---|---|---|\n"
        "| run_a | 5 | 120.20 |\n"
        "| run_b | 5 | 90.80 |\n"
        "| run_c | 3 | 151.00 |\n\n"
        "run_b is the fastest (90.80 s). The runs have different row "
        "counts, so a pooled mean over all 13 rows would not be any "
        "run's mean.\n"),
    "m2_a3_read_file_two_spec_sheets": (
        "a3/answer.md",
        "# Circulation pumps in service\n\n"
        "| pump | max flow rate | max head |\n|---|---|---|\n"
        "| A (PA-200) | 18.4 L/min | 6.2 m |\n"
        "| B (PB-310) | 11.9 L/min | 9.8 m |\n\n"
        "Pump B has the higher maximum head (9.8 m against 6.2 m). Pump A "
        "has the higher flow rate.\n\n"
        "The folder also holds `pump_spare.txt` for a PX-900: the sheet "
        "says it is a spare and not in service, so it is excluded from "
        "the comparison.\n"),
    "m2_a5_read_document_the_current_note": (
        "d4/answer.md",
        "# Method note MN-14\n\n"
        "The digestion must not exceed **120 C**; above that selenium is "
        "lost. Digest at 95 C for 2 h.\n\n"
        "Results must be reported as **mg/kg dry mass**, not as mg/L.\n\n"
        "`method_note_old.txt` is the superseded 2019 note; its 150 C "
        "limit does not apply.\n"),
    "m2_a6_search_calcs_within_method": (
        "sv1/answer.md",
        "# Energies within each level of theory\n\n"
        "Ranked **within each method** -- a total energy compares only "
        "within one level of theory, never across methods.\n\n"
        "| method | lowest calculation | energy (Eh) |\n|---|---|---|\n"
        "| PBE0/def2-SVP | ethanol_pbe0 | -154.981234 |\n"
        "| B3LYP/def2-SVP | ethanol_b3lyp | -154.907142 |\n"
        "| PBE0/def2-TZVP | ethanol_pbe0_tzvp | -155.041188 |\n\n"
        "Best PBE0/def2-SVP minus best B3LYP/def2-SVP:\n"
        "-154.981234 - (-154.907142) = **0.074092 Eh** = **46.49 "
        "kcal/mol**.\n\n"
        "`methanol_pbe0` is a different molecule and is excluded from the "
        "ethanol comparison.\n"),
}

# The answers a broken capability produces. Each is the plausible wrong
# answer, not noise -- a test that only fails on garbage proves the
# verifier can read, not that it can discriminate.
_BAD = {
    "m2_a1_grep_file_where_is_it_defined": (
        "g1/answer.md",
        "# Files containing `_settle`\n\n"
        "Found `_settle` in: src/core.py, src/run.py and notes/README.md.\n"
        "There are two call sites: one in core.py and one in run.py.\n"),
    "m2_a2_bash_mean_cycle_time": (
        "b2/answer.md",
        "# Mean cycle time\n\n"
        "Averaging every row across the three runs gives a mean cycle "
        "time of 116.00 s. run_a looks like the fastest run overall.\n"),
    "m2_a3_read_file_two_spec_sheets": (
        "a3/answer.md",
        "# Pumps\n\n"
        "Pump A flows at 24.0 L/min. That is the highest flow rate of the "
        "sheets in the folder.\n"),
    "m2_a5_read_document_the_current_note": (
        "d4/answer.md",
        "# Method note\n\n"
        "Do not exceed 150 C. Report results in mg/L.\n"),
    "m2_a6_search_calcs_within_method": (
        "sv1/answer.md",
        "# Lowest energies\n\n"
        "The lowest calculation overall is methanol_pbe0 at -115.501000 "
        "Eh, so it wins. The ethanol runs are at -154.98 and -154.91 Eh, "
        "a difference of about 0.07 Eh, so PBE0 is lower.\n"),
}


@pytest.fixture(scope="module")
def tasks_by_id():
    tasks = [t for t in load_tasks() if t.id.startswith("m2_")]
    assert len(tasks) == 5, (
        f"expected the five m2 tasks, found {len(tasks)} -- the census "
        "counts on them existing")
    return {t.id: t for t in tasks}


def _workspace(tmp_path: Path, task) -> Path:
    """A workspace with the task's real fixture, built by its own setup."""
    ws = tmp_path / "ws"
    ws.mkdir(parents=True, exist_ok=True)
    if task.setup:
        ok, out = run_setup(task.setup, ws, root=ROOT)
        assert ok, f"setup did not build the fixture: {out}"
    return ws


@pytest.mark.parametrize("task_id", sorted(_GOOD))
def test_a_correct_answer_passes(tmp_path, tasks_by_id, task_id):
    task = tasks_by_id[task_id]
    ws = _workspace(tmp_path, task)
    rel, body = _GOOD[task_id]
    (ws / rel).write_text(body, encoding="utf-8")
    ok, out = run_acceptance(task.verify, ws, root=ROOT)
    assert ok is True, f"{task_id} rejected a correct answer:\n{out}"


@pytest.mark.parametrize("task_id", sorted(_BAD))
def test_the_wrong_answer_fails(tmp_path, tasks_by_id, task_id):
    """The half that makes the test a test.

    This is the answer the capability's failure mode produces. If it
    passed, the verifier would be measuring the presence of an
    answer.md, not the capability.
    """
    task = tasks_by_id[task_id]
    ws = _workspace(tmp_path, task)
    rel, body = _BAD[task_id]
    (ws / rel).write_text(body, encoding="utf-8")
    ok, out = run_acceptance(task.verify, ws, root=ROOT)
    assert ok is False, (
        f"{task_id} PASSED the wrong answer, so it measures nothing:\n"
        f"{body}\n-- verifier said:\n{out}")
    # And it failed for its own reason, not because it could not be
    # called: a usage error would fail every attempt regardless of merit.
    assert "usage" not in out.lower(), (
        f"{task.verify} answered with a usage error, so it is not callable "
        f"the way the runner calls it:\n{out}")


def test_a_missing_deliverable_is_a_failure_not_a_usage_error(tmp_path,
                                                             tasks_by_id):
    """An empty workspace must grade as a failed attempt, never as a
    script that could not run."""
    for task_id, task in sorted(tasks_by_id.items()):
        ws = _workspace(tmp_path / task_id, task)
        ok, out = run_acceptance(task.verify, ws, root=ROOT)
        assert ok is False, f"{task_id} passed with no answer written:\n{out}"
        assert "usage" not in out.lower(), (
            f"{task.verify} is not callable in the runner's form:\n{out}")


def test_the_setup_builds_what_the_tasks_ask_about(tmp_path, tasks_by_id):
    """The fixture the questions refer to must exist. A setup that
    silently built nothing would make every attempt an unmeasured one
    reported as a model failure."""
    ws = _workspace(tmp_path, tasks_by_id["m2_a1_grep_file_where_is_it_defined"])
    manifest = ws / "m2_fixture_manifest.json"
    assert manifest.is_file(), "the setup wrote no manifest"
    built = json.loads(manifest.read_text(encoding="utf-8"))["built"]
    assert sorted(built) == ["A1", "A2", "A3", "A4", "A5", "A6"]
    for rel in ("g1/src/core.py", "b2/run_a/cycle_times.csv",
                "a3/specs/pump_a.txt", "d4/method_note.docx",
                "calc_archive/calc/ethanol_pbe0/ethanol_pbe0.out"):
        assert (ws / rel).is_file(), f"setup did not write {rel}"


def test_the_docx_fixture_is_a_real_document(tmp_path, tasks_by_id):
    """The read_document task is only a measurement if its fixture is a
    document format a plain read cannot handle -- otherwise read_file
    would answer it and the task would measure read_file."""
    ws = _workspace(tmp_path, tasks_by_id["m2_a5_read_document_the_current_note"])
    docx = ws / "d4" / "method_note.docx"
    raw = docx.read_bytes()
    assert raw[:2] == b"PK", "the .docx is not a zip container"
    # The decisive property: the limit is NOT in the bytes a naive text
    # read would find, so an agent must actually parse the container.
    assert b"120 C" not in raw, (
        "the .docx stores its text uncompressed, so a plain read could "
        "recover it and the task would not measure read_document")


def test_every_python_file_the_setup_writes_is_valid_python(tmp_path, tasks_by_id):
    """The tree an agent searches must be a tree that could run.

    A real benchmark run found the ``g1`` fixture's ``src/core.py`` would
    not parse: its docstring sat at column 0 instead of under the ``def``,
    so ``ast.parse`` raised IndentationError and the ``from core import
    step, _settle`` in ``run.py`` would fail at runtime.  Nothing caught
    it because every verifier grades the answer TEXT.

    The task survives that -- it asks where a name is written, and a grep
    does not compile.  But a fixture that cannot be imported makes the
    wrong answers unrepresentative: an agent that tried to USE the tree
    rather than read it would get a different failure from the one the
    task is about.  The files are built by the setup script, so the
    property is checkable here rather than by eye.
    """
    ws = _workspace(tmp_path, tasks_by_id["m2_a1_grep_file_where_is_it_defined"])
    python_files = sorted(p for p in ws.rglob("*.py") if "m2_capability" not in p.name)
    assert python_files, "the setup wrote no Python file to check"
    broken = []
    for path in python_files:
        try:
            ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        except SyntaxError as exc:
            broken.append(f"{path.relative_to(ws)}: {exc.msg} (line {exc.lineno})")
    assert not broken, (
        "the fixture tree an agent must search is not valid Python, so its "
        "wrong answers are not the ones a real failure produces:\n  "
        + "\n  ".join(broken))
