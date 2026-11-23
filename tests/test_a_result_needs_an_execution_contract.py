"""A report is not a result: the calculation has to have happened.

Input: workspaces as a run leaves them behind, read by
``run_evidence.collect``, judged by ``run_evidence.contract`` and by the
grader that the benchmark acceptance path calls
(``accept/w20_common.grade_run``).

Output measured: the four flags of the execution contract, and the
verdict of the grader with the contract mounted.

Why this exists, measured elsewhere and reproduced here in the first
test: a rubric can score a report whose calculation never finished. In
the published Quantum ESPRESSO study (arXiv:2611.10629) an execution
audit voided 1.6 % of runs, and 20 of the 26 voided runs had a raw score
of 0.85 or better; their example is a seven-point curve that met the
numerical tolerance while only three of the seven calculations
completed. This repository had the same hole: ``all_steps_terminated``
asks whether every step that EXISTS ended normally, which three finished
calculations out of seven satisfy.

The flags are True (shown), False (contradicted by the files) or None
(the workspace holds no evidence either way). None is never a pass.
"""
from __future__ import annotations

import os
import time

from delfin.agent import run_evidence
from delfin.agent.pack.benchmark.accept import w20_common


BANNER = """
                 *****************************************
                 *                 ORCA                  *
                 *****************************************

FINAL SINGLE POINT ENERGY      -76.026760543001

                             ****ORCA TERMINATED NORMALLY****
TOTAL RUN TIME: 0 days 0 hours 0 minutes 12.0 seconds
"""

UNFINISHED = """
                 *****************************************
                 *                 ORCA                  *
                 *****************************************
SCF ITERATIONS
  1    -75.9    0.02
  2    -76.0    0.00
"""


def _step(ws, name: str, *, output: str = BANNER) -> None:
    """One ORCA step in the workspace: its input and its output."""
    (ws / f"{name}.inp").write_text("! PBE0 def2-TZVP\n*xyz 0 1\nH 0 0 0\n*\n")
    if output is not None:
        (ws / f"{name}.out").write_text(output)


# ---------------------------------------------------------------------
# the hole this change closes
# ---------------------------------------------------------------------

def test_three_of_seven_calculations_is_not_a_complete_grid(tmp_path):
    """The measured failure mode: the points that ran all ended normally,
    and the task needed seven.

    The control in the same test is the point: the OLD judgement
    (``all_steps_terminated``) is satisfied by this workspace, so the
    report would have been graded on three sevenths of its evidence.
    """
    for i in range(3):
        _step(tmp_path, f"eos_{i}")
    ev = run_evidence.collect(tmp_path)

    # control: the pre-existing check cannot see the four that are missing
    assert w20_common.all_steps_terminated(ev) is True, (
        "every step present ended normally -- that is exactly why this "
        "check alone cannot judge completeness")

    con = run_evidence.contract(ev, expected_steps=7)
    assert con.grid_complete is False
    assert con.failed() is True
    assert any("requires 7" in p for p in con.problems), con.problems

    graded = w20_common.grade_run(ev, expected_steps=7)
    assert graded["verdict"] == "fail", graded
    assert any("execution contract" in p for p in graded["problems"])


def test_the_same_workspace_passes_when_the_task_needed_three(tmp_path):
    """Control: the grid is complete when three is what the task asked
    for. The flag judges the task's requirement, not a fixed number."""
    for i in range(3):
        _step(tmp_path, f"eos_{i}")
    ev = run_evidence.collect(tmp_path)
    con = run_evidence.contract(ev, expected_steps=3)
    assert con.grid_complete is True
    assert con.failed() is False
    assert w20_common.grade_run(ev, expected_steps=3)["verdict"] == "pass"


def test_without_a_declared_requirement_completeness_is_not_judged(tmp_path):
    """No requirement declared: the flag is None and says so. It does not
    become True because nothing contradicted it."""
    for i in range(3):
        _step(tmp_path, f"eos_{i}")
    ev = run_evidence.collect(tmp_path)
    con = run_evidence.contract(ev)
    assert con.grid_complete is None
    assert con.certified() is False
    assert "grid_complete" in con.unjudged()
    assert any("declared no required number" in n for n in con.notes), con.notes
    # ... and an unjudged flag does not fail the run either: it is
    # reported, which is what lets a reader see what was not checked.
    graded = w20_common.grade_run(ev)
    assert graded["verdict"] == "pass"
    assert "grid_complete" in graded["unjudged"]


# ---------------------------------------------------------------------
# engine_invoked
# ---------------------------------------------------------------------

def test_a_result_claimed_with_no_engine_output_is_void(tmp_path):
    """A workspace with a report and nothing else: no flag can be shown,
    and the claim is contradicted rather than unproven."""
    (tmp_path / "result.json").write_text('{"energy_hartree": -76.02676}')
    (tmp_path / "answer.txt").write_text("The energy is -76.02676 Hartree.\n")
    ev = run_evidence.collect(tmp_path)
    con = run_evidence.contract(ev, expected_steps=1)
    assert con.engine_invoked is False
    assert con.failed() is True
    assert w20_common.grade_run(ev, expected_steps=1)["verdict"] == "fail"


def test_an_input_without_an_output_does_not_show_a_calculation(tmp_path):
    """The input file is the plan, not the evidence."""
    (tmp_path / "opt.inp").write_text("! PBE0 def2-TZVP\n*xyz 0 1\nH 0 0 0\n*\n")
    ev = run_evidence.collect(tmp_path)
    con = run_evidence.contract(ev, expected_steps=1)
    assert con.engine_invoked is False
    assert con.grid_complete is False


def test_an_unfinished_output_is_not_a_finished_step(tmp_path):
    """An output that stops mid-SCF: the step is contradicted, not
    unknown."""
    _step(tmp_path, "opt", output=UNFINISHED)
    ev = run_evidence.collect(tmp_path)
    con = run_evidence.contract(ev, expected_steps=1)
    assert con.grid_complete is False
    assert con.engine_invoked is True, (
        "the engine did write an output; what failed is the grid")


# ---------------------------------------------------------------------
# wallclock_plausible
# ---------------------------------------------------------------------

def test_an_output_older_than_its_input_did_not_come_from_it(tmp_path):
    """An output carried in from elsewhere: the input that is supposed to
    have produced it is newer than the output."""
    _step(tmp_path, "opt")
    old = time.time() - 7200
    os.utime(tmp_path / "opt.out", (old, old))
    ev = run_evidence.collect(tmp_path)
    con = run_evidence.contract(ev, expected_steps=1)
    assert con.wallclock_plausible is False
    assert any("older than" in p for p in con.problems), con.problems
    assert w20_common.grade_run(ev, expected_steps=1)["verdict"] == "fail"


def test_more_run_time_than_the_files_span_is_contradicted(tmp_path):
    """The engine states six hours; the workspace's files were all
    written within a minute of each other."""
    _step(tmp_path, "opt", output=BANNER.replace(
        "0 days 0 hours 0 minutes 12.0 seconds",
        "0 days 6 hours 0 minutes 0.0 seconds"))
    ev = run_evidence.collect(tmp_path)
    con = run_evidence.contract(ev, expected_steps=1)
    assert con.wallclock_plausible is False
    # The message names the step rather than a total: the check takes the
    # LONGEST step and never the sum, because DELFIN runs engine jobs in
    # parallel and a sum would fail an honest run.
    assert any("did not fit in this workspace" in p for p in con.problems), \
        con.problems
    assert any("opt.out" in p for p in con.problems), con.problems


def test_a_stated_run_time_inside_the_span_is_plausible(tmp_path):
    """Control: twelve seconds of stated run time against files written
    minutes apart is no contradiction."""
    _step(tmp_path, "opt")
    old = time.time() - 600
    os.utime(tmp_path / "opt.inp", (old, old))
    ev = run_evidence.collect(tmp_path)
    con = run_evidence.contract(ev, expected_steps=1)
    assert con.wallclock_plausible is True


def test_without_a_stated_run_time_the_clock_is_not_judged(tmp_path):
    """An xtb log carries no ORCA run-time line and this workspace has no
    input beside it, so the clock cannot be judged -- and is reported as
    unjudged rather than as fine."""
    (tmp_path / "xtb_opt.out").write_text(
        "normal termination of xtb\nrc=0\n")
    ev = run_evidence.collect(tmp_path)
    con = run_evidence.contract(ev, expected_steps=1)
    assert con.wallclock_plausible is None
    assert "wallclock_plausible" in con.unjudged()


# ---------------------------------------------------------------------
# no_egress
# ---------------------------------------------------------------------

def test_egress_is_unjudged_unless_the_records_are_handed_over(tmp_path):
    """The files of a workspace do not record network activity. The flag
    says so instead of returning True."""
    _step(tmp_path, "opt")
    con = run_evidence.contract(run_evidence.collect(tmp_path),
                                expected_steps=1)
    assert con.no_egress is None
    assert any("do not record network" in n for n in con.notes), con.notes


def test_an_empty_egress_log_certifies_and_a_record_voids(tmp_path):
    """An empty list means the log was read and holds nothing, which is
    evidence; a record in the run's window voids the run."""
    _step(tmp_path, "opt")
    ev = run_evidence.collect(tmp_path)
    clean = run_evidence.contract(ev, expected_steps=1, egress_events=[])
    assert clean.no_egress is True
    assert clean.certified() is True, clean.flags

    dirty = run_evidence.contract(
        ev, expected_steps=1,
        egress_events=[{"kind": "outbound", "detail": "https://example.invalid"}])
    assert dirty.no_egress is False
    assert dirty.failed() is True
    graded = w20_common.grade_run(
        ev, expected_steps=1, egress_events=[{"kind": "outbound"}])
    assert graded["verdict"] == "fail"
    assert any("isolated run" in p for p in graded["problems"])


# ---------------------------------------------------------------------
# the honesty rule of the flags
# ---------------------------------------------------------------------

def test_an_unjudged_flag_is_never_read_as_a_pass():
    """The rule that makes the contract worth having: absence of
    evidence is not compliance."""
    con = run_evidence.Contract()
    assert con.certified() is False
    assert con.failed() is False
    assert sorted(con.unjudged()) == [
        "engine_invoked", "grid_complete", "no_egress",
        "wallclock_plausible"]
