"""An abandoned trial is not a failed calculation.

CONTROL.txt is where a DELFIN pipeline states its own plan:
``calc_initial``, ``oxidation_steps``, ``reduction_steps``. Those stages
ARE the grid, and the pipeline spells them ``initial``, ``ox_step_1``,
``red_step_2`` -- as a file or as the directory that owns the stage.

Two defects followed from the contract not reading that plan.

FIRST, ``grid_complete`` was unjudged for every pipeline run, because
``expected_steps`` had to be passed by a caller and the session-end
stage passed none. The number was in the run's own CONTROL all along.

SECOND, and worse, the flag counted EVERY step in the workspace. Measured
on a real 33-minute run: 14 ORCA outputs, 2 of them the production
TD-DFT calculations that terminated normally, 12 abandoned probes from
the agent trying settings. The contract reported "12 step(s) did not end
normally" and refused the grid. An agent that tries five settings and
keeps the one that converges is doing good science, and a check that
fires on correct work is a check somebody switches off.

So the grid is judged over the DECLARED stages, exploration is reported
beside it rather than counted against it, and where CONTROL declares
nothing every step is judged exactly as before -- without a declaration
there is no basis for the distinction and inventing one would be a
guess.

The guard that matters most is the third test: a CONTROL that declares
stages no file in the workspace carries must NOT be judged over the
empty set, which would come out complete.

The outputs here are written by this test. The repository cannot hold
real engine output -- LICENSE_GUARD.md, LGPL-3.0 against the ORCA EULA
-- which is also why the committed example pipelines carry their stage
directories and no calculations.
"""
from __future__ import annotations

import pytest

from delfin.agent import run_evidence
from delfin.config import declared_redox_steps, parse_control_text

_CONTROL = """\
input_file=input.txt
charge=+2
calc_initial=yes
oxidation_steps=
reduction_steps=1,2
method=
"""

_DONE = ("FINAL SINGLE POINT ENERGY   -76.12345678\n"
         "TOTAL RUN TIME: 0 days 0 hours 0 minutes 12.0 seconds\n"
         "                   ****ORCA TERMINATED NORMALLY****\n")
_CUT = ("               *** Restarting incremental Fock matrix formation\n")


def _stage(ws, name, text=_DONE, *, as_dir=True):
    """One stage, written the way the pipeline writes it."""
    folder = (ws / f"{name}_OCCUPIER") if as_dir else ws
    folder.mkdir(parents=True, exist_ok=True)
    (folder / f"{name}.inp").write_text("! PBE0 def2-SVP\n", encoding="utf-8")
    (folder / f"{name}.out").write_text(text, encoding="utf-8")


def _pipeline(tmp_path, control=_CONTROL):
    ws = tmp_path / "run"
    ws.mkdir()
    (ws / "CONTROL.txt").write_text(control, encoding="utf-8")
    return ws


def _contract(ws, **kw):
    return run_evidence.contract(run_evidence.collect(ws), **kw)


# --- the reader, pinned to the pipeline's own convention --------------

def test_the_reader_agrees_with_the_pipelines_own_parsing():
    """The same three lines live in cli_imag.py and cli_report.py. A
    third copy is how three answers to one question start to disagree,
    so the convention is asserted here against those callers' rule:
    calc_initial defaults to yes, the step keys are comma-separated
    indices, and a non-numeric entry is skipped rather than guessed."""
    cfg = parse_control_text(_CONTROL)
    assert declared_redox_steps(cfg) == ["initial", "red_step_1",
                                         "red_step_2"]
    cfg_ox = parse_control_text(
        "calc_initial=no\noxidation_steps=1,2\nreduction_steps=\n")
    assert declared_redox_steps(cfg_ox) == ["ox_step_1", "ox_step_2"]
    cfg_junk = parse_control_text(
        "calc_initial=yes\nreduction_steps=1,abc,3\n")
    assert declared_redox_steps(cfg_junk) == ["initial", "red_step_1",
                                              "red_step_3"]


def test_an_empty_control_declares_nothing_but_the_initial():
    assert declared_redox_steps({}) == ["initial"]
    assert declared_redox_steps({"calc_initial": "no"}) == []


# --- the grid is the declared one -------------------------------------

def test_the_required_number_comes_from_control(tmp_path):
    """No caller has to pass it, which is why it was never judged."""
    ws = _pipeline(tmp_path)
    for name in ("initial", "red_step_1", "red_step_2"):
        _stage(ws, name)
    con = _contract(ws)
    assert con.grid_complete is True
    assert any("taken from CONTROL" in n for n in con.notes)


def test_abandoned_exploration_does_not_fail_the_grid(tmp_path):
    """The measured case: the declared stages finished and a handful of
    trials were abandoned beside them."""
    ws = _pipeline(tmp_path)
    for name in ("initial", "red_step_1", "red_step_2"):
        _stage(ws, name)
    for probe in ("pA", "pB", "pC"):
        (ws / f"{probe}.inp").write_text("! PBE0\n", encoding="utf-8")
        (ws / f"{probe}.out").write_text(_CUT, encoding="utf-8")
    con = _contract(ws)
    assert con.grid_complete is True, con.problems
    assert not any("did not end normally" in p for p in con.problems)


def test_the_abandoned_trials_are_still_reported(tmp_path):
    """Reported, not hidden: a reader must see that three trials were
    left unfinished, even though they do not count against the grid."""
    ws = _pipeline(tmp_path)
    for name in ("initial", "red_step_1", "red_step_2"):
        _stage(ws, name)
    (ws / "pA.inp").write_text("! PBE0\n", encoding="utf-8")
    (ws / "pA.out").write_text(_CUT, encoding="utf-8")
    con = _contract(ws)
    note = " ".join(con.notes)
    assert "outside the declared stages" in note
    assert "pA.out" in note
    assert "abandoned is not a failed calculation" in note


def test_a_declared_stage_that_failed_still_fails_the_grid(tmp_path):
    """The flag must keep its teeth: the whole point is that a stage the
    run promised and did not deliver is a failure."""
    ws = _pipeline(tmp_path)
    _stage(ws, "initial")
    _stage(ws, "red_step_1")
    _stage(ws, "red_step_2", text=_CUT)
    con = _contract(ws)
    assert con.grid_complete is False
    assert any("did not end normally" in p for p in con.problems)
    assert any("red_step_2" in p for p in con.problems)


def test_a_missing_declared_stage_fails_the_grid(tmp_path):
    """Two of three stages ran. Every file present being fine is what
    two of three calculations also satisfies."""
    ws = _pipeline(tmp_path)
    _stage(ws, "initial")
    _stage(ws, "red_step_1")
    con = _contract(ws)
    assert con.grid_complete is False
    assert any("requires 3 calculation(s)" in p for p in con.problems)


# --- the guard that keeps this honest ---------------------------------

def test_a_declaration_nothing_matches_judges_every_step(tmp_path):
    """THE important guard. If the declared names matched nothing, the
    grid would be judged over an empty set and come out complete -- a
    workspace with a different naming would certify itself."""
    ws = _pipeline(tmp_path)
    for probe in ("alpha", "beta"):
        (ws / f"{probe}.inp").write_text("! PBE0\n", encoding="utf-8")
        (ws / f"{probe}.out").write_text(_CUT, encoding="utf-8")
    con = _contract(ws)
    assert con.grid_complete is False, (
        "a workspace whose files match no declared stage must not pass")
    assert any("carries their names" in n for n in con.notes)
    assert any("did not end normally" in p for p in con.problems)


def test_without_a_declaration_every_step_is_judged(tmp_path):
    """No CONTROL, so no basis to tell a trial from a stage. Behaves
    exactly as before this change."""
    ws = tmp_path / "bare"
    ws.mkdir()
    (ws / "a.inp").write_text("! PBE0\n", encoding="utf-8")
    (ws / "a.out").write_text(_CUT, encoding="utf-8")
    con = _contract(ws, expected_steps=1)
    assert con.grid_complete is False
    assert any("did not end normally" in p for p in con.problems)


def test_a_caller_that_passes_the_number_still_decides(tmp_path):
    """An explicit expectation is not overridden by CONTROL: a benchmark
    task that requires four calculations says so, and the run's own
    CONTROL does not get to lower the bar."""
    ws = _pipeline(tmp_path)
    for name in ("initial", "red_step_1", "red_step_2"):
        _stage(ws, name)
    con = _contract(ws, expected_steps=4)
    assert con.grid_complete is False
    assert any("requires 4 calculation(s)" in p for p in con.problems)


def test_a_stage_written_as_a_file_is_recognised_too(tmp_path):
    """The pipeline writes a stage either as red_step_1.out or inside
    red_step_1_OCCUPIER/; both belong to the stage."""
    ws = _pipeline(tmp_path)
    for name in ("initial", "red_step_1", "red_step_2"):
        _stage(ws, name, as_dir=False)
    con = _contract(ws)
    assert con.grid_complete is True, con.problems


def test_a_workspace_with_no_control_declares_nothing(tmp_path):
    """An absent CONTROL is not a plan. calc_initial defaults to "yes"
    for the pipeline, so asking the reader about an empty mapping would
    answer ["initial"] -- one stage, for a workspace that has no
    pipeline and no such file. An agent session would then be measured
    against a plan nobody wrote."""
    ws = tmp_path / "agentrun"
    ws.mkdir()
    (ws / "a.inp").write_text("! PBE0\n", encoding="utf-8")
    (ws / "a.out").write_text(_DONE, encoding="utf-8")
    con = _contract(ws)
    assert not any("taken from CONTROL" in n for n in con.notes)
    assert not any("CONTROL declares" in n for n in con.notes)
    # and with nothing declared and no expectation given, completeness
    # is not judged rather than guessed
    assert con.grid_complete is None


# --- blindness to the result is worse than a strict verdict ------------

def test_an_output_whose_input_has_another_name_is_still_a_step(tmp_path):
    """The measured case, and the most serious of the lot.

    Steps were anchored on .inp files alone. On a real run the agent
    edited `tddft_pbe0.inp` and wrote `tddft_pbe0_v2.out`, so BOTH
    production calculations -- the only two that terminated normally --
    were missing from the evidence, and the contract reported "12 of 12
    step(s) did not end normally" over a run that had delivered its
    results.

    That is not a strict judgement. It is blindness to the result, and
    it would void a correct run.
    """
    ws = tmp_path / "run"
    ws.mkdir()
    (ws / "job.inp").write_text("! PBE0 def2-SVP\n", encoding="utf-8")
    (ws / "job.out").write_text(_CUT, encoding="utf-8")        # abandoned
    (ws / "job_v2.out").write_text(_DONE, encoding="utf-8")    # the result
    ev = run_evidence.collect(ws)
    outputs = {s.get("output") for s in ev.steps}
    assert "job_v2.out" in outputs, (
        "the successful output is invisible to the contract")
    done = [s for s in ev.steps if s.get("terminated_normally") is True]
    assert [s["output"] for s in done] == ["job_v2.out"]


def test_such_a_step_says_its_input_is_not_known(tmp_path):
    """It can be judged and it cannot be re-executed, and the record has
    to say which: a step whose input cannot be pointed at must not look
    as though it could be run again."""
    ws = tmp_path / "run"
    ws.mkdir()
    (ws / "only_output.out").write_text(_DONE, encoding="utf-8")
    ev = run_evidence.collect(ws)
    step = next(s for s in ev.steps if s.get("output") == "only_output.out")
    assert step["input"] is None
    assert step.get("input_known") is False
    assert step["terminated_normally"] is True


def test_a_paired_output_is_not_collected_twice(tmp_path):
    """The .inp-anchored pass runs first, so an output that has its own
    input is one step and not two."""
    ws = tmp_path / "run"
    ws.mkdir()
    (ws / "job.inp").write_text("! PBE0\n", encoding="utf-8")
    (ws / "job.out").write_text(_DONE, encoding="utf-8")
    ev = run_evidence.collect(ws)
    assert [s.get("output") for s in ev.steps].count("job.out") == 1


def test_the_failure_message_names_what_succeeded_as_well(tmp_path):
    """"12 step(s) did not end normally" was literally true of a run
    whose production calculations both finished. A sentence listing only
    the failures invites the wrong conclusion from a reader who cannot
    see the directory."""
    ws = tmp_path / "run"
    ws.mkdir()
    (ws / "probe.inp").write_text("! PBE0\n", encoding="utf-8")
    (ws / "probe.out").write_text(_CUT, encoding="utf-8")
    (ws / "result_v2.out").write_text(_DONE, encoding="utf-8")
    con = run_evidence.contract(run_evidence.collect(ws), expected_steps=1)
    msg = " ".join(con.problems)
    assert "1 of 2" in msg
    assert "the ones that did" in msg and "result_v2.out" in msg
