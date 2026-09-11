"""A CONTROL edit reaches the jobs it changes, and only those.

A recalc keeps the jobs it finds finished.  Which of them an edit makes out
of date is told against the CONTROL they were computed with: a completed run
records it, and the Recalc tab records the file it replaces when no run has.
Measured before: the Recalc tab submitted a classic recalc, which keeps every
finished job whatever the edit; and OCCUPIER kept its inputs in any recalc,
so a new functional reached none of its jobs.  Now, after an edit that
reaches the ORCA inputs, OCCUPIER writes its inputs anew; one that says what
the finished job ran keeps the job, one that says something else is a new
calculation.  Without an edit nothing is written anew, as before.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin import recalc_control, smart_recalc
from delfin.config import set_control_value
from delfin.define import TEMPLATE
from delfin.workflows import pipeline

CONTROL = set_control_value(set_control_value(TEMPLATE, "charge", "0"), "method", "OCCUPIER")


# ------------------------------------------------------------------ what changed

def test_resources_are_not_a_change_and_the_level_of_theory_is():
    assert not recalc_control.ControlChange(keys=recalc_control.changed_keys(
        CONTROL, set_control_value(set_control_value(CONTROL, "PAL", "96"), "maxcore", "9000"))).computation
    change = recalc_control.ControlChange(keys=recalc_control.changed_keys(
        CONTROL, set_control_value(CONTROL, "functional", "B3LYP")))
    assert change.computation and not change.structure


def test_what_builds_the_structure_is_a_structure_change():
    for key, value in (("smiles_converter", "QUICK"), ("MANTA_SEEDS", "7"), ("global_optimizer", "CREST"),
                       ("charge", "1"), ("solvent", "toluene")):
        keys = recalc_control.changed_keys(CONTROL, set_control_value(CONTROL, key, value))
        assert recalc_control.ControlChange(keys=keys).structure, key


def test_an_edit_reaches_only_the_occupier_inputs_written_from_it():
    def change(**values):
        text = CONTROL
        for key, value in values.items():
            text = set_control_value(text, key, value)
        return recalc_control.ControlChange(keys=recalc_control.changed_keys(CONTROL, text))

    esd = change(ESD_LINEW="60", TDDFT_nroots="20", IMAG_max_rounds="3")
    assert esd.computation and not esd.reaches_occupier_frequency_jobs and not esd.reaches_occupier_fobs
    freq = change(freq_type="numFREQ")
    assert freq.reaches_occupier_frequency_jobs and not freq.reaches_occupier_fobs
    functional = change(functional="B3LYP")
    assert functional.reaches_occupier_fobs and functional.reaches_occupier_frequency_jobs


@pytest.mark.parametrize("target, fobs, frequency_jobs", [
    ("initial", False, True), ("ox_step_1", False, True), ("input5", True, False),
    ("red_step_2_OCCUPIER", True, True), ("all", True, True), ("S1", False, False),
])
def test_an_override_reaches_the_jobs_it_names(target, fobs, frequency_jobs):
    change = recalc_control.ControlChange(keys=recalc_control.changed_keys(
        CONTROL, CONTROL + f"\nkeyword:{target}=[TightSCF]\n"))

    assert (change.reaches_occupier_fobs, change.reaches_occupier_frequency_jobs) == (fobs, frequency_jobs)


def test_a_key_an_older_file_lacks_counts_with_its_default():
    older = "\n".join(line for line in CONTROL.splitlines() if not line.startswith("MANTA_"))

    assert recalc_control.changed_keys(older, CONTROL) == frozenset()


def test_a_completed_run_is_what_the_next_recalc_is_told_against(tmp_path):
    recalc_control.record_completed_run(tmp_path, CONTROL, "N 0 0 0\n")

    assert recalc_control.change_since_last_run(tmp_path, CONTROL, "N 0 0 0\n") == \
        recalc_control.ControlChange()
    edited = recalc_control.change_since_last_run(tmp_path, set_control_value(CONTROL, "functional", "TPSS"),
                                                  "N 0 0 0\n")
    assert edited.keys == {"functional"}
    moved = recalc_control.change_since_last_run(tmp_path, CONTROL, "N 0 0 0.5\n")
    assert moved.input_changed and moved.structure


def test_without_a_record_the_change_is_unknown(tmp_path):
    assert recalc_control.change_since_last_run(tmp_path, CONTROL) is None


def test_the_recalc_tab_records_what_it_replaces_only_when_no_run_did(tmp_path):
    assert recalc_control.remember_before_edit(tmp_path, CONTROL)
    assert json.loads((tmp_path / recalc_control.RECORD_NAME).read_text())["by"] == "recalc edit"

    recalc_control.record_completed_run(tmp_path, set_control_value(CONTROL, "functional", "TPSS"))
    assert not recalc_control.remember_before_edit(tmp_path, CONTROL)
    assert "TPSS" in recalc_control.previous_run(tmp_path)["control"]


# ------------------------------------------------------------ the built structure

def _handed_on(job: Path, smiles: str | None = None) -> None:
    (job / "start.txt").write_text("N 0 0 0\n")
    (job / "initial_OCCUPIER").mkdir()
    if smiles:
        (job / "guppy_input.txt").write_text(smiles + "\n")


@pytest.fixture
def recalc(monkeypatch):
    monkeypatch.setenv("DELFIN_RECALC", "1")


def test_a_recalc_keeps_the_structure_when_nothing_that_builds_it_changed(tmp_path, recalc):
    _handed_on(tmp_path)
    config = {"_recalc_change": recalc_control.ControlChange(keys=frozenset({"functional"}))}

    assert pipeline._structure_kept_by_recalc(config, tmp_path, None)


def test_it_builds_it_again_when_the_edit_reaches_it(tmp_path, recalc):
    _handed_on(tmp_path)
    config = {"_recalc_change": recalc_control.ControlChange(keys=frozenset({"MANTA_SEEDS"})),
              "_recalc_rebuild_structure": True}

    assert not pipeline._structure_kept_by_recalc(config, tmp_path, None)


def test_a_classic_recalc_keeps_it_whatever_the_edit(tmp_path, recalc):
    _handed_on(tmp_path)
    config = {"_recalc_change": recalc_control.ControlChange(keys=frozenset({"MANTA_SEEDS"}))}

    assert pipeline._structure_kept_by_recalc(config, tmp_path, None)


def test_without_a_record_manta_keeps_the_structure_it_built_for_this_smiles(tmp_path, recalc):
    _handed_on(tmp_path, smiles="CCO")
    config = {"smiles_converter": "MANTA"}

    assert pipeline._structure_kept_by_recalc(config, tmp_path, "CCO")
    assert not pipeline._structure_kept_by_recalc(config, tmp_path, "CCN")


def test_without_a_record_an_xyz_input_is_read_again_as_before(tmp_path, recalc):
    _handed_on(tmp_path)

    assert not pipeline._structure_kept_by_recalc({}, tmp_path, None)


def test_a_first_run_builds_it(tmp_path, monkeypatch):
    monkeypatch.setenv("DELFIN_RECALC", "0")
    _handed_on(tmp_path)
    config = {"_recalc_change": recalc_control.ControlChange()}

    assert not pipeline._structure_kept_by_recalc(config, tmp_path, None)


# --------------------------------------------------------------- an input anew

INPUT = "! PBE0 def2-SVP OPT\n%pal nprocs 4 end\n%maxcore 2000\n* xyz 0 1\nN 0 0 0\n*\n"


def test_an_input_written_anew_that_says_the_same_leaves_the_old_file(tmp_path):
    inp = tmp_path / "input.inp"
    inp.write_text(INPUT)
    smart_recalc.store_fingerprint(inp)
    before = smart_recalc.snapshot_input(inp)
    inp.write_text(INPUT.replace("nprocs 4", "nprocs 40").replace("  ", " "))

    assert smart_recalc.settle_rewritten_input(inp, before, tmp_path)
    assert inp.read_text() == INPUT
    assert smart_recalc.fingerprint_unchanged(inp)


def test_one_that_says_something_else_stays(tmp_path):
    inp = tmp_path / "input.inp"
    inp.write_text(INPUT)
    before = smart_recalc.snapshot_input(inp)
    inp.write_text(INPUT.replace("PBE0", "B3LYP"))

    assert not smart_recalc.settle_rewritten_input(inp, before, tmp_path)
    assert "B3LYP" in inp.read_text()


def test_the_overrides_are_merged_before_comparing(tmp_path):
    (tmp_path / "CONTROL.txt").write_text("keyword:input=[TightOpt]\n")
    inp = tmp_path / "input.inp"
    inp.write_text(INPUT.replace("OPT", "TightOpt"))   # the finished run merged it
    before = smart_recalc.snapshot_input(inp)
    inp.write_text(INPUT)                             # written anew from CONTROL

    assert smart_recalc.settle_rewritten_input(inp, before, tmp_path)


def test_a_new_calculation_does_not_inherit_the_old_ones_retries(tmp_path):
    from delfin.orca_recovery import OrcaErrorType, RetryStateTracker

    inp = tmp_path / "input5.inp"
    inp.write_text(INPUT)
    (tmp_path / "input5.retry1.inp").write_text(INPUT)
    (tmp_path / "input5.retry1.inp.fprint").write_text("x\n")
    tracker = RetryStateTracker(tmp_path / ".delfin_recovery_state.json")
    for _ in range(3):
        tracker.increment_attempt("input5", OrcaErrorType.GEOMETRY_NOT_CONVERGED)

    smart_recalc.forget_earlier_attempts(inp)

    assert not list(tmp_path.glob("input5.retry*"))
    assert RetryStateTracker(tmp_path / ".delfin_recovery_state.json").should_retry(
        "input5", OrcaErrorType.GEOMETRY_NOT_CONVERGED, 3)


# ------------------------------------------------------- OCCUPIER after an edit

from delfin import occupier_flat_extraction as ofe  # noqa: E402

OK = "\nFINAL SINGLE POINT ENERGY      {e}\n\n                             ****ORCA TERMINATED NORMALLY****\n"
SEQUENCE = [
    {"index": 1, "m": 1, "BS": ""},
    {"index": 2, "m": 3, "BS": "", "from": 1},
    {"index": 3, "m": 3, "BS": "3,1", "from": 2},
]


def _occupier_stage(tmp_path, monkeypatch, functional_now: str):
    """A finished stage written with PBE0; the input generator now writes *functional_now*."""
    monkeypatch.setenv("DELFIN_RECALC", "1")
    monkeypatch.setenv("DELFIN_SMART_RECALC", "1")
    stage = tmp_path / "initial_OCCUPIER"
    stage.mkdir()

    def write(functional, name, mult):
        return f"! {functional} def2-SVP OPT\n%pal nprocs 4 end\n* xyz 0 {mult}\nN 0 0 0\n*\n"

    for entry in SEQUENCE:
        idx = entry["index"]
        stem = "input" if idx == 1 else f"input{idx}"
        (stage / f"{stem}.inp").write_text(write("PBE0", stem, entry["m"]))
        (stage / f"{stem}.xyz").write_text("1\n\nN 0 0 0\n")
        (stage / ("output.out" if idx == 1 else f"output{idx}.out")).write_text(OK.format(e=-56.4 + idx / 100))
        smart_recalc.store_fingerprint(stage / f"{stem}.inp")
    (stage / "input3.retry1.inp").write_text("stale retry\n")
    (stage / "OCCUPIER.txt").write_text("from the first run\n")

    def generator(src_idx, inp_name, charge, mult, *args, **kwargs):
        Path(inp_name).write_text(write(functional_now, inp_name, mult))

    ran = []

    def orca(inp, out, **kwargs):
        ran.append(Path(inp).name)
        Path(out).write_text(OK.format(e=-57.0))
        return True

    monkeypatch.setattr(ofe, "read_and_modify_file_OCCUPIER", generator)
    monkeypatch.setattr(ofe, "run_orca_with_intelligent_recovery", orca)
    monkeypatch.setattr(ofe, "_wait_for_geometry_source", lambda *a, **k: None)
    monkeypatch.setattr(ofe, "check_and_warn_competing_processes", lambda *a, **k: None)
    jobs, best_id = ofe._create_occupier_fob_jobs(
        folder_name="initial_OCCUPIER", folder_path=stage, stage_prefix="initial",
        sequence=SEQUENCE, sequence_label="even_seq", total_cores=8,
        global_config={"OCCUPIER_compare": "FSPE", "_recalc_regenerate_fob_inputs": True,
                       "enable_auto_recovery": "yes"},
        ensure_setup=lambda: stage, source_folder=None, metals=[], metal_basisset=None,
        main_basisset="def2-SVP", solvent="", occ_results={}, stage_charge=0,
    )
    return stage, {j.job_id: j for j in jobs}, ran


def test_after_an_edit_that_leaves_the_fob_inputs_alone_every_fob_is_kept(tmp_path, monkeypatch):
    stage, jobs, ran = _occupier_stage(tmp_path, monkeypatch, functional_now="PBE0")
    before = {p.name: p.read_bytes() for p in stage.glob("input*.inp")}

    assert not any(jobs[f"initial_fob_{i}"].precomplete_check() for i in (1, 2, 3))
    for i in (1, 2, 3):
        jobs[f"initial_fob_{i}"].work(8)

    assert ran == []
    assert {p.name: p.read_bytes() for p in stage.glob("input*.inp")} == before


def test_after_an_edit_that_changes_them_every_fob_is_computed_again(tmp_path, monkeypatch):
    stage, jobs, ran = _occupier_stage(tmp_path, monkeypatch, functional_now="B3LYP")

    for i in (1, 2, 3):
        jobs[f"initial_fob_{i}"].work(8)

    assert sorted(ran) == ["input.inp", "input2.inp", "input3.inp"]
    assert "B3LYP" in (stage / "input3.inp").read_text()
    assert not (stage / "input3.retry1.inp").exists(), "the old calculation's retry is not this one's"


def test_the_frequency_job_starts_where_the_finished_one_started(tmp_path):
    from delfin.workflows.engine import occupier as engine

    own = tmp_path / "ox_step_1.xyz"
    own.write_text("1\nCoordinates from ORCA-job /s/.orca_iso_ox_step_1_1_2_x/ox_step_1 E -56.2\nN 0 0 0\n")
    imag = tmp_path / "red_step_1.xyz"
    imag.write_text("1\nCoordinates from ORCA-job red_step_1.imag1 E -56.2\nN 0 0 0\n")
    handed = tmp_path / "ox_step_2.xyz"
    handed.write_text("1\nCoordinates from ORCA-job input2 E -56.2\nN 0 0 0\n")

    assert engine._written_by_stage_job(own, "ox_step_1")
    assert engine._written_by_stage_job(imag, "red_step_1")
    assert not engine._written_by_stage_job(handed, "ox_step_2")


def test_a_smart_recalc_is_told_what_the_edit_reaches(tmp_path, monkeypatch):
    from delfin import cli

    monkeypatch.setenv("DELFIN_SMART_RECALC", "1")
    monkeypatch.setenv("DELFIN_RECALC_CONTROL_EDITED", "0")   # restored after the test
    recalc_control.record_completed_run(tmp_path, CONTROL, "N 0 0 0\n")
    config = {}
    cli._note_control_change(config, tmp_path, set_control_value(CONTROL, "freq_type", "numFREQ"), "N 0 0 0\n")

    assert config["_recalc_regenerate_main_inputs"] and not config["_recalc_regenerate_fob_inputs"]
    assert not config["_recalc_rebuild_structure"]
    assert smart_recalc.control_edited()


def test_a_classic_recalc_is_not(tmp_path, monkeypatch):
    from delfin import cli

    monkeypatch.setenv("DELFIN_SMART_RECALC", "0")
    monkeypatch.setenv("DELFIN_RECALC_CONTROL_EDITED", "0")
    recalc_control.record_completed_run(tmp_path, CONTROL, "N 0 0 0\n")
    config = {}
    cli._note_control_change(config, tmp_path, set_control_value(CONTROL, "functional", "B3LYP"), "N 0 0 0\n")

    assert "_recalc_regenerate_fob_inputs" not in config
    assert not smart_recalc.control_edited()
