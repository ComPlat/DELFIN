"""An OCCUPIER stage keeps what its first configuration computed, and no stage writes into another's files.

Two things wrote over finished results, both measured on the archive:

* A stage's first configuration starts from ``input0.xyz`` and ends in
  ``input.xyz`` / ``input.gbw``, the names the stage's setup also uses to hand
  the stage its start.  Set up again -- a recalc that reruns any FoB of the
  stage -- the setup put the start back over the finished geometry, and the
  comparison then handed that unoptimised start to the stage's frequency job.
  120 archived stages hold a start geometry in input.xyz next to a complete
  output.out; reproduced on NH3.
* The ox/red stages get the parent's preferred orbitals as ``input.gbw``, a
  link.  ORCA's results are copied back from its isolated directory, and the
  copy wrote through the link: the stage's orbitals went into the parent's
  file.  All 898 archived parent orbital files a later stage linked to were
  rewritten after their own job had ended.
"""

from __future__ import annotations

import os
from pathlib import Path

from delfin import copy_helpers, orca, thread_safe_helpers

OK = "\n                             ****ORCA TERMINATED NORMALLY****\n"
START = "N 0.0 0.0 0.0\nH 0.0 0.94 0.0\nH 0.81 -0.47 0.0\nH -0.81 -0.47 0.0\n"
FINISHED = ("4\nCoordinates from ORCA-job /scratch/run/initial_OCCUPIER/.orca_iso_input_1_2_x/input E -56.43\n"
            "N 0.0 0.0 0.1\nH 0.0 1.0 0.0\nH 0.87 -0.5 0.0\nH -0.87 -0.5 0.0\n")


def _job(tmp_path: Path) -> Path:
    job = tmp_path / "job"
    job.mkdir()
    (job / "start.txt").write_text(START)
    (job / "CONTROL.txt").write_text("charge=0\nmethod=OCCUPIER\ninput_file=input.txt\n")
    return job


def test_a_stage_set_up_again_keeps_its_finished_first_geometry(tmp_path):
    job = _job(tmp_path)
    stage = job / "initial_OCCUPIER"
    stage.mkdir()
    (stage / "input.xyz").write_text(FINISHED)
    (stage / "output.out").write_text(OK)

    copy_helpers.prepare_occ_folder_only_setup("initial_OCCUPIER", 0, parent_dir=job)

    assert (stage / "input.xyz").read_text() == FINISHED
    assert "N 0.0 0.0 0.0" in (stage / "input0.xyz").read_text()
    assert not (stage / "input.txt").exists()


def test_a_new_stage_is_handed_its_start_as_before(tmp_path):
    job = _job(tmp_path)

    stage = copy_helpers.prepare_occ_folder_only_setup("initial_OCCUPIER", 0, parent_dir=job)

    assert (stage / "input.xyz").read_text() == (stage / "input0.xyz").read_text()
    assert (stage / "input.xyz").read_text().startswith("4\n\nN 0.0 0.0 0.0")


def test_an_unfinished_first_configuration_is_handed_the_start_again(tmp_path):
    job = _job(tmp_path)
    stage = job / "initial_OCCUPIER"
    stage.mkdir()
    (stage / "input.xyz").write_text(FINISHED)
    (stage / "output.out").write_text("... killed at the walltime\n")

    copy_helpers.prepare_occ_folder_only_setup("initial_OCCUPIER", 0, parent_dir=job)

    assert "N 0.0 0.0 0.0" in (stage / "input.xyz").read_text()


def _redox_stage(tmp_path, monkeypatch, *, finished: bool):
    job = _job(tmp_path)
    parent = job / "initial_OCCUPIER"
    parent.mkdir()
    (parent / "input.gbw").write_bytes(b"parent orbitals")
    (job / "input_initial_OCCUPIER.xyz").write_text(FINISHED.replace("input E", "input E -1"))
    stage = job / "ox_step_1_OCCUPIER"
    stage.mkdir()
    if finished:
        (stage / "input.xyz").write_text(FINISHED.replace("initial_OCCUPIER", "ox_step_1_OCCUPIER"))
        (stage / "input.gbw").write_bytes(b"the stage's own orbitals")
        (stage / "output.out").write_text(OK)
    monkeypatch.setattr(thread_safe_helpers, "read_occupier_file_threadsafe",
                        lambda *a, **k: (1, "", 1, None))
    thread_safe_helpers.prepare_occ_folder_2_only_setup("ox_step_1_OCCUPIER", "initial_OCCUPIER", 1,
                                                        {"OCCUPIER_method": "manually"}, job)
    return stage


def test_a_redox_stage_set_up_again_keeps_its_own_geometry_and_orbitals(tmp_path, monkeypatch):
    stage = _redox_stage(tmp_path, monkeypatch, finished=True)

    assert "ox_step_1_OCCUPIER" in (stage / "input.xyz").read_text()
    assert not (stage / "input.gbw").is_symlink()
    assert (stage / "input.gbw").read_bytes() == b"the stage's own orbitals"


def test_a_new_redox_stage_is_still_handed_the_parents_orbitals(tmp_path, monkeypatch):
    stage = _redox_stage(tmp_path, monkeypatch, finished=False)

    assert (stage / "input.gbw").read_bytes() == b"parent orbitals"


def _isolated_run(tmp_path, monkeypatch, make_link):
    parent = tmp_path / "initial_OCCUPIER"
    parent.mkdir()
    (parent / "input.gbw").write_bytes(b"parent orbitals")
    stage = tmp_path / "ox_step_1_OCCUPIER"
    stage.mkdir()
    make_link(parent / "input.gbw", stage / "input.gbw")
    inp = stage / "input.inp"
    inp.write_text("! PBE0 def2-SVP OPT\n* xyz 1 2\nN 0 0 0\n*\n")

    def fake(orca_path, input_file_path, output_log, timeout=None, scratch_subdir=None,
             working_dir=None, extra_env=None):
        (Path(working_dir) / "input.gbw").write_bytes(b"the stage's own orbitals")
        Path(output_log).write_text(OK)
        return True

    monkeypatch.setattr(orca, "_run_orca_subprocess", fake)
    assert orca._run_orca_isolated("/bin/true", inp, stage / "output.out")
    return parent, stage


def test_a_result_copied_back_does_not_write_through_a_symlink(tmp_path, monkeypatch):
    parent, stage = _isolated_run(tmp_path, monkeypatch,
                                  lambda src, dst: dst.symlink_to(os.path.relpath(src, dst.parent)))

    assert (parent / "input.gbw").read_bytes() == b"parent orbitals"
    assert not (stage / "input.gbw").is_symlink()
    assert (stage / "input.gbw").read_bytes() == b"the stage's own orbitals"


def test_nor_through_a_hard_link(tmp_path, monkeypatch):
    parent, stage = _isolated_run(tmp_path, monkeypatch, lambda src, dst: os.link(src, dst))

    assert (parent / "input.gbw").read_bytes() == b"parent orbitals"
    assert (stage / "input.gbw").read_bytes() == b"the stage's own orbitals"
