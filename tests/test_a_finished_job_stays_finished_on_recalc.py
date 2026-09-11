"""A smart recalc computes what changed, and a finished job whose input did not change stays finished.

Measured on a finished formaldehyde ESD run from the dashboard, recalc with
nothing changed: S0 was computed again, and then every state and rate after
it.  The S0 input is written anew on every run, from initial.xyz when there
is one -- and after an ESD run initial.xyz is S0's own optimised geometry,
copied there for the redox steps.  So S0's input was its own last result.
298 of 300 archived initial.xyz carry "Coordinates from ORCA-job S0".

The fingerprint also held every input byte, the cores (%pal) included, and
absolute dependency paths: a job that got other cores, or a job folder that
moved, ran again.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path

from delfin import smart_recalc
from delfin.esd_input_generator import written_by_s0

INPUT = """! PBE0 def2-SVP OPT
%moinp "S0.gbw"
%pal nprocs 4 end
%maxcore 2000
* xyzfile 0 1 S0.xyz
"""


def _job(tmp_path: Path) -> Path:
    job = tmp_path / "job"
    job.mkdir()
    (job / "S1.inp").write_text(INPUT)
    (job / "S0.gbw").write_bytes(b"\0" * 100)
    (job / "S0.xyz").write_text("1\n\nH 0 0 0\n")
    return job


def test_the_cores_a_job_gets_do_not_change_its_fingerprint(tmp_path):
    job = _job(tmp_path)
    smart_recalc.store_fingerprint(job / "S1.inp")

    (job / "S1.inp").write_text(INPUT.replace("nprocs 4", "nprocs 16").replace("%maxcore 2000", "%maxcore 3800"))
    assert smart_recalc.fingerprint_unchanged(job / "S1.inp")

    (job / "S1.inp").write_text(INPUT.replace("OPT", "TightOpt"))
    assert not smart_recalc.fingerprint_unchanged(job / "S1.inp")


def test_a_job_folder_that_moved_keeps_its_fingerprint(tmp_path):
    job = _job(tmp_path)
    smart_recalc.store_fingerprint(job / "S1.inp")
    moved = tmp_path / "archive" / "job"
    moved.parent.mkdir()
    shutil.copytree(job, moved)          # copy2: sizes and mtimes kept

    assert smart_recalc.fingerprint_unchanged(moved / "S1.inp")


def test_a_changed_dependency_still_counts(tmp_path):
    job = _job(tmp_path)
    smart_recalc.store_fingerprint(job / "S1.inp")
    os.utime(job / "S0.gbw", ns=(1, 1))

    assert not smart_recalc.fingerprint_unchanged(job / "S1.inp")


def test_a_fingerprint_an_earlier_delfin_wrote_still_counts(tmp_path):
    job = _job(tmp_path)
    (job / "S1.inp.fprint").write_text(smart_recalc._legacy_fingerprint(job / "S1.inp") + "\n")

    assert smart_recalc.fingerprint_unchanged(job / "S1.inp")


def test_a_recalc_rewrites_an_old_fingerprint_so_the_folder_can_move_after_it(tmp_path, monkeypatch):
    # measured: a job skipped on its first recalc kept the old sidecar, and the
    # folder copied elsewhere afterwards computed S0 and every job after it again
    monkeypatch.setenv("DELFIN_RECALC", "1")
    monkeypatch.setenv("DELFIN_SMART_RECALC", "1")
    job = _job(tmp_path)
    (job / "S1.out").write_text("...\n****ORCA TERMINATED NORMALLY****\n")
    (job / "S1.inp.fprint").write_text(smart_recalc._legacy_fingerprint(job / "S1.inp") + "\n")

    assert smart_recalc.should_skip(job / "S1.inp", job / "S1.out", required_outputs=[])
    assert (job / "S1.inp.fprint").read_text().strip() == smart_recalc.compute_fingerprint(job / "S1.inp")

    moved = tmp_path / "moved" / "job"
    moved.parent.mkdir()
    shutil.copytree(job, moved)
    assert smart_recalc.should_skip(moved / "S1.inp", moved / "S1.out", required_outputs=[])


def test_looking_at_a_fingerprint_does_not_write_one(tmp_path):
    # the dashboard asks fingerprint_unchanged() about folders it only shows
    job = _job(tmp_path)
    old = smart_recalc._legacy_fingerprint(job / "S1.inp") + "\n"
    (job / "S1.inp.fprint").write_text(old)

    assert smart_recalc.fingerprint_unchanged(job / "S1.inp")
    assert (job / "S1.inp.fprint").read_text() == old


def test_s0_does_not_start_from_its_own_result(tmp_path):
    own = tmp_path / "initial.xyz"
    own.write_text("4\nCoordinates from ORCA-job S0 E -114.290398910078\nC 0 0 0\nO 0 0 1.2\nH 0 1 0\nH 0 -1 0\n")
    isolated = tmp_path / "iso.xyz"
    isolated.write_text("1\nCoordinates from ORCA-job /scratch/.orca_iso_S0_1_2_x/S0 E -1.0\nH 0 0 0\n")
    upstream = tmp_path / "up.xyz"
    upstream.write_text("1\nCoordinates from ORCA-job initial E -56.44\nH 0 0 0\n")

    assert written_by_s0(own) and written_by_s0(isolated)
    assert not written_by_s0(upstream)
    assert not written_by_s0(tmp_path / "missing.xyz")


def test_the_s0_input_of_a_finished_run_is_written_the_same_again(tmp_path, monkeypatch):
    from delfin.config import _load_template_defaults
    from delfin.esd_input_generator import create_state_input

    monkeypatch.chdir(tmp_path)
    (tmp_path / "start.txt").write_text("C 0.0 0.0 -0.529\nO 0.0 0.0 0.677\nH 0.0 0.937 -1.117\nH 0.0 -0.937 -1.117\n")
    esd = tmp_path / "ESD"
    esd.mkdir()
    config = {**_load_template_defaults(), "functional": "PBE0", "main_basisset": "def2-SVP", "PAL": 4,
              "maxcore": 1000, "solvent": "water", "ESD_modul": "yes", "ESD_modus": "TDDFT"}

    first = Path(create_state_input("S0", esd, 0, "water", [], "def2-SVP", "def2-TZVP", config)).read_text()
    # the run is over: S0 optimised, and its geometry copied for the redox steps
    (tmp_path / "initial.xyz").write_text("4\nCoordinates from ORCA-job S0 E -114.29\nC 0 0 -0.52\nO 0 0 0.68\n"
                                          "H 0 0.94 -1.12\nH 0 -0.94 -1.12\n")
    again = Path(create_state_input("S0", esd, 0, "water", [], "def2-SVP", "def2-TZVP", config)).read_text()

    assert again == first
