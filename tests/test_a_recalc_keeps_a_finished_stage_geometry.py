"""A recalc does not put OCCUPIER's start geometry over a stage's finished result.

OCCUPIER hands the configuration it chose to ``<stage>.xyz``, where the
stage's frequency job starts from it and writes its optimised geometry over
it.  A recalc repeats the hand-over and keeps the finished job, so the start
geometry came back over the result -- including the geometry IMAG had freed
of an imaginary mode.  In the archive, 293 OCCUPIER stage geometries hold a
FoB's geometry although the stage's job had finished; 291 are in folders that
were run more than once.
"""

from __future__ import annotations

from pathlib import Path

from delfin import occupier_flat_extraction as ofe

ENDING = "\nFINAL SINGLE POINT ENERGY      {e}\n\n                             ****ORCA TERMINATED NORMALLY****\n"


def _handover(tmp_path: Path, monkeypatch, *, stage_done: bool, header: str, recalc: bool = True):
    monkeypatch.setenv("DELFIN_RECALC", "1" if recalc else "0")
    folder = tmp_path / "ox_step_1_OCCUPIER"
    folder.mkdir()
    (folder / "input2.xyz").write_text("1\nCoordinates from ORCA-job input2 E -56.2\nN 0 0 0.1\n")
    (tmp_path / "ox_step_1.xyz").write_text(f"1\n{header}\nN 0 0 0.2\n")
    (tmp_path / "ox_step_1.out").write_text(ENDING.format(e=-56.21) if stage_done else "... killed\n")
    monkeypatch.setattr(ofe, "read_occupier_file", lambda *a, **k: (2, "", 2, None))
    ofe._update_runtime_cache("ox_step_1_OCCUPIER", folder, {}, {})
    return (tmp_path / "ox_step_1.xyz").read_text()


def test_a_recalc_does_not_put_the_start_geometry_over_the_finished_result(tmp_path, monkeypatch):
    kept = _handover(tmp_path, monkeypatch, stage_done=True,
                     header="Coordinates from ORCA-job /scratch/.orca_iso_ox_step_1_7_8_x/ox_step_1 E -56.21")

    assert "N 0 0 0.2" in kept


def test_an_unfinished_stage_job_still_starts_from_occupiers_choice(tmp_path, monkeypatch):
    handed = _handover(tmp_path, monkeypatch, stage_done=False,
                       header="Coordinates from ORCA-job ox_step_1 E -56.19")

    assert "N 0 0 0.1" in handed


def test_a_stage_geometry_from_anywhere_else_is_replaced_as_before(tmp_path, monkeypatch):
    handed = _handover(tmp_path, monkeypatch, stage_done=True,
                       header="Coordinates from ORCA-job input3 E -56.20")

    assert "N 0 0 0.1" in handed


def test_a_full_run_hands_over_what_occupier_chose_this_time(tmp_path, monkeypatch):
    handed = _handover(tmp_path, monkeypatch, stage_done=True, recalc=False,
                       header="Coordinates from ORCA-job ox_step_1 E -56.21")

    assert "N 0 0 0.1" in handed
