"""An OCCUPIER recalc that has to run one FoB again starts it, and compares it with the rest.

A FoB whose output is complete is marked done before scheduling, and its work
-- the only place its energy was recorded -- never runs.  A FoB started from
it waited for that energy: in the archive, Jerome's CoHPP was resubmitted six
times, and every time initial FoB 6 (from FoB 5) and red_step_2 FoB 3 (from
FoB 2) waited until the walltime ended the job.  Measured here on NH3: FoB 5,
whose optimisation had given up, waited on FoB 4 with no ORCA running.

And the stage's comparison was marked done whenever OCCUPIER.txt existed, so
the FoB that did run again was never compared.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin import occupier_flat_extraction as ofe
from delfin import smart_recalc

ENDING = "\nFINAL SINGLE POINT ENERGY      {e}\n\n                             ****ORCA TERMINATED NORMALLY****\n"
SEQUENCE = [
    {"index": 1, "m": 1, "BS": ""},
    {"index": 2, "m": 3, "BS": "", "from": 1},
    {"index": 3, "m": 3, "BS": "3,1", "from": 2},
]


class _Started(Exception):
    pass


def _stage(tmp_path: Path, monkeypatch, *, third_done: bool):
    monkeypatch.setenv("DELFIN_RECALC", "1")
    monkeypatch.setenv("DELFIN_SMART_RECALC", "1")
    folder = tmp_path / "initial_OCCUPIER"
    folder.mkdir()
    for idx, energy in ((1, -56.40), (2, -56.30), (3, -56.35)):
        stem = "input" if idx == 1 else f"input{idx}"
        out = folder / ("output.out" if idx == 1 else f"output{idx}.out")
        (folder / f"{stem}.inp").write_text(f"! PBE0 def2-SVP OPT\n* xyz 0 {SEQUENCE[idx - 1]['m']}\nN 0 0 0\n*\n")
        (folder / f"{stem}.xyz").write_text("1\n\nN 0 0 0\n")
        if idx < 3 or third_done:
            out.write_text(ENDING.format(e=energy))
        else:
            out.write_text("... killed at the walltime\n")
        smart_recalc.store_fingerprint(folder / f"{stem}.inp")
    (folder / "OCCUPIER.txt").write_text("from the first run\n")

    jobs, best_id = ofe._create_occupier_fob_jobs(
        folder_name="initial_OCCUPIER", folder_path=folder, stage_prefix="initial",
        sequence=SEQUENCE, sequence_label="even_seq", total_cores=8,
        global_config={"OCCUPIER_compare": "FSPE"}, ensure_setup=lambda: folder,
        source_folder=None, metals=[], metal_basisset=None, main_basisset="def2-SVP",
        solvent="", occ_results={}, stage_charge=0,
    )
    return folder, {job.job_id: job for job in jobs}, best_id


def test_a_fob_started_from_a_finished_one_does_not_wait_for_it(tmp_path, monkeypatch):
    folder, jobs, _ = _stage(tmp_path, monkeypatch, third_done=False)
    # the scheduler asks every job in the order they were built
    assert jobs["initial_fob_1"].precomplete_check()
    assert jobs["initial_fob_2"].precomplete_check()
    assert not jobs["initial_fob_3"].precomplete_check()

    seen = {}

    def _wait(folder_, source_idx, *, completion_check, **_kw):
        seen["source"], seen["done"] = source_idx, completion_check()
        raise _Started

    monkeypatch.setattr(ofe, "_wait_for_geometry_source", _wait)
    monkeypatch.setattr(ofe, "resolve_sequences_for_delta", lambda *a, **k: {"even_seq": SEQUENCE}, raising=False)
    with pytest.raises(_Started):
        jobs["initial_fob_3"].work(4)

    assert seen == {"source": 2, "done": True}


def test_the_comparison_is_made_again_when_one_fob_ran_again(tmp_path, monkeypatch):
    _, jobs, best_id = _stage(tmp_path, monkeypatch, third_done=False)
    for fob in ("initial_fob_1", "initial_fob_2", "initial_fob_3"):
        jobs[fob].precomplete_check()

    assert not jobs[best_id].precomplete_check()


def test_a_stage_where_every_fob_is_kept_keeps_its_comparison(tmp_path, monkeypatch):
    _, jobs, best_id = _stage(tmp_path, monkeypatch, third_done=True)
    monkeypatch.setattr(ofe, "_update_runtime_cache", lambda *a, **k: None)
    for fob in ("initial_fob_1", "initial_fob_2", "initial_fob_3"):
        assert jobs[fob].precomplete_check()

    assert jobs[best_id].precomplete_check()

