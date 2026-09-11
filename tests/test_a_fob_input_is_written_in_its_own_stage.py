"""A FoB's input is written into its own stage, whatever the process's working directory.

The FoB input writer found its stage by changing the process's working
directory -- one directory for all threads -- under a lock the stages'
frequency jobs do not take when they change it themselves.  In a dry run of
an archived job of Jerome's the frequency jobs of ox_step_1 and red_step_1
started in the same second as FoB 5 of red_step_3: its input landed in the
job's own folder, and the stage kept the old one.  The writer now reads and
writes the stage's files by path.
"""

from __future__ import annotations

import os

from delfin.config import _load_template_defaults
from delfin.occupier import read_and_modify_file_OCCUPIER


def test_the_input_goes_to_the_stage_when_the_process_is_elsewhere(tmp_path, monkeypatch):
    stage = tmp_path / "red_step_3_OCCUPIER"
    stage.mkdir()
    (stage / "input0.xyz").write_text("4\n\nN 0.0 0.0 0.0\nH 0.0 0.94 0.0\nH 0.81 -0.47 0.0\nH -0.81 -0.47 0.0\n")
    elsewhere = tmp_path / "job_root"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)          # as a frequency job running alongside leaves it
    config = {**_load_template_defaults(), "functional": "PBE0", "main_basisset": "def2-SVP",
              "PAL": 4, "maxcore": 1000, "solvent": "water"}

    read_and_modify_file_OCCUPIER(0, "input5.inp", 0, 1, "water", [], None, "def2-SVP", config, "",
                                  work_dir=stage)

    assert (stage / "input5.inp").read_text().startswith("!")
    assert not (elsewhere / "input5.inp").exists()
    assert os.getcwd() == str(elsewhere)
