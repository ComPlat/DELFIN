"""calc_policy reads charge, spin, cores and limits the way the programs do.

Input: delfin/agent/calc_policy.py on its own, without the gate. Each case
is a spelling the programs accept and the agent writes: xtb's ``-P`` and
``.UHF``, ORCA's ``* xyz``/``* xyzfile`` line, ``%pal``/``PALn``, Slurm's
six time formats, ``#SBATCH`` directives overridden on the command line.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import calc_policy as C

WATER = "3\n\nO 0 0 0\nH 0 0.76 0.59\nH 0 -0.76 0.59\n"


def test_the_first_atom_after_the_orca_geometry_line_is_counted(tmp_path):
    # A regex that let `\s*` cross the line break read the first atom as a
    # file name and counted 28 electrons for butane's 34.
    inp = tmp_path / "w.inp"
    inp.write_text("! HF def2-SVP\n* xyz 0 1\nO 0 0 0\nH 0 0.76 0.59\n"
                   "H 0 -0.76 0.59\n*\n")
    assert C.check_orca_input(inp)[2] == ""
    inp.write_text("! HF def2-SVP\n* xyz 0 2\nO 0 0 0\nH 0 0.76 0.59\n"
                   "H 0 -0.76 0.59\n*\n")
    assert "multiplicity 2" in C.check_orca_input(inp)[2]


def test_an_xyzfile_reference_is_read(tmp_path):
    (tmp_path / "w.xyz").write_text(WATER)
    inp = tmp_path / "w.inp"
    inp.write_text("! HF def2-SVP\n* xyzfile 1 2 w.xyz\n")
    assert C.check_orca_input(inp)[2] == ""      # 9 electrons, doublet


def test_cores_and_memory_of_an_orca_input():
    assert C.orca_resources("! PAL8 r2SCAN-3c\n%maxcore 3000\n") == (8, 3000)
    assert C.orca_resources("%pal\n nprocs 12\nend\n") == (12, None)


def test_an_unknown_element_is_not_a_pass(tmp_path):
    (tmp_path / "x.xyz").write_text("1\n\nXq 0 0 0\n")
    v = C.judge_command(["xtb x.xyz"], tmp_path, lambda s: False,
                        tool_wall_s=120, kind="workstation")
    assert v.calculation and not v.admit and "unknown element" in v.reason


def test_xtb_reads_uhf_from_the_file_it_reads(tmp_path, monkeypatch):
    monkeypatch.setattr(C, "_available_cpus", lambda: 8)
    (tmp_path / "w.xyz").write_text(WATER)
    (tmp_path / ".CHRG").write_text("1\n")
    (tmp_path / ".UHF").write_text("1\n")
    v = C.judge_command(["xtb w.xyz --opt"], tmp_path, lambda s: False,
                        tool_wall_s=120, kind="workstation")
    assert v.admit, v.reason


def test_a_directory_made_in_the_same_command_is_resolved(tmp_path,
                                                          monkeypatch):
    monkeypatch.setattr(C, "_available_cpus", lambda: 8)
    (tmp_path / "w.xyz").write_text(WATER)
    v = C.judge_command(["mkdir run", "cd run", "xtb ../w.xyz --opt"],
                        tmp_path, lambda s: s.startswith("mkdir"),
                        tool_wall_s=120, kind="workstation")
    assert v.admit, v.reason


@pytest.mark.parametrize("text,seconds", [
    ("30", 1800), ("30:15", 1815), ("02:00:00", 7200), ("1-12", 129600),
    ("1-00:30", 88200), ("3-00:00:00", 259200), ("UNLIMITED", None),
    ("nonsense", None)])
def test_slurm_time_formats(text, seconds):
    assert C.slurm_seconds(text) == seconds


def test_the_command_line_overrides_the_script(tmp_path):
    opts = C.sbatch_options("#!/bin/bash\n#SBATCH -t 01:00:00\n#SBATCH -n 4\n",
                            ["--time=03:00:00"])
    assert opts["--time"] == "03:00:00" and opts["--ntasks"] == "4"


def test_a_wrapped_and_redirected_command_is_taken_apart():
    run = C._parse_run("OMP_NUM_THREADS=2 timeout 30m nice -n 5 xtb a.xyz "
                       "--opt > out.log 2>&1")
    assert run.program == "xtb" and run.wall_s == 1800
    assert run.env == {"OMP_NUM_THREADS": "2"}
    assert run.argv == ["a.xyz", "--opt"]


def test_a_command_that_starts_nothing_is_no_calculation(tmp_path):
    v = C.judge_command(["ls -la", "grep x y"], tmp_path, lambda s: True,
                        kind="workstation")
    assert not v.calculation


def test_the_host_is_read_from_the_environment(monkeypatch):
    assert C.host_kind({"SLURM_JOB_ID": "42"}) == "slurm_job"
    monkeypatch.setattr(C.shutil, "which", lambda name: None)
    assert C.host_kind({}) == "workstation"
    monkeypatch.setattr(C.shutil, "which", lambda name: "/usr/bin/" + name)
    assert C.host_kind({}) == "login_node"


def test_paths_are_lexical(tmp_path):
    assert C._resolve(tmp_path / "run", "../w.xyz") == Path(
        str(tmp_path)) / "w.xyz"


def test_the_queue_is_asked_through_the_scheduler_client(monkeypatch):
    """One throttled, shared address for scheduler queries
    (delfin.scheduler_client); while its gate is closed, the last answer
    to the same question is used, and without one the count is unknown."""
    from delfin import scheduler_client as sc
    client = sc.SchedulerClient(min_interval_s=60, clock=lambda: 0.0,
                                runner=lambda argv: "101\n102\n")
    monkeypatch.setattr(sc, "default_client", client)
    assert C.queued_jobs() == 2
    assert client.queries_run == 1
    assert C.queued_jobs() == 2 and client.queries_run == 1   # cached
    fresh = sc.SchedulerClient(min_interval_s=60, clock=lambda: 0.0,
                               runner=lambda argv: "x\n")
    fresh._kind_last_run["squeue"] = 0.0      # another squeue just ran
    monkeypatch.setattr(sc, "default_client", fresh)
    assert C.queued_jobs() is None and fresh.queries_run == 0
