"""Red tests for run_evidence.collect (CB0A phase 1).

Drives the public path: delfin.agent.run_evidence.collect builds an
Evidence from files a real DELFIN run leaves behind. Fixtures are
constructed by hand in tmp_path — no calculation is started.

These tests must run against the UNFIXED code and fail with ImportError
until run_evidence.py exists (that is the red state).
"""

from __future__ import annotations

import pytest


def make_workspace(tmp_path):
    """A minimal pre-run DELFIN workspace: CONTROL.txt only."""
    ws = tmp_path / "run"
    ws.mkdir()
    (ws / "CONTROL.txt").write_text(
        "charge=0\nmultiplicity=1\n"
        "smiles=O\n"
        "PAL=4\nmaxcore=2000\n"
        "method=B3LYP\nbasis=def2-SVP\n"
        "auto=true\n"
    )
    return ws


def test_collect_reads_control_defaults_from_a_pre_run_workspace(tmp_path):
    from delfin.agent.run_evidence import collect

    ws = make_workspace(tmp_path)
    ev = collect(ws)
    # parse_control_text keeps CONTROL values as their literal strings;
    # evidence exposes the numeric view for comparison keys.
    assert ev.control["charge"] == 0
    assert ev.control["PAL"] == 4
    assert ev.steps == []  # no inputs yet: nothing claimed about steps


def test_collect_never_writes(tmp_path):
    from delfin.agent.run_evidence import collect

    ws = make_workspace(tmp_path)
    before = sorted(str(p) for p in ws.rglob("*"))
    collect(ws)
    after = sorted(str(p) for p in ws.rglob("*"))
    assert before == after


def test_collect_collects_orca_inp_per_step_with_charge_and_method(tmp_path):
    from delfin.agent.run_evidence import collect

    ws = make_workspace(tmp_path)
    stepdir = ws / "calc" / "classic_opt_freq_xtb"
    stepdir.mkdir(parents=True)
    (stepdir / "dft.inp").write_text(
        "! B3LYP def2-SVP OPT FREQ\n"
        "%cpcm smd true\n"
        "* xyz 0 1\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH -0.24 0.93 0.0\n*\n"
    )
    (stepdir / "dft.out").write_text(
        "ORCA TERMINATED NORMALLY\n"
    )
    ev = collect(ws)
    assert len(ev.steps) == 1
    step = ev.steps[0]
    assert step["charge"] == 0
    assert step["multiplicity"] == 1
    assert step["program"] == "orca"
    assert step["terminated_normally"] is True
    assert step["functional"] == "B3LYP"
    assert step["basis"] == "def2-SVP"


def test_collect_flags_output_terminated_abnormally(tmp_path):
    from delfin.agent.run_evidence import collect

    ws = make_workspace(tmp_path)
    stepdir = ws / "calc" / "step1"
    stepdir.mkdir(parents=True)
    (stepdir / "dft.inp").write_text(
        "* xyz 0 1\nO 0.0 0.0 0.0\nH 0.9 0.0 0.0\n*\n"
    )
    (stepdir / "dft.out").write_text(
        "SCF NOT CONVERGED AFTER 150 STEPS\nslurm deadline: walltime exhausted\n"
    )
    ev = collect(ws)
    assert ev.steps[0]["terminated_normally"] is False


def test_collect_reads_slurm_records_from_submit_log(tmp_path):
    from delfin.agent.run_evidence import collect

    ws = make_workspace(tmp_path)
    (ws / "submit.log").write_text(
        "Submitted batch job 12345678\n"
        "SBATCH --time=02:00:00 --cpus-per-task=8 --mem=64000MB\n"
    )
    ev = collect(ws)
    rec = ev.slurm_records[0]
    assert rec["job_id"] == "12345678"
    assert rec["requested_cores"] == 8
    assert rec["requested_time"] == "02:00:00"


def test_collect_ignores_agent_prose_files(tmp_path):
    from delfin.agent.run_evidence import collect

    ws = make_workspace(tmp_path)
    (ws / "notes.txt").write_text("charge=+7 multiplicity=99")
    ev = collect(ws)
    assert ev.control["charge"] == 0  # from CONTROL.txt, never the notes
