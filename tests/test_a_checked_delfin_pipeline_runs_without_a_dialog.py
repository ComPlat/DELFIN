"""A DELFIN pipeline whose CONTROL.txt checks out runs without a dialog.

Input: a DELFIN job, given as a CONTROL.txt with its geometry (input_file,
the Submit tab's coordinates, or a SMILES). The agent starts one in five
ways: the ops submit_calculation (mode delfin), the ops pipeline_run,
`delfin` in the shell, an sbatch script that runs `delfin`, and the
dashboard's `/submit`.

The same rule as for xtb, CREST and ORCA (user decision 2026-10-22), read
through DELFIN's own functions, without changing the pipeline:
- delfin.config.validate_control_text finds nothing (placeholders,
  required keys, functional, solvent ...);
- for OCCUPIER, delfin.config.occupier_sequence_problems finds nothing;
- the geometry's electrons at `charge` allow `multiplicity_0` (classic
  and manually; OCCUPIER takes its multiplicities from the parity's
  sequence);
- PAL (and PAL x maxcore) fit the host, or the partition and the queue
  cap of 4 when submitted; a job that requests fewer cores than PAL is a
  setup error.
Before, every one of these routes asked, however correct the setup.
"""
from __future__ import annotations

import re
from types import SimpleNamespace as W

import pytest

from delfin.agent import calc_policy as C

WATER = "O 0.0 0.0 0.0\nH 0.0 0.76 0.59\nH 0.0 -0.76 0.59\n"


def _control(pal=48, mult="1", method="classic", solvent="water",
             charge="0"):
    from delfin.define import TEMPLATE
    t = TEMPLATE.replace("charge=[CHARGE]", f"charge={charge}")
    t = t.replace("method=[classic|manually|OCCUPIER]", f"method={method}")
    t = t.replace("[SOLVENT]", solvent)
    t = re.sub(r"^smiles_converter=.*$", "smiles_converter=QUICK", t,
               flags=re.M)
    t = re.sub(r"^multiplicity_0=.*$", f"multiplicity_0={mult}", t,
               flags=re.M)
    t = re.sub(r"^PAL=.*$", f"PAL={pal}", t, flags=re.M)
    return t


@pytest.fixture
def cluster(monkeypatch):
    monkeypatch.setattr(C, "host_kind", lambda env=None: "login_node")
    monkeypatch.setattr(C, "partition_limits", lambda name: {
        "PartitionName": "cpu", "MaxTime": "3-00:00:00",
        "MaxCPUsPerNode": "64"})
    monkeypatch.setattr(C, "queued_jobs", lambda: 0)


def _folder(tmp_path, **kw):
    (tmp_path / "CONTROL.txt").write_text(_control(**kw))
    (tmp_path / "input.txt").write_text(WATER)
    return tmp_path


# -- the check ---------------------------------------------------------------

def test_a_checked_pipeline_is_submitted_unasked(tmp_path, cluster):
    v = C.judge_delfin_submission(_folder(tmp_path), 48, 4500, "12:00:00")
    assert v.calculation and v.admit, v.reason


def test_delfin_s_own_control_check_decides(tmp_path, cluster):
    v = C.judge_delfin_submission(_folder(tmp_path, solvent="[SOLVENT]"),
                                  48, 4500, "12:00:00")
    assert not v.admit and "SOLVENT" in v.reason


def test_a_multiplicity_the_electrons_cannot_have_is_refused(tmp_path,
                                                            cluster):
    v = C.judge_delfin_submission(_folder(tmp_path, mult="2"), 48, 4500,
                                  "12:00:00")
    assert not v.admit and "multiplicity 2" in v.reason


def test_occupier_takes_its_multiplicity_from_the_sequence(tmp_path,
                                                           cluster):
    v = C.judge_delfin_submission(
        _folder(tmp_path, method="OCCUPIER", mult=""), 48, 4500, "12:00:00")
    assert v.admit, v.reason


def test_a_job_smaller_than_its_pal_is_a_setup_error(tmp_path, cluster):
    v = C.judge_delfin_submission(_folder(tmp_path), 12, 4500, "12:00:00")
    assert not v.admit and "asks for 48 cores and the job requests 12" \
        in v.reason


def test_a_smiles_is_counted_like_a_geometry(tmp_path, cluster):
    v = C.judge_delfin_job(_control(mult="2"), 48, 4500, "12:00:00",
                           base_dir=tmp_path, name="CONTROL.txt",
                           geometry_text="O")
    assert not v.admit and "multiplicity 2" in v.reason
    v = C.judge_delfin_job(_control(), 48, 4500, "12:00:00",
                           base_dir=tmp_path, name="CONTROL.txt",
                           geometry_text="O")
    assert v.admit, v.reason


# -- the shell and sbatch ------------------------------------------------------

def _workstation(monkeypatch, cpus=8, memory_mb=32000):
    # The host's cores and memory are pinned: the template's maxcore of
    # 4500 MB at PAL 2 is 9000 MB, more than half of a 16 GB machine.
    monkeypatch.setattr(C, "host_kind", lambda env=None: "workstation")
    monkeypatch.setattr(C, "_available_cpus", lambda: cpus)
    monkeypatch.setattr(C, "_memory_mb", lambda: memory_mb)


def test_a_local_run_keeps_to_the_host(tmp_path, monkeypatch):
    _workstation(monkeypatch)
    _folder(tmp_path, pal=48)
    v = C.judge_command(["delfin"], tmp_path, lambda s: False,
                        tool_wall_s=600, kind="workstation")
    assert v.calculation and not v.admit and "48" in v.reason
    _folder(tmp_path, pal=2)
    v = C.judge_command(["delfin"], tmp_path, lambda s: False,
                        tool_wall_s=600, kind="workstation")
    assert v.admit, v.reason


def test_a_local_run_keeps_to_half_the_memory(tmp_path, monkeypatch):
    _workstation(monkeypatch, memory_mb=16000)
    _folder(tmp_path, pal=2)
    v = C.judge_command(["delfin"], tmp_path, lambda s: False,
                        tool_wall_s=600, kind="workstation")
    assert not v.admit and v.reason.startswith(
        "delfin (CONTROL.txt) asks for 9000 MB"), v.reason


def test_only_a_run_is_a_calculation(tmp_path):
    _folder(tmp_path)
    for cmd in ("delfin --purge", "delfin cleanup", "delfin --define"):
        v = C.judge_command([cmd], tmp_path, lambda s: False,
                            kind="workstation")
        assert not v.calculation, cmd


def test_a_job_script_that_runs_delfin_is_checked(tmp_path, cluster):
    _folder(tmp_path)
    (tmp_path / "job.sh").write_text(
        "#!/bin/bash\n#SBATCH --ntasks=48\n#SBATCH --time=12:00:00\n"
        "cd $SLURM_SUBMIT_DIR\ndelfin\n")
    v = C.judge_command(["sbatch job.sh"], tmp_path, lambda s: False,
                        kind="login_node")
    assert v.admit, v.reason
    (tmp_path / "job.sh").write_text(
        "#!/bin/bash\n#SBATCH --ntasks=8\ndelfin\n")
    v = C.judge_command(["sbatch job.sh"], tmp_path, lambda s: False,
                        kind="login_node")
    assert not v.admit and "48" in v.reason


# -- the dashboard's /submit ---------------------------------------------------

def _refs(control, coords, batch=""):
    return {"control_widget": W(value=control),
            "coords_widget": W(value=coords),
            "smiles_batch_widget": W(value=batch),
            "job_type_widget": W(value="12:00:00"),
            "custom_time_widget": W(value=12)}


def test_the_submit_tab_sends_a_checked_job_unasked(cluster):
    from delfin.dashboard.tab_agent import _submit_tab_verdict
    v = _submit_tab_verdict(_refs(_control(), WATER))
    assert v is not None and v.admit, v and v.reason


def test_the_submit_tab_asks_for_a_wrong_spin(cluster):
    from delfin.dashboard.tab_agent import _submit_tab_verdict
    v = _submit_tab_verdict(_refs(_control(mult="2"), WATER))
    assert v is not None and not v.admit and "multiplicity 2" in v.reason


def test_a_batch_keeps_its_dialog(cluster):
    from delfin.dashboard.tab_agent import _submit_tab_verdict
    assert _submit_tab_verdict(_refs(_control(), WATER,
                                     batch="a, O\nb, CCO\n")) is None


def test_the_shell_gate_starts_a_checked_run(tmp_path, monkeypatch):
    """Through the bash gate _dispatch uses: `delfin` was an acting
    command, asked every time; checked, it starts, and other forms keep
    their own decision."""
    from delfin.agent import api_client as A
    _workstation(monkeypatch)
    _folder(tmp_path, pal=2)

    def gate(cmd):
        perms = A.KitToolPermissions(workspace=tmp_path, mode="default",
                                     agent_role="solo_agent")
        perms.confirm_callback = None
        return A._doc_executor._gate_bash_call("bash", {"command": cmd},
                                               perms, [], [], [])
    assert gate("delfin") is None
    assert gate("delfin --purge") is not None
    _folder(tmp_path, pal=48)
    out = gate("delfin")
    assert out is not None and "delfin (CONTROL.txt) asks for 48 cores" in out
