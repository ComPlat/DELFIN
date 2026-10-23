"""A calculation whose setup checks out and whose resources fit runs unasked.

Input: the agent's bash gate and its MCP gate, for commands and tool calls
that start xtb, CREST, ORCA or an sbatch job (delfin/agent/calc_policy.py).

User decision 2026-10-22: a calculation starts without a dialog when its
setup passes the check (input present, electron count consistent with
charge and spin, ORCA input free of validator errors) and its resources
are moderate and compliant. On a login node that means at most two cores
and a 30 min wall-time limit locally, and anything bigger as an sbatch job
within the partition's limits, with at most four of the user's jobs
queued. A calculation that fails the check is asked about with the reason,
and refused with the reason when nobody can be asked.

Before: in a DELFIN source tree `xtb <file>.xyz` with any flags ran unasked
(`-P 64` on a login node, an odd electron count); everywhere else every xtb
run and every sbatch asked, however small. The host is chosen per test, so
this file means the same on every machine.
"""
from __future__ import annotations

import pytest

from delfin.agent import api_client as A
from delfin.agent import calc_policy as C

BUTANE = """14

C   -1.9   0.0   0.1
C   -0.6   0.7  -0.2
C    0.6  -0.2   0.1
C    1.9   0.5  -0.2
H   -2.7   0.6  -0.2
H   -2.0  -0.9  -0.5
H   -2.0  -0.2   1.2
H   -0.6   1.6   0.4
H   -0.6   1.0  -1.3
H    0.6  -1.1  -0.5
H    0.6  -0.5   1.2
H    2.7  -0.1   0.1
H    2.0   1.4   0.4
H    2.0   0.7  -1.3
"""


@pytest.fixture
def ws(tmp_path, monkeypatch):
    (tmp_path / "b.xyz").write_text(BUTANE)
    atoms = "\n".join(BUTANE.splitlines()[2:])
    (tmp_path / "b.inp").write_text(
        f"! r2SCAN-3c Opt\n%pal nprocs 16 end\n%maxcore 2000\n"
        f"* xyz 0 1\n{atoms}\n*\n")
    (tmp_path / "doublet.inp").write_text(
        f"! r2SCAN-3c Opt\n* xyz 0 2\n{atoms}\n*\n")
    (tmp_path / "job.sh").write_text(
        "#!/bin/bash\n#SBATCH --ntasks=16\n#SBATCH --time=02:00:00\n"
        "cd $SLURM_SUBMIT_DIR\n$ORCA_DIR/orca b.inp > b.out\n")
    monkeypatch.setattr(C, "_available_cpus", lambda: 8)
    monkeypatch.setattr(C, "partition_limits", lambda name: {
        "PartitionName": name or "cpu", "MaxTime": "3-00:00:00",
        "MaxNodes": "4", "MaxCPUsPerNode": "64"})
    monkeypatch.setattr(C, "queued_jobs", lambda: 1)
    return tmp_path


def _host(monkeypatch, kind):
    monkeypatch.setattr(C, "host_kind", lambda env=None: kind)


def _perms(ws, mode="default", confirm=None):
    perms = A.KitToolPermissions(workspace=ws, mode=mode,
                                 agent_role="solo_agent")
    perms.confirm_callback = confirm
    return perms


def _bash(ws, cmd, mode="default", confirm=None, background=False,
          timeout_s=None):
    args = {"command": cmd}
    if timeout_s is not None:
        args["timeout_s"] = timeout_s
    name = "bash_background" if background else "bash"
    # The function _dispatch calls for both tools, so the wall-time limit
    # the gate sees is the one the run will have.
    return A._doc_executor._gate_bash_call(name, args, _perms(ws, mode,
                                                              confirm), [])


# -- runs unasked ----------------------------------------------------------

def test_a_checked_small_xtb_run_starts_unasked(ws, monkeypatch):
    _host(monkeypatch, "workstation")
    assert _bash(ws, "xtb b.xyz --opt --chrg 0 --uhf 0 -P 2 > xtb.out") is None


def test_on_a_login_node_a_two_core_run_with_a_foreground_limit_starts(
        ws, monkeypatch):
    _host(monkeypatch, "login_node")
    assert _bash(ws, "xtb b.xyz --opt -P 2") is None


def test_on_a_login_node_a_background_run_needs_its_own_limit(ws, monkeypatch):
    _host(monkeypatch, "login_node")
    assert _bash(ws, "crest b.xyz -T 2", background=True,
                 timeout_s=1800) is None
    assert _bash(ws, "timeout 1800 crest b.xyz -T 2",
                 background=True) is None
    out = _bash(ws, "crest b.xyz -T 2", background=True)
    assert out is not None and "30 min" in out


def test_a_job_within_the_partition_is_submitted_unasked(ws, monkeypatch):
    _host(monkeypatch, "login_node")
    assert _bash(ws, "sbatch job.sh") is None


# -- not unasked, and the reason is given ----------------------------------

def test_a_spin_the_electrons_cannot_have_is_refused_with_the_reason(
        ws, monkeypatch):
    _host(monkeypatch, "workstation")
    out = _bash(ws, "xtb b.xyz --opt --uhf 1")
    assert out is not None and "unpaired" in out


def test_an_orca_input_with_the_wrong_multiplicity_is_refused(ws, monkeypatch):
    _host(monkeypatch, "workstation")
    out = _bash(ws, "orca doublet.inp > d.out")
    assert out is not None and "multiplicity 2" in out


def test_a_wide_run_on_a_login_node_is_sent_to_sbatch(ws, monkeypatch):
    _host(monkeypatch, "login_node")
    out = _bash(ws, "xtb b.xyz --opt -P 16")
    assert out is not None and "sbatch" in out


def test_a_fifth_job_is_asked_about(ws, monkeypatch):
    _host(monkeypatch, "login_node")
    monkeypatch.setattr(C, "queued_jobs", lambda: 4)
    out = _bash(ws, "sbatch job.sh")
    assert out is not None and "at most 4" in out


def test_a_job_beyond_the_partition_time_is_asked_about(ws, monkeypatch):
    _host(monkeypatch, "login_node")
    out = _bash(ws, "sbatch --time=5-00:00:00 job.sh")
    assert out is not None and "partition limit" in out


def test_the_dialog_shows_why_it_asks(ws, monkeypatch):
    _host(monkeypatch, "login_node")
    shown = []

    def confirm(name, args, preview):
        shown.append(preview)
        return False
    _bash(ws, "xtb b.xyz --opt -P 16", confirm=confirm)
    assert shown and "sbatch" in shown[0]


def test_plan_mode_starts_no_calculation(ws, monkeypatch):
    _host(monkeypatch, "workstation")
    assert _bash(ws, "xtb b.xyz --opt -P 2", mode="plan") is not None


def test_a_chain_with_another_command_is_not_admitted_by_its_calculation(
        ws, monkeypatch):
    _host(monkeypatch, "workstation")
    assert _bash(ws, "xtb b.xyz --opt -P 2 && rm -rf /tmp/x") is not None


def test_a_delfin_tree_no_longer_runs_any_xtb_unasked():
    from pathlib import Path
    import delfin
    tree = Path(delfin.__file__).resolve().parents[1]
    perms = A.KitToolPermissions(workspace=tree, mode="default")
    assert not perms.matches_bash_auto_allow("xtb mol.xyz --opt -P 64")
    assert perms.matches_bash_auto_allow("xtb --version")


# -- the ops tools ---------------------------------------------------------

def _mcp(ws, name, args, mode="default"):
    return A._doc_executor._gate_mcp_tool(f"mcp__delfin-ops__{name}", args,
                                          _perms(ws, mode))


def test_a_checked_orca_submission_is_sent_unasked(ws, monkeypatch):
    _host(monkeypatch, "login_node")
    sub = ws / "sub"
    sub.mkdir()
    (sub / "b.inp").write_text((ws / "b.inp").read_text())
    assert _mcp(ws, "submit_calculation", {
        "folder": str(sub), "mode": "orca", "pal": 16, "maxcore": 2000,
        "time_limit": "12:00:00"}) is None


def test_a_submission_of_a_wrong_input_is_refused_with_the_reason(
        ws, monkeypatch):
    _host(monkeypatch, "login_node")
    sub = ws / "sub"
    sub.mkdir()
    (sub / "d.inp").write_text((ws / "doublet.inp").read_text())
    out = _mcp(ws, "submit_calculation", {
        "folder": str(sub), "mode": "orca", "pal": 16, "maxcore": 2000,
        "time_limit": "12:00:00"})
    assert out is not None and "multiplicity" in out


def test_a_pipeline_submission_still_asks(ws, monkeypatch):
    _host(monkeypatch, "login_node")
    assert _mcp(ws, "submit_calculation", {
        "folder": str(ws), "mode": "delfin"}) is not None


def test_a_local_orca_run_on_a_login_node_is_not_unasked(ws, monkeypatch):
    # run_orca_input has no wall-time limit of its own.
    _host(monkeypatch, "login_node")
    assert _mcp(ws, "run_orca_input",
                {"input_file": str(ws / "b.inp")}) is not None


# -- what the run reports ---------------------------------------------------

def test_the_run_names_what_was_started_where(tmp_path, monkeypatch):
    """Through the real executor: the result carries the check's line,
    whether or not xtb is installed on this machine (an absent program
    still returns a JSON result with its exit code)."""
    import json
    _host(monkeypatch, "workstation")
    monkeypatch.setattr(C, "_available_cpus", lambda: 8)
    (tmp_path / "h2.xyz").write_text("2\n\nH 0 0 0\nH 0 0 0.74\n")
    out = json.loads(A._doc_executor.execute(
        "bash", {"command": "xtb h2.xyz --sp -P 1 > h2.out 2>&1"},
        _perms(tmp_path)))
    assert "calculation_check" in out, out
    assert "xtb h2.xyz, 1 core(s) on this machine" in out["calculation_check"]
