"""Adversarial tests for SL1 partition choice and time parsing.

Reviewer nacht-s14, package SL1 (wave 25). These tests attack the
UNCHANGED owned code; a red test is a FINDING for the builder
(nacht-s13), a green test documents defensible behaviour.

No scheduler is called: sbatch/scontrol are injectable via the
``sbatch=`` / ``scontrol_output=`` parameters of the owned modules, or
monkeypatched where no seam exists.
"""
import subprocess

import pytest

from delfin import slurm_submit
from delfin.slurm_submit import (
    choose_partitions,
    partition_accepts,
    partition_list,
    sbatch_command,
    time_limit_seconds,
)


# ------------------------------------------------------------------ #
# partition_list: dedup, whitespace, junk
# ------------------------------------------------------------------ #
def test_partition_list_splits_and_dedups():
    assert partition_list("cpu, cpu_il,cpu") == ("cpu", "cpu_il")


def test_partition_list_empty_and_none():
    assert partition_list("") == ()
    assert partition_list(None) == ()
    assert partition_list([]) == ()


# ------------------------------------------------------------------ #
# partition_accepts: the answer SLURM gave decides, not our guess
# ------------------------------------------------------------------ #
def test_partition_accepts_true_on_test_only_ok(tmp_path, monkeypatch):
    script = tmp_path / "s.sh"
    script.write_text("#!/bin/bash\n", encoding="utf-8")

    captured = {}

    def fake_run(cmd, **kwargs):
        captured["cmd"] = cmd
        return subprocess.CompletedProcess(cmd, 0, stdout="job to start at 12:00",
                                           stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)
    assert partition_accepts("cpu", tmp_path, script) is True
    # --test-only must be in the command so NOTHING is submitted
    assert "--test-only" in captured["cmd"]
    assert all(not c.startswith("--parsable") for c in captured["cmd"])


def test_partition_accepts_false_on_allocation_failure(tmp_path, monkeypatch):
    script = tmp_path / "s.sh"
    script.write_text("#!/bin/bash\n", encoding="utf-8")

    def fake_run(cmd, **kwargs):
        return subprocess.CompletedProcess(
            cmd, 1, stdout="sbatch: error: CPU count invalid", stderr="allocation failure")

    monkeypatch.setattr(subprocess, "run", fake_run)
    assert partition_accepts("cpu", tmp_path, script) is False


def test_partition_accepts_none_when_sbatch_missing(tmp_path, monkeypatch):
    script = tmp_path / "s.sh"
    script.write_text("#!/bin/bash\n", encoding="utf-8")

    def fake_run(cmd, **kwargs):
        raise FileNotFoundError("no sbatch here")

    monkeypatch.setattr(subprocess, "run", fake_run)
    assert partition_accepts("cpu", tmp_path, script) is None


# ------------------------------------------------------------------ #
# choose_partitions: unknown-availability candidates must NOT be listed
# ------------------------------------------------------------------ #
def test_choose_partitions_drops_unprobed_candidate(tmp_path, monkeypatch):
    """When sbatch cannot be asked at all (None), the candidate must be
    dropped -- listing a partition SLURM may refuse makes sbatch reject
    the WHOLE job."""
    script = tmp_path / "s.sh"
    script.write_text("#!/bin/bash\n", encoding="utf-8")

    def fake_run(cmd, **kwargs):
        raise FileNotFoundError("no sbatch")

    monkeypatch.setattr(subprocess, "run", fake_run)
    chosen = choose_partitions(script, tmp_path, partitions=("cpu", "cpu_il"))
    assert chosen == ""


def test_choose_partitions_keeps_only_accepted(tmp_path, monkeypatch):
    script = tmp_path / "s.sh"
    script.write_text("#!/bin/bash\n", encoding="utf-8")

    def fake_run(cmd, **kwargs):
        name = [c for c in cmd if c.startswith("--partition=")][0]
        if name == "--partition=cpu":
            return subprocess.CompletedProcess(cmd, 0, "job to start at now", "")
        return subprocess.CompletedProcess(cmd, 1, "", "allocation failure")

    monkeypatch.setattr(subprocess, "run", fake_run)
    chosen = choose_partitions(script, tmp_path, partitions=("cpu", "cpu_il"))
    assert chosen == "cpu"


# ------------------------------------------------------------------ #
# sbatch_command: silent except (finding candidate)
# ------------------------------------------------------------------ #
def test_sbatch_command_still_appends_script_when_probe_broken(tmp_path, monkeypatch):
    script = tmp_path / "s.sh"
    script.write_text("#!/bin/bash\n", encoding="utf-8")
    monkeypatch.setattr(slurm_submit, "configured_partitions", lambda *a, **k: ("cpu", "cpu_il"))
    # corrupt read: partition_accepts will raise via monkeypatched subprocess
    def fake_run(cmd, **kwargs):
        raise RuntimeError("boom")

    monkeypatch.setattr(subprocess, "run", fake_run)
    cmd = sbatch_command("sbatch", script, tmp_path)
    assert str(script) in cmd
    # A job whose availability could not be probed must not carry a
    # partition list that only partial probes produced.
    assert "--partition=" not in " ".join(cmd)


def test_sbatch_command_script_that_names_partition_is_left_alone(tmp_path):
    script = tmp_path / "s.sh"
    script.write_text("#SBATCH --partition=gpu_h100\n#!/bin/bash\n", encoding="utf-8")
    cmd = sbatch_command("sbatch", script, tmp_path)
    assert "--partition=" not in cmd
    assert cmd[-1] == str(script)


# ------------------------------------------------------------------ #
# time_limit_seconds: SLURM formats
# ------------------------------------------------------------------ #
def test_time_limit_seconds_formats():
    assert time_limit_seconds("02:00:00") == 7200
    assert time_limit_seconds("90") == 5400          # bare = minutes
    assert time_limit_seconds("30:00") == 1800       # MM:SS
    assert time_limit_seconds("1-12:00:00") == 129600
    assert time_limit_seconds("1-12") == 129600      # D-HH pads


def test_time_limit_seconds_broken_input_raises():
    for bad in ("", "abc"):
        with pytest.raises((ValueError, IndexError)):
            time_limit_seconds(bad)


def test_time_limit_seconds_extra_field_is_not_silently_minutes():
    """'1:2:3:4' is not a SLURM time: it must be refused with a reason.

    REWRITTEN along its intent when the defect was fixed (SL1 phase 1,
    rule 6). The old version pinned the misbehaviour on the unchanged
    code: a 4-field colon form fell through to the bare-minutes branch
    and became 60 s (slurm_submit.py:190-195, reviewer F1). The new
    decision: >3 colon fields raise ValueError.
    """
    with pytest.raises(ValueError):
        time_limit_seconds("1:2:3:4")
