"""sbatch_command: a broken partition probe stays visible but never breaks.

Two directions (wave-17 lesson 5): a read error must not stop the
command from being built (the plain ``sbatch script`` form is the
documented fallback), and the failure must be LOGGED rather than
swallowed by a bare ``except Exception: pass``.
"""
from __future__ import annotations

import logging

from delfin.slurm_submit import sbatch_command


def test_read_error_does_not_break_the_command(tmp_path, caplog):
    missing = tmp_path / "not_here.sbatch"
    with caplog.at_level(logging.WARNING, logger="delfin.slurm_submit"):
        cmd = sbatch_command("sbatch", missing, tmp_path)
    assert cmd == ["sbatch", str(missing)]


def test_read_error_is_logged_not_swallowed(tmp_path, caplog):
    missing = tmp_path / "not_here.sbatch"
    with caplog.at_level(logging.WARNING, logger="delfin.slurm_submit"):
        sbatch_command("sbatch", missing, tmp_path)
    assert any("partition choice skipped" in r.message for r in caplog.records)


def test_working_script_still_carries_the_partition_choice(tmp_path, monkeypatch):
    script = tmp_path / "job.sbatch"
    script.write_text("#!/bin/bash\nsleep 1\n")
    # One configured partition only (env route, read in
    # configured_partitions): choose_partitions returns it without probing
    # the scheduler, so this runs offline.
    monkeypatch.setenv("DELFIN_SLURM_PARTITIONS", "cpu")
    cmd = sbatch_command("sbatch", script, tmp_path)
    assert "--partition=cpu" in cmd
    assert cmd[-1] == str(script)
