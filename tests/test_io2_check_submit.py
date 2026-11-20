"""IO2 phase 1: the shared submit-check helper.

Target: ``hpc_io_check.check_submit(script_path_or_text, env=None) ->
(refused, warned)`` returns the check's findings split by severity:
``refused`` a list of formatted finding strings, ``warned`` likewise.
One helper both submit routes (dashboard backend, agent bash gate) and
``cluster_tools.submit`` call, so their verdicts cannot drift.

Red control: before the helper exists these tests fail on import /
AttributeError. The wiring tests (test_io2_backend_wiring.py) fail
because backend_slurm never calls it.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from delfin import hpc_io_check


LUSTRE_MOUNTS = "mgmt@o2ib:fs /pfs lustre rw 0 0\n"

WALK_SCRIPT = "#!/bin/bash\ndu -sh /pfs/data6/workspace/\n"
SCRATCH_SCRIPT = (
    "#!/bin/bash\n"
    "export ORCA_SCRATCH=/pfs/data6/user_scratch/\n"
)


def _lustre(monkeypatch):
    monkeypatch.setattr(hpc_io_check.MountsSource, "read",
                        lambda self: LUSTRE_MOUNTS)


def test_returns_refused_and_warned_tuple(monkeypatch, tmp_path: Path):
    _lustre(monkeypatch)
    script = tmp_path / "run.sh"
    script.write_text(SCRATCH_SCRIPT, encoding="utf-8")
    refused, warned = hpc_io_check.check_submit(script, env={"HOME": "/home/x"})
    assert refused == []
    assert len(warned) == 1
    assert "scratch_on_lustre" in warned[0]


def test_refuse_becomes_refused_entry(monkeypatch, tmp_path: Path):
    _lustre(monkeypatch)
    script = tmp_path / "run.sh"
    script.write_text(WALK_SCRIPT, encoding="utf-8")
    refused, warned = hpc_io_check.check_submit(script, env={"HOME": "/home/x"})
    assert len(refused) == 1
    assert "tree_walk_on_lustre" in refused[0]
    assert warned == []


def test_accepts_text_directly(monkeypatch):
    _lustre(monkeypatch)
    refused, _ = hpc_io_check.check_submit(
        WALK_SCRIPT, env={"HOME": "/home/x"})
    assert len(refused) == 1


def test_no_lustre_mounts_means_clean(monkeypatch, tmp_path: Path):
    monkeypatch.setattr(hpc_io_check.MountsSource, "read",
                        lambda self: "tmpfs / tmpfs rw 0 0\n")
    script = tmp_path / "run.sh"
    script.write_text(WALK_SCRIPT, encoding="utf-8")
    refused, warned = hpc_io_check.check_submit(script, env={"HOME": "/home/x"})
    assert refused == []
    assert warned == []


def test_missing_script_file_is_an_error_not_a_pass(monkeypatch, tmp_path: Path):
    """A path that cannot be read must not silently pass the check."""
    _lustre(monkeypatch)
    with pytest.raises(OSError):
        hpc_io_check.check_submit(tmp_path / "does_not_exist.sh",
                                  env={"HOME": "/home/x"})
