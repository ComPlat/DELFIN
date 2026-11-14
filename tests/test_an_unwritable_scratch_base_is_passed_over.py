"""An ORCA scratch base that cannot be written is passed over, not fatal.

Inside a SLURM job the scratch base is the job's node-local tmp directory
(SLURM_TMPDIR, or /scratch/slurm_tmpdir/job_<id> when it exists). An agent
session's command sandbox mounts that directory read-only, so every ORCA
step started from an agent session inside a job failed before ORCA ran:
"[Errno 30] Read-only file system" creating the run folder -- seen in a live
three-step chain (xtb_opt -> orca_opt -> orca_sp) on 2026-11-14.
Now each base is tried in its order and one that cannot hold the run folder
is skipped with a warning; the next one (finally the user's own directory
under the system temp directory) is used. A writable base keeps its place.
"""

from __future__ import annotations

import os
import stat
import tempfile

import pytest

from delfin import orca


@pytest.fixture
def fresh_scratch(tmp_path, monkeypatch):
    for var in ("ORCA_SCRDIR", "ORCA_TMPDIR", "DELFIN_SCRATCH", "SLURM_TMPDIR",
                "SLURM_JOB_ID", "DELFIN_RUN_TOKEN"):
        monkeypatch.delenv(var, raising=False)
    system_tmp = tmp_path / "system_tmp"
    system_tmp.mkdir()
    monkeypatch.setattr(tempfile, "tempdir", str(system_tmp))
    monkeypatch.setattr(orca, "_RUN_SCRATCH_DIR", None)
    yield tmp_path
    orca._RUN_SCRATCH_DIR = None


def _read_only(path):
    path.mkdir()
    path.chmod(stat.S_IRUSR | stat.S_IXUSR)
    return path


needs_non_root = pytest.mark.skipif(
    hasattr(os, "geteuid") and os.geteuid() == 0,
    reason="root writes into any directory")


@needs_non_root
def test_a_read_only_job_tmpdir_falls_back_to_the_users_own_dir(fresh_scratch, monkeypatch):
    job_tmp = _read_only(fresh_scratch / "job_tmp")
    monkeypatch.setenv("SLURM_TMPDIR", str(job_tmp))
    try:
        run_dir = orca._ensure_orca_scratch_dir()
    finally:
        job_tmp.chmod(0o755)
    assert run_dir.is_dir()
    assert run_dir.parent.name == orca._user_scratch_name()
    assert run_dir.parent.parent == fresh_scratch / "system_tmp"
    assert not any(job_tmp.iterdir())


@needs_non_root
def test_the_next_writable_base_in_order_wins(fresh_scratch, monkeypatch):
    blocked = _read_only(fresh_scratch / "blocked")
    usable = fresh_scratch / "usable"
    monkeypatch.setenv("ORCA_SCRDIR", str(blocked))
    monkeypatch.setenv("DELFIN_SCRATCH", str(usable))
    monkeypatch.setenv("SLURM_TMPDIR", str(fresh_scratch / "job_tmp"))
    try:
        run_dir = orca._ensure_orca_scratch_dir()
    finally:
        blocked.chmod(0o755)
    assert run_dir.parent == usable


def test_a_writable_job_tmpdir_is_still_used(fresh_scratch, monkeypatch):
    job_tmp = fresh_scratch / "job_tmp"
    job_tmp.mkdir()
    monkeypatch.setenv("SLURM_TMPDIR", str(job_tmp))
    run_dir = orca._ensure_orca_scratch_dir()
    assert run_dir.parent == job_tmp
    assert run_dir.is_dir()
