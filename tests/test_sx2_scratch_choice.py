"""The ORCA scratch choice inside an agent session's command cage.

Inside a cage the session's scratch directory is bound read-write into
the cage and exported as TMPDIR and DELFIN_SESSION_SCRATCH
(api_client._session_scratch_dir), while the SLURM job's node-local tmp
directory (SLURM_TMPDIR, /scratch/slurm_tmpdir/job_<id>) is mounted
READ-ONLY -- a base the candidate list prefers after DELFIN_SCRATCH. The
cage marker is DELFIN_AGENT_SESSION (approval_guard.ENV_MARKER), set on
every agent child. So when a calculation is started from a command of an
agent session, the session's scratch must be preferred over the job's
read-only tmp dir, and a DELFIN_SCRATCH on Lustre the cage cannot write
(the workaround SCC asked us to stop) must not win either. Outside an
agent session nothing changes: DELFIN_SCRATCH stays first, SLURM_TMPDIR
keeps its place.
"""

from __future__ import annotations

import os
import stat
import tempfile

import pytest

from delfin import orca


CAGE_MARKER = "DELFIN_AGENT_SESSION"


@pytest.fixture
def fresh_scratch(tmp_path, monkeypatch):
    for var in ("ORCA_SCRDIR", "ORCA_TMPDIR", "DELFIN_SCRATCH", "SLURM_TMPDIR",
                "SLURM_JOB_ID", "DELFIN_RUN_TOKEN", CAGE_MARKER,
                "DELFIN_SESSION_SCRATCH"):
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
def test_in_a_cage_the_session_scratch_beats_the_read_only_job_tmp(
        fresh_scratch, monkeypatch):
    """The session scratch is writable inside the cage; the job's node-local
    tmp is not. The session scratch must win even though the job's tmp
    comes earlier in the classic candidate order."""
    session_scratch = fresh_scratch / "session_scratch"
    session_scratch.mkdir()
    monkeypatch.setenv(CAGE_MARKER, str(os.getpid()))
    monkeypatch.setenv("DELFIN_SESSION_SCRATCH", str(session_scratch))
    monkeypatch.setenv("SLURM_TMPDIR", str(_read_only(fresh_scratch / "job_tmp")))
    monkeypatch.setenv("SLURM_JOB_ID", "4711")
    run_dir = orca._ensure_orca_scratch_dir()
    assert run_dir.parent == session_scratch


@needs_non_root
def test_in_a_cage_the_session_scratch_beats_an_unwritable_delfin_scratch(
        fresh_scratch, monkeypatch):
    """A DELFIN_SCRATCH on Lustre that the cage cannot write (the workaround
    the operator had to use) must not be chosen either; the session scratch
    still wins."""
    session_scratch = fresh_scratch / "session_scratch"
    session_scratch.mkdir()
    monkeypatch.setenv(CAGE_MARKER, str(os.getpid()))
    monkeypatch.setenv("DELFIN_SESSION_SCRATCH", str(session_scratch))
    monkeypatch.setenv("DELFIN_SCRATCH", str(_read_only(fresh_scratch / "lustre")))
    monkeypatch.setenv("SLURM_TMPDIR", str(_read_only(fresh_scratch / "job_tmp")))
    monkeypatch.setenv("SLURM_JOB_ID", "4711")
    run_dir = orca._ensure_orca_scratch_dir()
    assert run_dir.parent == session_scratch


def test_without_the_cage_marker_the_session_variable_changes_nothing(
        fresh_scratch, monkeypatch):
    """DELFIN_SESSION_SCRATCH set WITHOUT the cage marker is ignored: the
    variable only exists inside the cage, and a stray value outside must
    not reorder the classic candidates."""
    delfin_scratch = fresh_scratch / "delfin_scratch"
    delfin_scratch.mkdir()
    session_scratch = fresh_scratch / "session_scratch"
    session_scratch.mkdir()
    monkeypatch.setenv("DELFIN_SESSION_SCRATCH", str(session_scratch))
    monkeypatch.setenv("DELFIN_SCRATCH", str(delfin_scratch))
    run_dir = orca._ensure_orca_scratch_dir()
    assert run_dir.parent == delfin_scratch


def test_outside_a_cage_a_writable_delfin_scratch_still_wins(fresh_scratch,
                                                             monkeypatch):
    """No cage marker, no session variable: DELFIN_SCRATCH keeps first
    place."""
    delfin_scratch = fresh_scratch / "delfin_scratch"
    delfin_scratch.mkdir()
    monkeypatch.setenv("DELFIN_SCRATCH", str(delfin_scratch))
    run_dir = orca._ensure_orca_scratch_dir()
    assert run_dir.parent == delfin_scratch


@needs_non_root
def test_outside_a_cage_the_job_tmpdir_keeps_its_place_over_the_session_var(
        fresh_scratch, monkeypatch):
    """Without the cage marker a writable job tmp dir still ranks above a
    DELFIN_SESSION_SCRATCH value -- the classic SLURM behaviour is
    preserved."""
    job_tmp = fresh_scratch / "job_tmp"
    job_tmp.mkdir()
    session_scratch = fresh_scratch / "session_scratch"
    session_scratch.mkdir()
    monkeypatch.setenv("DELFIN_SESSION_SCRATCH", str(session_scratch))
    monkeypatch.setenv("SLURM_TMPDIR", str(job_tmp))
    run_dir = orca._ensure_orca_scratch_dir()
    assert run_dir.parent == job_tmp
