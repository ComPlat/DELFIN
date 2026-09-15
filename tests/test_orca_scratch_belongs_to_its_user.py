"""The default ORCA scratch directory is the user's own.

It was one directory for everybody under the system temp directory. On a
shared login node the first user to run created it 0755, and every other
user's ORCA run failed with "Permission denied" creating its run folder
inside -- observed on 2026-09-15 with /scratch/delfin_orca_scratch owned by
another account.
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
    monkeypatch.setattr(tempfile, "tempdir", str(tmp_path))
    monkeypatch.setattr(orca, "_RUN_SCRATCH_DIR", None)
    yield tmp_path
    orca._RUN_SCRATCH_DIR = None


@pytest.mark.skipif(hasattr(os, "geteuid") and os.geteuid() == 0,
                    reason="root writes into any directory")
def test_another_users_scratch_directory_does_not_block_a_run(fresh_scratch):
    theirs = fresh_scratch / "delfin_orca_scratch"
    theirs.mkdir()
    theirs.chmod(stat.S_IRUSR | stat.S_IXUSR | stat.S_IRGRP | stat.S_IXGRP
                 | stat.S_IROTH | stat.S_IXOTH)
    try:
        run_dir = orca._ensure_orca_scratch_dir()
    finally:
        theirs.chmod(0o755)
    assert run_dir.is_dir()
    assert run_dir.parent.name == orca._user_scratch_name()
    assert run_dir.parent.parent == fresh_scratch


def test_the_name_is_safe_for_any_user_name(monkeypatch):
    import getpass
    monkeypatch.setattr(getpass, "getuser", lambda: "ka ew/7404")
    assert orca._user_scratch_name() == "delfin_orca_scratch-ka_ew_7404"
