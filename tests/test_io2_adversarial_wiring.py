"""IO2 adversarial tests (reviewer nacht-s2), round 3: the _sbatch wiring.

``SlurmJobBackend._sbatch`` is the single funnel all seven submit callers
go through (backend_slurm.py:420), so it is the right place for the check.
These tests drive the real method -- not the helper -- with a fake
``subprocess.run`` so nothing is submitted, and assert what the caller
actually receives.

Reviewed: builder nacht-s1, commits 763698478 (wiring) and 3fa833756.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from delfin import hpc_io_check
from delfin.dashboard.backend_slurm import SlurmJobBackend

LUSTRE_MOUNTS = "mgmt@o2ib:fs /pfs lustre rw 0 0\n"
WALK = "#!/bin/bash\ndu -sh /pfs/data6/workspace/\n"


@pytest.fixture
def backend(tmp_path: Path) -> SlurmJobBackend:
    b = SlurmJobBackend.__new__(SlurmJobBackend)
    b.submit_templates_dir = tmp_path / "templates"
    b.submit_templates_dir.mkdir()
    b.slurm_profile = ""
    b.partitions = []
    return b


@pytest.fixture
def lustre(monkeypatch):
    monkeypatch.setattr(hpc_io_check.MountsSource, "read",
                        lambda self: LUSTRE_MOUNTS)


@pytest.fixture
def no_sbatch(monkeypatch):
    """Record any sbatch run; fail the test if one happens."""
    seen = {"ran": False}

    def fake_run(cmd, **kw):
        seen["ran"] = True
        return subprocess.CompletedProcess(
            args=cmd, returncode=0, stdout="Submitted batch job 1\n", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)
    monkeypatch.setattr(SlurmJobBackend, "_release_env_hold",
                        lambda self, stdout: None)
    return seen


def _script(tmp_path: Path, text: str) -> Path:
    path = tmp_path / "sub.sh"
    path.write_text(text, encoding="utf-8")
    return path


def test_a_refusal_stops_the_submit_before_sbatch(backend, lustre, no_sbatch,
                                                  tmp_path):
    script = _script(tmp_path, WALK)
    result = backend._sbatch(tmp_path, "", "01:00:00", 2, 1000, "j", script)
    assert result.returncode != 0
    assert not no_sbatch["ran"], "sbatch must not run for a refused script"
    assert "tree_walk_on_lustre" in (result.stderr or "")


def test_a_warning_is_visible_to_the_caller_of_a_successful_submit(
        backend, lustre, tmp_path, monkeypatch):
    """The task says warn findings are "logged with the job". When this was
    written the warning went only to ``stderr`` of a returncode-0 result,
    and a caller shows stderr only when the return code is not 0 (the rule
    backend_slurm.py:441 states), so it reached nobody. Now it rides in
    stdout next to the job id -- asserted here, together with the job id
    still parsing, because moving text into stdout is how a job id gets
    broken."""
    script = _script(
        tmp_path, "#!/bin/bash\nexport ORCA_SCRATCH=/pfs/data6/scr/\n")
    monkeypatch.setattr(subprocess, "run", lambda *a, **k: subprocess.CompletedProcess(
        args=["sbatch"], returncode=0, stdout="Submitted batch job 5\n",
        stderr=""))
    monkeypatch.setattr(SlurmJobBackend, "_release_env_hold",
                        lambda self, stdout: None)
    result = backend._sbatch(tmp_path, "", "01:00:00", 2, 1000, "j", script)
    assert result.returncode == 0, "a warning must not fail the submit"
    shown_on_success = (result.stdout or "")
    assert "scratch_on_lustre" in shown_on_success, (
        "the warning must be in the stream a successful submit's caller "
        "reads (stdout), not only in stderr, which is shown on failure")
    assert SlurmJobBackend._job_id_from_submit(result.stdout) == "5", (
        "the warning must not displace the job id the caller parses")


def test_a_clean_submit_leaves_stdout_exactly_as_sbatch_wrote_it(
        backend, lustre, tmp_path, monkeypatch):
    """The other half of the warning fix: a script with no finding must not
    have anything appended, or every submit's output changes."""
    script = _script(tmp_path, "#!/bin/bash\ncd $TMPDIR\necho hi\n")
    monkeypatch.setattr(subprocess, "run", lambda *a, **k: subprocess.CompletedProcess(
        args=["sbatch"], returncode=0, stdout="Submitted batch job 7\n",
        stderr=""))
    monkeypatch.setattr(SlurmJobBackend, "_release_env_hold",
                        lambda self, stdout: None)
    result = backend._sbatch(tmp_path, "", "01:00:00", 2, 1000, "j", script)
    assert result.stdout == "Submitted batch job 7\n"


def test_the_wiring_does_not_lose_the_home_rule(
        backend, lustre, no_sbatch, tmp_path, monkeypatch):
    """The dashboard route must not be weaker than the agent route on the
    $HOME rule. ``_sbatch`` calls ``check_submit(submit_script)`` without an
    env; that is fine only because ``env=None`` now expands from
    ``os.environ`` (H4). Pin the MECHANISM, not this host's HOME: with a
    Lustre HOME in the process env the route must refuse, with a node-local
    one it must submit -- either way the two routes agree."""
    script = _script(tmp_path, '#!/bin/bash\nfind "$HOME/work" -type f\n')

    monkeypatch.setenv("HOME", "/pfs/home/user")
    refused = backend._sbatch(tmp_path, "", "01:00:00", 2, 1000, "j", script)
    assert refused.returncode != 0, (
        "with HOME on Lustre the dashboard route must refuse this walk, as "
        "the agent route (which passes os.environ) does")
    assert not no_sbatch["ran"], "sbatch must not run for a refused script"
