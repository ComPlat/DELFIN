"""IO2 adversarial tests (reviewer nacht-s2), round 5: the REAL consumers of
the submit result.

D1 (19042b858) moved a warn finding into the stdout of a successful
``_sbatch``. That is the right stream -- stderr is only shown on failure --
but stdout is not free: 15 dashboard call sites parse the job id out of it
with ``stdout.strip().split()[-1]``, i.e. the LAST whitespace token.

    tab_submit.py            714, 873, 950, 1566, 1751, 1889, 2009, 2086
    tab_calculations_browser.py 5407, 5473, 10822, 10997
    tab_orca_builder.py      2543
    tab_recalc.py            156
    ts_goat_panel.py         1066

With the warning appended, the last token is a word from the warning text,
so the dashboard hands a bogus "job id" to auto-watch and to squeue. The
backend's own ``_job_id_from_submit`` (a regex) is unaffected -- which is
exactly why a test that only asks the backend would not catch this.

Measured on builder nacht-s1's 19042b858.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from delfin import hpc_io_check
from delfin.dashboard.backend_slurm import SlurmJobBackend

LUSTRE_MOUNTS = "mgmt@o2ib:fs /pfs lustre rw 0 0\n"
WARN_SCRIPT = "#!/bin/bash\nexport ORCA_SCRATCH=/pfs/data6/scr/\n"
CLEAN_SCRIPT = "#!/bin/bash\ncd $TMPDIR\necho hi\n"


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


def _submit(backend, tmp_path, monkeypatch, text):
    script = tmp_path / "s.sh"
    script.write_text(text, encoding="utf-8")
    monkeypatch.setattr(subprocess, "run", lambda *a, **k: subprocess.CompletedProcess(
        args=["sbatch"], returncode=0,
        stdout="Submitted batch job 4242\n", stderr=""))
    monkeypatch.setattr(SlurmJobBackend, "_release_env_hold",
                        lambda self, stdout: None)
    return backend._sbatch(tmp_path, "", "01:00:00", 2, 1000, "j", script)


def test_a_successful_submit_still_yields_its_job_id_to_the_dashboard(
        backend, lustre, tmp_path, monkeypatch):
    """The exact expression 15 call sites use, on the real returned stdout.
    On 19042b858 a warned submit yields "$TMPDIR)" instead of "4242": the
    job would be watched and queried under a job id that does not exist."""
    result = _submit(backend, tmp_path, monkeypatch, WARN_SCRIPT)
    assert result.returncode == 0
    stdout = result.stdout or ""
    dashboard_id = stdout.strip().split()[-1] if stdout.strip() else "(unknown)"
    assert dashboard_id == "4242", (
        "the dashboard parses the job id as the last token of stdout; a "
        "warning appended after it displaces the id")


def test_the_warning_is_still_reported(backend, lustre, tmp_path, monkeypatch):
    """The reason D1 exists: the warning must be somewhere the caller of a
    successful submit actually reads. Both may be true at once -- the id is
    last, or the warning is carried in its own field."""
    result = _submit(backend, tmp_path, monkeypatch, WARN_SCRIPT)
    assert "scratch_on_lustre" in (result.stdout or "") + (result.stderr or "")


def test_a_clean_submit_yields_its_job_id(backend, lustre, tmp_path, monkeypatch):
    """Control: without a finding the expression works, so a failure above
    is caused by the warning and nothing else."""
    result = _submit(backend, tmp_path, monkeypatch, CLEAN_SCRIPT)
    assert result.stdout.strip().split()[-1] == "4242"


def test_the_backends_own_parser_agrees(backend, lustre, tmp_path, monkeypatch):
    """The two parsers must not disagree: the regex finds 4242 either way,
    which is how this defect stayed invisible to a backend-only test."""
    result = _submit(backend, tmp_path, monkeypatch, WARN_SCRIPT)
    assert SlurmJobBackend._job_id_from_submit(result.stdout) == "4242"
