"""IO2 phase 2: the dashboard's SLURM backend checks the submit script
BEFORE sbatch runs it.

Real call site: SlurmJobBackend.submit_delfin / submit_orca / ... all
funnel into _sbatch, which appends the script to the sbatch argv and
runs it. The check must run inside _sbatch (one place, all seven
callers covered): a "refuse" finding returns a CompletedProcess with
returncode 1 and the finding text in stderr (what the dashboard shows
where a failed submit shows); "warn" findings ride in the result
metadata without stopping the submit. DELFIN's own generated scripts
must pass clean — a refusal there is a template bug, not a reason to
weaken the rule.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from delfin import hpc_io_check
from delfin.dashboard.backend_slurm import SlurmJobBackend


LUSTRE_MOUNTS = "mgmt@o2ib:fs /pfs lustre rw 0 0\n"


@pytest.fixture()
def backend(tmp_path: Path) -> SlurmJobBackend:
    b = SlurmJobBackend.__new__(SlurmJobBackend)
    b.submit_templates_dir = tmp_path / "templates"
    b.submit_templates_dir.mkdir()
    b.slurm_profile = ""
    b.partitions = []
    return b


def _bad_script(tmp_path: Path, name: str = "bad.sh") -> Path:
    script = tmp_path / name
    script.write_text("#!/bin/bash\ndu -sh /pfs/data6/workspace/\n",
                      encoding="utf-8")
    return script


def test_sbatch_refuses_a_lustre_walking_script(backend, monkeypatch,
                                                tmp_path):
    monkeypatch.setattr(hpc_io_check.MountsSource, "read",
                        lambda self: LUSTRE_MOUNTS)
    script = _bad_script(tmp_path)
    called = {}
    monkeypatch.setattr(subprocess, "run",
                        lambda *a, **k: called.setdefault("run", True))
    result = backend._sbatch(tmp_path, "", "01:00:00", 2, 1000,
                             "j", script)
    assert result.returncode == 1
    assert "tree_walk_on_lustre" in (result.stderr or "")
    assert "run" not in called  # sbatch never executed


def test_sbatch_warn_rides_and_submit_proceeds(backend, monkeypatch,
                                               tmp_path):
    monkeypatch.setattr(hpc_io_check.MountsSource, "read",
                        lambda self: LUSTRE_MOUNTS)
    script = tmp_path / "warn.sh"
    script.write_text(
        "#!/bin/bash\nexport ORCA_SCRATCH=/pfs/data6/scr/\n",
        encoding="utf-8")
    import subprocess as sp
    monkeypatch.setattr(
        sp, "run",
        lambda *a, **k: sp.CompletedProcess(args=a[0], returncode=0,
                                            stdout="Submitted batch job 7\n",
                                            stderr=""))
    result = backend._sbatch(tmp_path, "", "01:00:00", 2, 1000,
                             "j", script)
    assert result.returncode == 0
    assert "scratch_on_lustre" in (result.stderr or "") + (result.stdout or "")


def test_sbatch_clean_script_submits_normally(backend, monkeypatch,
                                              tmp_path):
    monkeypatch.setattr(hpc_io_check.MountsSource, "read",
                        lambda self: LUSTRE_MOUNTS)
    script = tmp_path / "clean.sh"
    script.write_text(
        "#!/bin/bash\ncp -a inputs \"$TMPDIR\"/ && \"$TMPDIR\"/run\n",
        encoding="utf-8")
    import subprocess as sp
    seen = {}

    def fake_run(cmd, **kw):
        seen["cmd"] = cmd
        return sp.CompletedProcess(args=cmd, returncode=0,
                                   stdout="Submitted batch job 9\n",
                                   stderr="")

    monkeypatch.setattr(sp, "run", fake_run)
    monkeypatch.setattr(SlurmJobBackend, "_release_env_hold",
                        lambda self, stdout: None)
    result = backend._sbatch(tmp_path, "", "01:00:00", 2, 1000,
                             "j", script)
    assert result.returncode == 0
    assert str(script) in seen["cmd"]


def test_delphins_own_templates_pass_clean(backend, monkeypatch):
    """Every script the backend can generate must pass the check with no
    refusal. A refusal here is a DELFIN template bug to fix, not a reason
    to weaken the rule. The HOME is on Lustre so the home rule fires —
    the template must not rely on a script reading into it."""
    monkeypatch.delenv("DELFIN_SQUEUE_MIN_INTERVAL", raising=False)
    real_dir = Path(hpc_io_check.DEFAULT_MOUNTS)
    if not real_dir.is_file():
        pytest.skip("no /proc/mounts on this host")
    templates = Path(__file__).resolve().parents[1] / "delfin" / "submit_templates"
    for template in sorted(templates.glob("*.sh")):
        refused, warned = hpc_io_check.check_submit(
            template, env={"HOME": "/pfs/home/nobody"})
        assert refused == [], \
            f"DELFIN's own template {template.name} is refused: {refused}"
