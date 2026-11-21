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


# --- cd-relative sbatch (operator live run L4 class, 2026-11-21) ---

def _cd_env(monkeypatch):
    _lustre(monkeypatch)
    # HOME not on Lustre; the walk finding must come from the script
    # that the cd target resolves to, not from anything else.
    return {"HOME": "/home/x"}


def test_sbatch_gate_catches_walk_after_cd_and(monkeypatch, tmp_path: Path):
    """`cd <dir> && sbatch walk.sh`: the script must be resolved against
    the cd target, not the process cwd (unpatched: FileNotFoundError ->
    allow, the submit reaches sbatch)."""
    env = _cd_env(monkeypatch)
    sub = tmp_path / "sub"
    sub.mkdir()
    (sub / "walk.sh").write_text(WALK_SCRIPT, encoding="utf-8")
    verdict, text, may_run = hpc_io_check.sbatch_gate(
        f"cd {sub} && sbatch walk.sh", env=env)
    assert verdict == "refuse", (verdict, text)
    assert "tree_walk_on_lustre" in text
    assert may_run is False


def test_sbatch_gate_catches_walk_after_cd_semicolon(monkeypatch,
                                                     tmp_path: Path):
    """Same hole, `;` separator: `cd <dir>; sbatch walk.sh`."""
    env = _cd_env(monkeypatch)
    sub = tmp_path / "sub"
    sub.mkdir()
    (sub / "walk.sh").write_text(WALK_SCRIPT, encoding="utf-8")
    verdict, text, may_run = hpc_io_check.sbatch_gate(
        f"cd {sub}; sbatch walk.sh", env=env)
    assert verdict == "refuse", (verdict, text)
    assert may_run is False


def test_sbatch_gate_subshell_cd(monkeypatch, tmp_path: Path):
    """Subshell form: `(cd <dir> && sbatch walk.sh)` — the cd applies
    inside the subshell, the gate must still resolve the script."""
    env = _cd_env(monkeypatch)
    sub = tmp_path / "sub"
    sub.mkdir()
    (sub / "walk.sh").write_text(WALK_SCRIPT, encoding="utf-8")
    verdict, text, _ = hpc_io_check.sbatch_gate(
        f"(cd {sub} && sbatch walk.sh)", env=env)
    assert verdict == "refuse", (verdict, text)


def test_sbatch_gate_cd_target_missing_script_refused_narrowly(
        monkeypatch, tmp_path: Path):
    """A relative script that names nothing FROM THE CD TARGET is not
    silently passed: refuse narrowly (sbatch would fail, but the gate
    cannot tell what the script would have done — fail closed)."""
    env = _cd_env(monkeypatch)
    sub = tmp_path / "sub"
    sub.mkdir()
    verdict, text, may_run = hpc_io_check.sbatch_gate(
        f"cd {sub} && sbatch walk.sh", env=env)
    assert verdict == "refuse", (verdict, text)
    assert "not found from the cd target" in text
    assert may_run is False


def test_sbatch_gate_clean_script_after_cd_still_allowed(
        monkeypatch, tmp_path: Path):
    """No new refusal: a CLEAN relative script under the cd target
    still runs (ordinary behaviour must keep working)."""
    env = _cd_env(monkeypatch)
    sub = tmp_path / "sub"
    sub.mkdir()
    (sub / "clean.sh").write_text("#!/bin/bash\necho ok\n",
                                  encoding="utf-8")
    verdict, text, may_run = hpc_io_check.sbatch_gate(
        f"cd {sub} && sbatch clean.sh", env=env)
    assert verdict == "allow", (verdict, text)
    assert may_run is True


def test_sbatch_gate_cd_to_scratch_subdir_clean(monkeypatch, tmp_path: Path):
    """`cd <dir> && sbatch /abs/walk.sh`: the ABSOLUTE script is checked
    as before; a cd in the same line must not disable that check."""
    env = _cd_env(monkeypatch)
    walk = tmp_path / "walk.sh"
    walk.write_text(WALK_SCRIPT, encoding="utf-8")
    verdict, text, _ = hpc_io_check.sbatch_gate(
        f"cd {tmp_path / 'sub'} && sbatch {walk}", env=env)
    assert verdict == "refuse", (verdict, text)
