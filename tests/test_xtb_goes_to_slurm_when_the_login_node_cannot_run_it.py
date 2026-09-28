"""xtb runs through SLURM when the login node cannot run it.

The model endpoint is only reachable from the login node, so the bench
runner stays there - but the cluster rule is that compute happens on
compute nodes. With DELFIN_CHEM_SLURM=1 the chemistry acceptance
scripts submit their xtb call via `sbatch --wait` instead of running it
in-process, and read the log back. The local path stays the default
so nothing changes where SLURM is not wanted.
"""

import importlib.util
import os
import sys
from pathlib import Path

import pytest

_SCRIPT = (Path(__file__).resolve().parents[1] /
           "delfin" / "agent" / "pack" / "benchmark" / "accept" /
           "chem_opt_is_a_minimum.py")


def _load():
    spec = importlib.util.spec_from_file_location("chem_opt", _SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["chem_opt"] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture()
def fake_sbatch(tmp_path, monkeypatch):
    """Record the batch script and pretend sbatch --wait ran the xtb
    command, producing a log with a marker line."""
    calls = []

    def fake_run(cmd, cwd=None, capture_output=True, text=True, timeout=None):
        calls.append((cmd, cwd))
        if cmd and str(cmd[0]).endswith("sbatch"):
            base = Path(str(cwd)) if cwd else Path(".")
            script = base / str(cmd[-1])
            body = script.read_text()
            log = script.parent / "m2-chem-accept-123.log"
            log.write_text("FAKE-SBATCH-RAN\n" + body)
            return type("P", (), {"returncode": 0, "stdout": str(log),
                                  "stderr": ""})()
        return type("P", (), {"returncode": 0, "stdout": "", "stderr": ""})()

    import subprocess
    monkeypatch.setattr(subprocess, "run", fake_run)
    return calls


def test_env_var_routes_xtb_through_sbatch(tmp_path, monkeypatch,
                                           fake_sbatch):
    monkeypatch.setenv("DELFIN_CHEM_SLURM", "1")
    mod = _load()
    out = mod._run(["/fake/xtb", "check.xyz", "--opt"], tmp_path)
    assert "FAKE-SBATCH-RAN" in out, (
        "with DELFIN_CHEM_SLURM=1 the xtb call must go through sbatch "
        f"--wait; got: {out[:200]}")


def test_without_env_var_xtb_runs_locally(tmp_path, monkeypatch,
                                          fake_sbatch):
    monkeypatch.delenv("DELFIN_CHEM_SLURM", raising=False)
    mod = _load()
    mod._run(["/fake/xtb", "check.xyz", "--opt"], tmp_path)
    submitted = [c for c in fake_sbatch
                 if "sbatch" in " ".join(map(str, c[0]))]
    assert not submitted, (
        "without the env var nothing may be submitted")
