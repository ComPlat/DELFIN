"""An ML-potential job computes something, on the GPU it was given.

submit_mlp sent jobs with DELFIN_MODE=mlp, and the job runner had no such
mode: every one ended in "Unknown mode: mlp". The MLP installer checked
PyTorch by finding a torch directory, so a remnant without its libraries
counted as installed, and it installed the CPU-only build, so a job given a
GPU computed on the CPU anyway.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from delfin.dashboard import local_runner
from delfin.tools._types import StepResult, StepStatus

REPO = Path(__file__).resolve().parents[1]
WATER = "3\nwater\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH -0.24 0.93 0.0\n"


@pytest.fixture
def mlp_job(tmp_path, monkeypatch):
    xyz = tmp_path / "water.xyz"
    xyz.write_text(WATER)
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    monkeypatch.chdir(run_dir)
    for key in ("DELFIN_MLP_TASK", "DELFIN_MLP_DEVICE", "DELFIN_MLP_FMAX", "DELFIN_MLP_STEPS"):
        monkeypatch.delenv(key, raising=False)
    monkeypatch.setenv("DELFIN_XYZ_FILE", str(xyz))
    monkeypatch.setenv("DELFIN_MLP_BACKEND", "ani2x")
    monkeypatch.setenv("DELFIN_CHARGE", "0")
    monkeypatch.setenv("DELFIN_MULT", "1")
    monkeypatch.setenv("DELFIN_PAL", "4")
    monkeypatch.setattr(local_runner, "_mlp_device", lambda: ("cuda", "NVIDIA H100"))
    calls = []

    def fake_run_step(step_name, *, geometry=None, cores=1, work_dir=None, **kwargs):
        calls.append({"step": step_name, "cores": cores, **kwargs})
        work_dir.mkdir(parents=True, exist_ok=True)
        out = work_dir / "optimized.xyz"
        out.write_text(WATER)
        return StepResult(step_name, StepStatus.SUCCESS, geometry=out if step_name == "mlp_optimize" else None,
                          work_dir=work_dir, data={"energy_eV": -2079.1234, "converged": True, "n_steps": 12})

    monkeypatch.setattr("delfin.tools.run_step", fake_run_step)
    return run_dir, calls


def test_an_mlp_job_runs_the_optimisation_on_its_gpu(mlp_job):
    run_dir, calls = mlp_job

    assert local_runner._run_mode("mlp") == 0

    assert calls == [{"step": "mlp_optimize", "cores": 4, "backend": "ani2x", "charge": 0,
                      "mult": 1, "device": "cuda", "fmax": 0.05, "steps": 200}]
    result = json.loads((run_dir / "mlp_result.json").read_text())
    assert result["status"] == "success" and result["device"] == "cuda" and result["gpu"] == "NVIDIA H100"
    assert result["energy_eV"] == -2079.1234
    assert (run_dir / "water_mlp_opt.xyz").is_file()


def test_a_single_point_is_a_single_point(mlp_job, monkeypatch):
    run_dir, calls = mlp_job
    monkeypatch.setenv("DELFIN_MLP_TASK", "single-point")

    assert local_runner._run_mode("mlp") == 0
    assert calls[0]["step"] == "mlp_single_point" and "fmax" not in calls[0]


def test_a_missing_backend_says_how_to_install_it(mlp_job, monkeypatch, capsys):
    run_dir, _ = mlp_job

    def missing(*args, **kwargs):
        raise ModuleNotFoundError("No module named 'torchani'")

    monkeypatch.setattr("delfin.tools.run_step", missing)

    assert local_runner._run_mode("mlp") == 1
    said = capsys.readouterr().out
    assert "python -m delfin.installer --install ani2x" in said
    assert json.loads((run_dir / "mlp_result.json").read_text())["status"] == "failed"


def test_bad_input_is_refused_before_anything_runs(mlp_job, monkeypatch):
    _, calls = mlp_job
    monkeypatch.setenv("DELFIN_MLP_TASK", "dance")
    assert local_runner._run_mode("mlp") == 1
    monkeypatch.setenv("DELFIN_MLP_TASK", "optimize")
    monkeypatch.setenv("DELFIN_XYZ_FILE", "/nonexistent/mol.xyz")
    assert local_runner._run_mode("mlp") == 1
    assert calls == []


def test_the_device_is_the_cpu_when_asked_or_when_there_is_no_gpu(monkeypatch):
    monkeypatch.setenv("DELFIN_MLP_DEVICE", "cpu")
    assert local_runner._mlp_device() == ("cpu", "")


def _shell_function(script: Path, name: str) -> str:
    text = script.read_text(encoding="utf-8")
    return name + "() {" + text.split(name + "() {", 1)[1].split("\n}\n", 1)[0] + "\n}\n"


def _fake_python(tmp_path: Path, imports: bool, found: bool) -> Path:
    fake = tmp_path / "python"
    fake.write_text(
        "#!/bin/sh\n"
        'case "$*" in\n'
        f'  *find_spec*) exit {0 if found else 1} ;;\n'
        f'  "-c import torch") exit {0 if imports else 1} ;;\n'
        '  *"-m pip"*) echo "$*" >> "$(dirname "$0")/pip.log"; exit 0 ;;\n'
        "esac\n"
        "exit 1\n")
    fake.chmod(0o755)
    return fake


@pytest.mark.parametrize("variant,expect_index", [("", False), ("cpu", True)])
def test_a_torch_that_does_not_import_is_installed_again(tmp_path, variant, expect_index):
    """A torch directory without its libraries was taken as installed."""
    script = REPO / "delfin" / "mlp_tools" / "install_mlp_tools.sh"
    fake = _fake_python(tmp_path, imports=False, found=True)
    body = ('log() { :; }\nwarn() { echo "WARN $*" >&2; }\n'
            + _shell_function(script, "python_has_module")
            + _shell_function(script, "check_pytorch")
            + f'check_pytorch "{fake}"\n')
    done = subprocess.run(["bash", "-c", body], capture_output=True, text=True, timeout=30,
                          env={"PATH": os.environ["PATH"], "LOG_DIR": str(tmp_path),
                               "DELFIN_TORCH_VARIANT": variant})

    pip = (tmp_path / "pip.log").read_text().splitlines()
    assert "does not import" in done.stderr
    assert "--force-reinstall" in pip[0], "the partial install is replaced, not trusted"
    assert all(("download.pytorch.org/whl/cpu" in line) == expect_index for line in pip), pip


def test_a_working_torch_is_left_alone(tmp_path):
    script = REPO / "delfin" / "mlp_tools" / "install_mlp_tools.sh"
    fake = _fake_python(tmp_path, imports=True, found=True)
    body = ('log() { :; }\nwarn() { :; }\n' + _shell_function(script, "python_has_module")
            + _shell_function(script, "check_pytorch") + f'check_pytorch "{fake}"\n')
    subprocess.run(["bash", "-c", body], capture_output=True, text=True, timeout=30,
                   env={"PATH": os.environ["PATH"], "LOG_DIR": str(tmp_path)})
    assert not (tmp_path / "pip.log").exists()
