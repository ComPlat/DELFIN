"""Venv staging is the default, and its fallbacks are what the job log says.

The template's Section 7 decided on DELFIN_STAGE_VENV with ":-0": a job sent
with plain sbatch (no dashboard, no profile env) ran its Python from the
network HOME, the very I/O the staging exists to avoid. The dashboard keeps
injecting DELFIN_STAGE_VENV=1 for every site; the template itself now agrees.

The fallbacks are not optional extras. Not enough room on the local disk,
or no tar that can be packed, must say so in the job log and run the venv
from where it is -- a calculation is never lost over where Python lives.
"""

import os
import pathlib
import subprocess

REPO = pathlib.Path(__file__).resolve().parents[1]
TEMPLATE = REPO / "delfin" / "submit_templates" / "submit_delfin.sh"


def _section7() -> str:
    text = TEMPLATE.read_text(encoding="utf-8")
    start = text.index("venv_cache_key() {")
    end = text.index("# ======================================================================\n# Section 8")
    return text[start:end]


def _bash(tmp_path, body, extra_env=None):
    script = tmp_path / "run.sh"
    script.write_text("set -euo pipefail\n" + _section7() + "\n" + body + "\n")
    env = {
        "HOME": str(tmp_path / "home"),
        "PATH": os.environ["PATH"],
        "DELFIN_VENV_CACHE_DIR": str(tmp_path / "cache"),
        "SLURM_JOB_ID": "77",
        "STAGE_BASE": str(tmp_path / "stage"),
    }
    env.update(extra_env or {})
    return subprocess.run(
        ["bash", str(script)], capture_output=True, text=True, timeout=120, env=env
    )


def _venv(tmp_path, packages=("numpy-1.26.4.dist-info",)):
    venv = tmp_path / "home" / "venv"
    site = venv / "lib" / "python3.11" / "site-packages"
    for package in packages:
        (site / package).mkdir(parents=True)
    (venv / "bin").mkdir()
    (venv / "bin" / "python").write_text("#!/bin/sh\n")
    (venv / "bin" / "python").chmod(0o755)
    (venv / "pyvenv.cfg").write_text("home = /usr/bin\n")
    return venv, site


def _stage(tmp_path, free_kb=None):
    stage = tmp_path / "stage"
    stage.mkdir()
    if free_kb is not None:
        shim = tmp_path / "shim"
        shim.mkdir()
        real_df = subprocess.run(["bash", "-c", "command -v df"], capture_output=True, text=True).stdout.strip()
        (shim / "df").write_text(
            f'#!/bin/sh\nif [ "$1" = "-Pk" ]; then echo "Filesystem 1024-blocks Used Available Capacity Mounted"\n'
            f'  echo "dev 1000000 100000 {free_kb} 1% {stage}"\n  exit 0\nfi\nexec "{real_df}" "$@"\n'
        )
        (shim / "df").chmod(0o755)
        return stage, shim
    return stage, None


def test_staging_is_the_default_without_any_env_set(tmp_path):
    """Plain sbatch, nothing exported: the venv still comes off the local disk."""
    venv, _ = _venv(tmp_path)
    stage, _ = _stage(tmp_path)
    done = _bash(tmp_path, f'true', extra_env={"DELFIN_VENV": str(venv)})

    assert done.returncode == 0, done.stderr
    assert "venv loaded from local SSD." in done.stdout, done.stdout


def test_too_little_room_says_so_and_runs_the_venv_from_home(tmp_path):
    venv, _ = _venv(tmp_path)
    stage, shim = _stage(tmp_path, free_kb=1)
    env_extra = {"DELFIN_VENV": str(venv), "PATH": f"{shim}:{os.environ['PATH']}"}
    done = _bash(tmp_path, "true", extra_env=env_extra)

    assert done.returncode == 0, done.stderr
    assert "running it from" in done.stdout or "running the venv from" in done.stdout, done.stdout
    assert "venv loaded from local SSD." not in done.stdout
    assert "Using venv directly" in done.stdout


def test_a_tar_that_cannot_be_packed_says_so_and_keeps_running(tmp_path):
    venv, _ = _venv(tmp_path, packages=("a-1.dist-info",))
    stage, _ = _stage(tmp_path)
    # DELFIN_VENV_TAR points at a nonexistent tar: packing is skipped, the
    # decision block finds no tar file and must fall back with a WARNING.
    done = _bash(tmp_path, "true", extra_env={"DELFIN_VENV": str(venv), "DELFIN_VENV_TAR": str(tmp_path / "cache" / "nope.tar")})

    assert done.returncode == 0, done.stderr
    assert "WARNING" in done.stdout
    assert "venv loaded from local SSD." not in done.stdout
    assert "Using venv directly" in done.stdout


def test_delfin_stage_venv_0_turns_staging_off(tmp_path):
    venv, _ = _venv(tmp_path)
    stage, _ = _stage(tmp_path)
    done = _bash(tmp_path, "true", extra_env={"DELFIN_VENV": str(venv), "DELFIN_STAGE_VENV": "0"})

    assert done.returncode == 0, done.stderr
    assert "venv loaded from local SSD." not in done.stdout
    assert "Using venv directly" in done.stdout
