"""The launcher finds DELFIN behind a linked console script or in an
environment named delfin or delfin-voila, prefers one that can serve the
dashboard, and names an installation that is too old or lacks voila instead
of saying "Not found" or failing inside tmux.

A site-wide install is often a virtual environment under /opt whose console
scripts are linked into /usr/local/bin. `delfin-voila` is then on PATH, but
no python lies beside the link, and the server's plain python3 does not have
DELFIN: the discovery reported "Not found" although the person could run
delfin-voila by hand. The link is now followed and the script's #! line read.

When the only DELFIN found lacks the agent dashboard (an install from before
it), the message names that installation, so the person knows what to update
rather than hunting for something that "is not there".
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1] / "tools" / "windows-launcher"


def _site_install(tmp_path: Path, *, with_agent: bool) -> tuple[Path, Path]:
    """/opt/<venv> with DELFIN in its site-packages and delfin-voila linked
    into a bin directory of its own, like /usr/local/bin."""
    venv = tmp_path / "opt" / "delfin311"
    subprocess.run([sys.executable, "-m", "venv", "--without-pip", str(venv)], check=True)
    site = next((venv / "lib").glob("python3*/site-packages"))
    modules = ["delfin/__init__.py", "delfin/cli_voila.py"]
    if with_agent:
        modules += ["delfin/agent/__init__.py", "delfin/agent/where.py",
                    "delfin/dashboard/__init__.py", "delfin/dashboard/session.py", "voila/__init__.py"]
    for module in modules:
        path = site / module
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("")
    (site / "delfin/cli_voila.py").write_text("def main(argv=None):\n    print('--strict-port')\n")
    script = venv / "bin" / "delfin-voila"
    script.write_text(f"#!{venv / 'bin' / 'python3'}\nfrom delfin.cli_voila import main\nmain()\n")
    # The link is followed to the environment, whose python is taken first.
    python = venv / "bin" / "python"
    script.chmod(0o755)
    local_bin = tmp_path / "usr-local-bin"
    local_bin.mkdir()
    (local_bin / "delfin-voila").symlink_to(script)
    return local_bin, python


def _run(tmp_path: Path, local_bin: Path):
    home = tmp_path / "home"
    home.mkdir(exist_ok=True)
    payload = tmp_path / "copied" / "remote_launcher.py"
    payload.parent.mkdir()
    payload.write_text('import sys; print("CHOSEN=" + sys.executable)\n')
    env = {k: v for k, v in os.environ.items()
           if k not in {"VIRTUAL_ENV", "CONDA_PREFIX", "PYTHONPATH", "DELFIN_REPO", "MAMBA_ROOT_PREFIX"}}
    env.update(PATH=f"{local_bin}:/usr/bin:/bin", HOME=str(home))
    return subprocess.run(["bash", str(ROOT / "remote_bootstrap.sh"), "0", "", "", str(payload), "dashboard"],
                          env=env, cwd=str(home), capture_output=True, text=True, timeout=60)


def test_a_linked_console_script_leads_to_its_environment(tmp_path):
    local_bin, python = _site_install(tmp_path, with_agent=True)

    result = _run(tmp_path, local_bin)

    assert result.returncode == 0, result.stderr
    assert f"DELFIN: Environment found: {python}" in result.stdout
    assert "CHOSEN=" in result.stdout


def test_an_install_without_the_agent_is_named_as_too_old(tmp_path):
    local_bin, python = _site_install(tmp_path, with_agent=False)

    result = _run(tmp_path, local_bin)

    assert result.returncode == 1
    assert "too old for this launcher" in result.stderr
    assert str(python) in result.stderr
    assert "Not found" not in result.stderr


def _environment(prefix: Path, *, voila: bool = True) -> Path:
    """A Python environment at prefix with the DELFIN modules the launcher
    needs, and voila (which serves the dashboard) unless told otherwise."""
    subprocess.run([sys.executable, "-m", "venv", "--without-pip", str(prefix)], check=True)
    site = next((prefix / "lib").glob("python3*/site-packages"))
    modules = ["delfin/__init__.py", "delfin/cli_voila.py", "delfin/agent/__init__.py",
               "delfin/agent/where.py", "delfin/dashboard/__init__.py", "delfin/dashboard/session.py"]
    for module in modules + (["voila/__init__.py"] if voila else []):
        path = site / module
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("")
    (site / "delfin/cli_voila.py").write_text("def main(argv=None):\n    print('--strict-port')\n")
    return prefix / "bin" / "python"


def test_a_micromamba_environment_named_delfin_is_found_without_a_location(tmp_path):
    # A server whose DELFIN lives in ~/micromamba/envs/delfin, not activated
    # by the login shell: before, the location had to be typed in by hand.
    python = _environment(tmp_path / "home" / "micromamba" / "envs" / "delfin")
    empty_bin = tmp_path / "bin"
    empty_bin.mkdir()

    result = _run(tmp_path, empty_bin)

    assert result.returncode == 0, result.stderr
    assert f"DELFIN: Environment found: {python}" in result.stdout


def test_two_environments_named_delfin_are_listed_not_guessed(tmp_path):
    first = _environment(tmp_path / "home" / "micromamba" / "envs" / "delfin")
    second = _environment(tmp_path / "home" / "miniconda3" / "envs" / "delfin")
    empty_bin = tmp_path / "bin"
    empty_bin.mkdir()

    result = _run(tmp_path, empty_bin)

    assert result.returncode == 1
    assert "Several DELFIN installations" in result.stderr
    assert str(first) in result.stderr and str(second) in result.stderr


def test_the_environment_that_can_serve_the_dashboard_wins(tmp_path):
    # Measured on a server: envs/delfin had DELFIN but no voila, envs/delfin-voila
    # had both. The first was taken, and the dashboard died inside tmux with a
    # message about --strict-port.
    _environment(tmp_path / "home" / "micromamba" / "envs" / "delfin", voila=False)
    complete = _environment(tmp_path / "home" / "micromamba" / "envs" / "delfin-voila")
    empty_bin = tmp_path / "bin"
    empty_bin.mkdir()

    result = _run(tmp_path, empty_bin)

    assert result.returncode == 0, result.stderr
    assert f"DELFIN: Environment found: {complete}" in result.stdout


def test_an_environment_without_voila_is_named_before_anything_starts(tmp_path):
    partial = _environment(tmp_path / "home" / "micromamba" / "envs" / "delfin", voila=False)
    empty_bin = tmp_path / "bin"
    empty_bin.mkdir()

    result = _run(tmp_path, empty_bin)

    assert result.returncode == 1
    assert "cannot serve the dashboard (voila is not installed)" in result.stderr
    assert str(partial) in result.stderr
    assert "CHOSEN=" not in result.stdout
    assert "--strict-port" not in result.stderr
