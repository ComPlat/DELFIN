"""A script in a subdirectory imports the INSTALLED package, not the
worktree's — import_origin names that before it turns a green fix red.

Reproduces the false red of 2026-09-26 (.gate/probe_login_shell_full.py):
python puts the SCRIPT's directory on sys.path[0], never the worktree, so
`import delfin` resolved through the editable install into the main
checkout. Here the same trap is staged with a throwaway package: a
mini source package in the workspace and a second, "installed" one that
wins the import unless the command does something about it.
"""
import os
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from delfin.agent import import_origin  # noqa: E402


@pytest.fixture()
def stage(tmp_path):
    """A workspace with a source package and an 'installed' twin."""
    ws = tmp_path / "workspace"
    (ws / "minipkg").mkdir(parents=True)
    (ws / "minipkg" / "__init__.py").write_text("VALUE = 'worktree'\n")
    installed = tmp_path / "installed"
    (installed / "minipkg").mkdir(parents=True)
    (installed / "minipkg" / "__init__.py").write_text("VALUE = 'installed'\n")
    # A subdirectory a probe script would live in.
    (ws / ".gate").mkdir()
    (ws / ".gate" / "probe.py").write_text(
        "import minipkg, sys; print(minipkg.VALUE); print(minipkg.__file__)\n")
    (ws / "probe_root.py").write_text(
        "import minipkg; print(minipkg.VALUE)\n")
    return {"ws": ws, "installed": installed}


def _run(args, cwd, extra_env=None):
    env = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}
    if extra_env:
        env.update(extra_env)
    return subprocess.run(
        [sys.executable] + args, cwd=cwd, env=env,
        capture_output=True, text=True, timeout=60)


def test_the_trap_is_real_subdir_script_imports_the_installed_copy(stage):
    """Control: without any fix the subdirectory script measures the wrong
    code. This is the false red, staged — it must fail BEFORE any fix and
    is what the note is for."""
    out = _run([".gate/probe.py"], cwd=stage["ws"],
               extra_env={"PYTHONPATH": str(stage["installed"])})
    assert out.returncode == 0, out.stderr
    first = out.stdout.splitlines()[0]
    # The trap, affirmed: the subdirectory script measures the INSTALLED
    # code even though the worktree holds a package of the same name. If
    # this assertion ever fails, the staging no longer reproduces the
    # false red and every other test here is measuring nothing.
    assert first == "installed", (
        f"staging broken: probe got {out.stdout!r} -- the trap this "
        "module guards against did not occur")


def test_shadowed_packages_names_the_installed_location(stage):
    shadow = import_origin.shadowed_packages(
        stage["ws"], extra_paths=[str(stage["installed"])])
    assert set(shadow) == {"minipkg"}
    assert "installed" in shadow["minipkg"]


def test_shadowed_packages_empty_when_worktree_wins(stage):
    # No installed twin: the probe's path extension is the worktree itself.
    shadow = import_origin.shadowed_packages(
        stage["ws"], extra_paths=[str(stage["ws"])])
    assert shadow == {}


def test_no_shadow_no_note(stage):
    shadow = import_origin.shadowed_packages(
        stage["ws"], extra_paths=[str(stage["ws"])])
    assert shadow == {}
    note = import_origin.note_for_command(
        "python3 .gate/probe.py", cwd=str(stage["ws"]), workspace=str(stage["ws"]))
    assert note == ""


class TestNoteFires:
    def _note(self, stage, cmd):
        return import_origin.note_for_command(
            cmd, cwd=str(stage["ws"]), workspace=str(stage["ws"]),
            _shadow={"minipkg": str(stage["installed"] / "minipkg" / "__init__.py")})

    def test_subdir_script(self, stage):
        note = self._note(stage, "python3 .gate/probe.py")
        assert note and "minipkg" in note and "PYTHONPATH" in note

    def test_absolute_subdir_script(self, stage):
        script = str(stage["ws"] / ".gate" / "probe.py")
        note = self._note(stage, f"python3 {script}")
        assert note and "minipkg" in note

    def test_python_with_flag_then_script(self, stage):
        note = self._note(stage, "python3 -B .gate/probe.py")
        assert note and "minipkg" in note


class TestNoteStaysQuiet:
    def _note(self, stage, cmd):
        return import_origin.note_for_command(
            cmd, cwd=str(stage["ws"]), workspace=str(stage["ws"]),
            _shadow={"minipkg": "/somewhere/installed/minipkg/__init__.py"})

    def test_python_dash_c(self, stage):
        assert self._note(stage, 'python3 -c "import minipkg"') == ""

    def test_python_dash_m(self, stage):
        assert self._note(stage, "python3 -m minipkg.cli") == ""

    def test_pytest(self, stage):
        assert self._note(stage, "python3 -m pytest tests/ -q") == ""
        assert self._note(stage, "pytest -q") == ""

    def test_pythonpath_already_set(self, stage):
        assert self._note(stage, "PYTHONPATH=. python3 .gate/probe.py") == ""

    def test_script_in_worktree_root(self, stage):
        assert self._note(stage, "python3 probe_root.py") == ""

    def test_non_python_command(self, stage):
        assert self._note(stage, "ls -la .gate") == ""

    def test_no_python_at_all(self, stage):
        assert self._note(stage, "grep -rn minipkg .gate") == ""
