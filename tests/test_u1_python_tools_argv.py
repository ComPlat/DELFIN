"""The py-group install runs the catalog command as argv, intact.

Input: a session interpreter whose path contains a space (a venv under
``~/Library/Application Support`` on macOS). Output: the argv handed to
subprocess.run. Splitting the display string on whitespace turned that
path into two words, so pip ran under a non-existent program.
"""

from __future__ import annotations

from delfin import installer


def test_a_session_python_with_a_space_stays_one_argument(monkeypatch):
    py = "/Users/me/Library/Application Support/delfin/venv/bin/python"
    monkeypatch.setattr(installer, "session_python", lambda: py)
    seen = []

    class _Done:
        returncode = 0
        stdout = stderr = ""

    def _run(argv, **kw):
        seen.append(argv)
        return _Done()

    monkeypatch.setattr(installer.subprocess, "run", _run)
    ok, _ = installer._install_python_tools([installer.find("pytest")])
    assert ok
    assert seen == [[py, "-m", "pip", "install", "pytest", "pytest-timeout"]]


def test_the_quoted_command_passes_the_gate_for_that_interpreter(monkeypatch):
    from delfin.agent.api_client import _refuse_unsafe_install
    py = "/Users/me/Library/Application Support/delfin/venv/bin/python"
    monkeypatch.setattr(installer, "session_python", lambda: py)
    cmd = installer.python_tools_install_command(installer.find("pytest"))
    assert cmd.startswith("'/Users/me/Library/Application Support/")
    assert _refuse_unsafe_install(cmd, py) == ""
