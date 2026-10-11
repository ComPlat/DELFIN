"""`!cmd` on the dashboard, through the same gate as on the terminal.

The terminal has had `!cmd` (repl_commands.SHELL_PREFIX) for a while: it
runs through `engine.run_gated_bash`, which is the agent's own dispatcher
and gate -- the same deny-list, secret scan, auto-allow list and approval
the model's bash calls go through. The dashboard had no `!` at all (goal
5 of the 2026-10-08 brief). One implementation now serves both, so the
gate cannot be forgotten on one of them.
"""

from __future__ import annotations

import json

from delfin.agent import repl as R
from delfin.agent.repl_commands import SHELL_PREFIX


class _Engine:
    def __init__(self, result):
        self._result = result
        self.calls = []

    def run_gated_bash(self, command):
        self.calls.append(command)
        if isinstance(self._result, Exception):
            raise self._result
        return self._result


def test_the_output_is_the_commands_not_the_envelope():
    eng = _Engine(json.dumps({"exit_code": 0, "stdout": "M README.md\n", "stderr": ""}))
    assert R.shell_escape(eng, "git status --short") == ["M README.md"]
    assert eng.calls == ["git status --short"]


def test_a_refusal_is_shown_as_one():
    eng = _Engine(json.dumps({"error": "blocked: rm -rf is on the deny-list"}))
    lines = R.shell_escape(eng, "rm -rf /")
    assert lines == ["refused: blocked: rm -rf is on the deny-list"]


def test_a_nonzero_exit_is_said():
    eng = _Engine(json.dumps({"exit_code": 2, "stdout": "", "stderr": "no such file"}))
    lines = R.shell_escape(eng, "ls nope")
    assert lines[0] == "no such file" and lines[-1] == "(exit 2)"


def test_a_backend_without_a_gate_has_no_shell():
    class _NoGate:
        pass
    lines = R.shell_escape(_NoGate(), "ls")
    assert len(lines) == 1 and "same gate" in lines[0]
    assert R.shell_escape(None, "ls")[0].startswith("! is not available")


def test_an_empty_command_does_nothing():
    eng = _Engine("{}")
    assert R.shell_escape(eng, "   ") == []
    assert eng.calls == []


def test_an_exception_is_shown_not_raised():
    lines = R.shell_escape(_Engine(RuntimeError("gate exploded")), "ls")
    assert lines == ["[error] gate exploded"]


def test_the_dashboard_is_wired_to_it():
    """The wiring only; the gate itself is `run_gated_bash`, tested where
    it lives. Source check because the tab is built around live widgets."""
    import inspect
    from delfin.dashboard import tab_agent

    src = inspect.getsource(tab_agent)
    assert "SHELL_PREFIX" in src and "shell_escape as _escape" in src
    assert "dashboard-shell-escape" in src, "it has to run off the UI thread"


def test_the_prefix_is_one_character_shared_by_both():
    assert SHELL_PREFIX == "!"
