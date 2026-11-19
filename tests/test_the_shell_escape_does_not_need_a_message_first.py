"""`!cmd` on the agent tab: before the first message, and on any backend.

Two defects, one shape -- the shell escape reached for a gate that was
not where it looked.

1. The gate is reached as ``engine.run_gated_bash``, so `!` needs an
   engine. On the terminal there always is one: the session is built
   before the input loop reads a line. On the agent tab the engine is
   built by the FIRST MESSAGE (``_ensure_engine``), and the shell-escape
   branch returns before reaching it. So a `!` typed into a freshly
   opened tab read ``state["engine"]`` -- initialised to None -- and was
   answered "! is not available on this backend ... and this one has
   none". That names the wrong cause twice: the backend did carry a
   gate, and nothing about it had to change. A person following the
   message goes looking at their provider configuration.

2. ``kit_permissions`` reads ``client._permissions``, which only the
   OpenAI-compatible client carries, so `!` was unavailable on the CLI
   backend -- the factory default -- although that backend gates the
   model's own tools through ``cli_gate``. See
   ``test_the_shell_escape_runs_on_the_cli_backend_too``.

Measured on 2026-11-19. Both branches had passed their tests, because
the only dashboard test for `!` asserted that two strings occur in the
module: the logic lived in a closure inside ``create_tab``, built around
live widgets, where nothing could drive it. The fix is the two
module-level functions exercised here -- the closure now holds the
widget plumbing and nothing else.
"""

from __future__ import annotations

import ast
import inspect
import json

from delfin.dashboard import tab_agent


class _Engine:
    """An engine whose gate answers; the gate itself is tested where it
    lives (``repl.shell_escape_result`` / ``engine.run_gated_bash``)."""

    def __init__(self, payload, workspace=None):
        self._payload = payload
        self.commands = []
        self.kit_permissions = workspace and _Perms(workspace)

    def run_gated_bash(self, command):
        self.commands.append(command)
        return json.dumps(self._payload)


class _Perms:
    def __init__(self, directory):
        self.agent_workspace_dir = str(directory)


# --------------------------------------------------------------------- 1

def test_an_absent_engine_is_built_rather_than_reported_as_a_backend():
    built = []

    def build():
        built.append(True)
        return "the new engine"

    assert tab_agent.shell_escape_engine(None, build) == "the new engine"
    assert built, (
        "a `!` sent before the first message must build the session the "
        "way an ordinary message does; otherwise it is refused with a "
        "sentence about the backend")


def test_an_engine_that_exists_is_not_rebuilt():
    """A second `!` must not start a second session."""
    def build():
        raise AssertionError("must not be called when an engine is in hand")

    assert tab_agent.shell_escape_engine("the session's engine", build) == (
        "the session's engine")


def test_a_build_that_fails_yields_no_engine_and_does_not_raise():
    """Into the UI thread, a raise would lose the typed command without a
    word. The caller renders the refusal instead."""
    def build():
        raise RuntimeError("no backend configured")

    assert tab_agent.shell_escape_engine(None, build) is None


# --------------------------------------------------------- the entry made

def test_the_output_of_the_command_reaches_the_chat(tmp_path):
    eng = _Engine({"exit_code": 0, "stdout": "hello\n", "stderr": ""})

    entry = tab_agent.shell_escape_entry(eng, "echo hello", tmp_path)

    assert eng.commands == ["echo hello"]
    assert entry["role"] == "system"
    assert "! echo hello" in entry["content"]
    assert "hello" in entry["content"]
    assert entry["shell"]["exit_code"] == 0


def test_a_refusal_is_shown_as_a_refusal(tmp_path):
    eng = _Engine({"error": "blocked: rm -rf is on the deny-list"})

    entry = tab_agent.shell_escape_entry(eng, "rm -rf /", tmp_path)

    assert "refused: blocked: rm -rf is on the deny-list" in entry["content"]


def test_a_nonzero_exit_is_said_and_stderr_kept(tmp_path):
    eng = _Engine({"exit_code": 2, "stdout": "", "stderr": "no such file"})

    entry = tab_agent.shell_escape_entry(eng, "ls nope", tmp_path)

    assert "no such file" in entry["content"]
    assert "(exit 2)" in entry["content"]


def test_without_an_engine_the_person_is_told_not_crashed(tmp_path):
    """The path reached when ``shell_escape_engine`` could build none."""
    entry = tab_agent.shell_escape_entry(None, "echo hi", tmp_path)

    assert "[error]" in entry["content"]
    assert entry["shell"]["exit_code"] is None


def test_output_too_long_for_the_chat_is_saved_whole(tmp_path):
    """The cap belongs to the chat, not to the command: the whole stream
    is kept in a file and the entry names it."""
    big = "x" * (tab_agent._SHELL_CHAT_CHARS + 500)
    eng = _Engine({"exit_code": 0, "stdout": big, "stderr": ""},
                  workspace=tmp_path)

    entry = tab_agent.shell_escape_entry(eng, "seq 1 100000", tmp_path)

    saved = entry["shell"]["saved"]["stdout"]
    assert saved.startswith(str(tmp_path)), saved
    assert len(open(saved, encoding="utf-8").read()) == len(big)
    assert len(entry["shell"]["stdout"]) == tab_agent._SHELL_CHAT_CHARS


def test_the_engines_own_workspace_wins_over_the_fallback(tmp_path):
    """Where long output lands is the session's workspace when it has
    one; the fallback is for a session that names none."""
    mine, other = tmp_path / "mine", tmp_path / "other"
    big = "y" * (tab_agent._SHELL_CHAT_CHARS + 10)
    eng = _Engine({"exit_code": 0, "stdout": big, "stderr": ""},
                  workspace=mine)

    entry = tab_agent.shell_escape_entry(eng, "cat big", other)

    assert str(mine) in entry["shell"]["saved"]["stdout"]


# ------------------------------------------------------- the wiring only

def _shell_branch_source() -> str:
    tree = ast.parse(inspect.getsource(tab_agent))
    send = next(n for n in ast.walk(tree)
                if isinstance(n, ast.FunctionDef) and n.name == "_on_send")
    for node in ast.walk(send):
        if (isinstance(node, ast.If)
                and "SHELL_PREFIX" in ast.unparse(node.test)
                and "startswith" in ast.unparse(node.test)):
            return "\n".join(ast.unparse(s) for s in node.body)
    raise AssertionError(
        "_on_send has no branch keyed on SHELL_PREFIX; the shell escape "
        "is gone from this surface")


def test_the_closure_calls_the_two_functions_tested_above():
    """All this check has to carry now: the branch is four statements of
    widget plumbing around them. What it USED to carry -- the engine, the
    gate call, the saving, the rendering -- is driven by the tests above."""
    body = _shell_branch_source()

    assert "shell_escape_engine(state.get('engine'), _ensure_engine)" in body
    assert "shell_escape_entry(" in body
    assert "threading.Thread" in body and "dashboard-shell-escape" in body, (
        "the gate may ask the person and the command may take a while; "
        "neither may happen on the UI thread")


def test_the_prefix_is_the_terminals(tmp_path):
    from delfin.agent.repl_commands import SHELL_PREFIX

    assert SHELL_PREFIX == "!"
