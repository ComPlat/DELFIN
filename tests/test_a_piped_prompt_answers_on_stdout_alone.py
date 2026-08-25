"""The headless contract: stdout is the answer, stderr is everything else.

`delfin-agent -p "..." > answer.txt` has to produce exactly the answer, so
nothing else may be written to stdout — not the banner, not tool lines,
not the session id. The JSON payload is pinned to its exact key set for
the same reason: it is a machine contract, and a key added or renamed
breaks callers silently.

The one-shot path deliberately routes through `cmd_run` rather than
reimplementing it, so these assertions cover both front doors at once.
"""

from __future__ import annotations

import io
import json
from pathlib import Path

import pytest

from delfin.agent import cli as agent_cli


_RESULT = {
    "text": "ANSWER",
    "tool_calls": [{"name": "bash", "input": {"command": "ls"}}],
    "input_tokens": 3,
    "output_tokens": 4,
    "error": "",
}


class _StubEngine:
    session_id = "sid-1"
    token_usage = {"input": 0, "output": 0}

    def stream_response(self, **kwargs):
        return "ANSWER"

    def export_state(self):
        return {}

    def record_cycle_outcome(self, *a, **k):
        return None


@pytest.fixture
def wired(monkeypatch, tmp_path):
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path))
    monkeypatch.setattr(agent_cli, "_build_engine", lambda args: _StubEngine())
    monkeypatch.setattr(agent_cli, "_run_once",
                        lambda engine, prompt, **kw: dict(_RESULT))
    monkeypatch.setattr(agent_cli, "_save_session",
                        lambda engine, root, **kw: "sid-1")
    monkeypatch.setattr(agent_cli, "_resume_or_create", lambda engine, args: "")


def _no_tty(monkeypatch, text: str = ""):
    stream = io.StringIO(text)
    stream.isatty = lambda: False           # type: ignore[method-assign]
    monkeypatch.setattr("sys.stdin", stream)


def test_a_printed_prompt_puts_only_the_answer_on_stdout(wired, monkeypatch, capsys):
    _no_tty(monkeypatch)
    rc = agent_cli.main(["-p", "do it"])
    out, err = capsys.readouterr()
    assert rc == 0
    assert out.strip() == "ANSWER"
    assert "sid-1" not in out
    assert "bash" not in out


def test_a_piped_prompt_is_read_from_stdin(wired, monkeypatch, capsys):
    _no_tty(monkeypatch, "do it from a pipe\n")
    seen: list[str] = []
    monkeypatch.setattr(
        agent_cli, "_run_once",
        lambda engine, prompt, **kw: (seen.append(prompt), dict(_RESULT))[1])
    rc = agent_cli.main([])
    out, _ = capsys.readouterr()
    assert rc == 0
    assert seen == ["do it from a pipe"]
    assert out.strip() == "ANSWER"


def test_the_json_payload_keeps_its_exact_shape(wired, monkeypatch, capsys):
    _no_tty(monkeypatch)
    rc = agent_cli.main(["-p", "do it", "--output-format", "json"])
    out, _ = capsys.readouterr()
    assert rc == 0
    payload = json.loads(out)
    assert set(payload) == {
        "text", "tool_calls", "input_tokens", "output_tokens",
        "error", "session_id",
    }
    assert payload["session_id"] == "sid-1"


def test_the_two_front_doors_report_the_same_thing(wired, monkeypatch, capsys):
    _no_tty(monkeypatch)
    agent_cli.main(["-p", "do it", "--output-format", "json"])
    via_chat = json.loads(capsys.readouterr().out)
    agent_cli.main(["run", "do", "it", "--json"])
    via_run = json.loads(capsys.readouterr().out)
    assert via_chat == via_run, (
        "chat -p and run diverged; the one-shot path is supposed to BE "
        "cmd_run, not a second implementation of it"
    )


def test_the_prompt_itself_never_lands_in_the_answer(monkeypatch):
    """`input(prompt)` writes the prompt to STDOUT.

    So an interactive session with a redirected stdout collected a `> ` for
    every turn, in the one stream that is supposed to carry the answer and
    nothing else. When stdout is a terminal the prompt still goes through
    input(), because readline needs it to redraw a wrapped line correctly.
    """
    from delfin.agent import repl

    out, err = io.StringIO(), io.StringIO()          # neither is a tty
    agent = repl.TerminalAgent(_StubEngine(), out=out, err=err)
    monkeypatch.setattr("builtins.input", lambda *a: "typed")

    assert agent._input("> ") == "typed"
    assert out.getvalue() == ""
    assert err.getvalue() == "> "


def test_one_turn_reports_the_keys_the_json_contract_promises(monkeypatch):
    """Pins `_run_once` itself, not just what main() does with it."""
    keys = agent_cli._run_once(_StubEngine(), "hello")
    assert set(keys) == {
        "text", "tool_calls", "input_tokens", "output_tokens", "error"}


def test_a_pipe_with_nothing_in_it_is_refused(wired, monkeypatch, capsys):
    _no_tty(monkeypatch, "")
    rc = agent_cli.main([])
    out, err = capsys.readouterr()
    assert rc == 2
    assert out == "", "the refusal belongs on stderr, stdout is the answer"
    assert "stdin" in err


def test_an_interactive_invocation_starts_the_session(wired, monkeypatch, capsys):
    """A terminal with nothing typed opens the loop and leaves on EOF.

    The banner belongs on stderr like every other piece of chrome, so a
    session that answers nothing still writes nothing to stdout.
    """
    stream = io.StringIO("")
    stream.isatty = lambda: True            # type: ignore[method-assign]
    monkeypatch.setattr("sys.stdin", stream)
    rc = agent_cli.main([])
    out, err = capsys.readouterr()
    assert rc == 0
    assert out == ""
    assert "delfin-agent" in err, "the banner says what you are looking at"
    assert "workspace" in err
