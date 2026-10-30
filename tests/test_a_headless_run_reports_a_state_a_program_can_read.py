"""The headless CLI as a state another program reads without guessing.

``delfin-agent chat -p ... --output-format json|stream-json`` is driven by
scripts and supervisors. What they need from one run, and what this file
pins:

* two stream events emitted where the gate decides, not reconstructed from
  text: ``permission_request`` (a call waiting for a human, under
  ``--supervised``) and ``permission_denied`` (a call the gate refused,
  with a stable ``reason_code``);
* a ``result`` object with ``schema_version``, ``stop_reason``,
  ``denials`` and ``files_changed``, one key list (``repl.TURN_KEYS``) for
  both formats;
* one exit code per stop reason;
* no credential in either new event.

The turns are real: the real engine, the real OpenAI-compatible tool loop
and the real permission gate, with only the model's HTTP stream scripted.
"""

from __future__ import annotations

import json
import threading
import time
from pathlib import Path

import pytest

from delfin.agent import api_client as _A
from delfin.agent import cli as agent_cli
from delfin.agent import mcp_client as M
from delfin.agent import turn_outcome as TO


@pytest.fixture(autouse=True)
def _uncaged(monkeypatch):
    """The refusal and approval events are driven by an off-list command;
    where the cage holds it, it runs unasked (_caged_command_may_run_unasked)
    and the verdict would depend on the host."""
    monkeypatch.setattr(_A, "_caged_command_may_run_unasked",
                        lambda *a, **k: False)


# --- scripted OpenAI-compatible stream -------------------------------------

class _FnDelta:
    def __init__(self, name=None, arguments=None):
        self.name = name
        self.arguments = arguments


class _TCDelta:
    def __init__(self, index, id=None, name=None, arguments=None):
        self.index = index
        self.id = id
        self.function = _FnDelta(name, arguments)


class _Delta:
    def __init__(self, content=None, tool_calls=None):
        self.content = content
        self.tool_calls = tool_calls


class _Choice:
    def __init__(self, delta, finish=None):
        self.delta = delta
        self.finish_reason = finish


class _Usage:
    prompt_tokens = 5
    completion_tokens = 3


class _Chunk:
    def __init__(self, choices, usage=None):
        self.choices = choices
        self.usage = usage


class _Stream:
    def __init__(self, chunks):
        self._chunks = chunks

    def __iter__(self):
        return iter(self._chunks)

    def close(self):
        pass


class _EmptyReg:
    def discover_all(self): return []
    def discover_resources(self): return []
    def discover_prompts(self): return []


def _call(name: str, args: dict, call_id: str = "c1") -> _Stream:
    return _Stream([_Chunk([_Choice(_Delta(tool_calls=[_TCDelta(
        0, id=call_id, name=name, arguments=json.dumps(args))]),
        finish="tool_calls")], usage=_Usage())])


def _answer(text: str = "done.") -> _Stream:
    return _Stream([_Chunk([_Choice(_Delta(content=text), finish="stop")],
                           usage=_Usage())])


# --- the headless run, end to end ------------------------------------------

@pytest.fixture
def headless(monkeypatch, tmp_path):
    """``run(argv, streams) -> (exit_code, stdout_lines, engine)``."""
    home = tmp_path / "home"
    ws = home / "ws"
    ws.mkdir(parents=True)
    # cmd_run moves into --cwd; monkeypatch puts the process back.
    monkeypatch.chdir(ws)
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))
    # The session-end outcome stage rewrites the packaged provider-profile
    # baseline for this provider (provider_profile.save_provider_profile);
    # a copy takes that write instead of the checkout.
    from delfin.agent import provider_profile as PP
    baseline = tmp_path / "learned_profiles.json"
    baseline.write_bytes(PP._DEFAULT_PATH.read_bytes())
    monkeypatch.setattr(PP, "_DEFAULT_PATH", baseline)

    from delfin.agent import model_capabilities as mc
    caps = mc.ModelCapabilities(model="m", provider="ollama",
                                context_window=200_000, supports_tools=True)
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: caps)
    monkeypatch.setattr(M, "get_registry", lambda *a, **k: _EmptyReg())
    monkeypatch.setattr(agent_cli, "_save_session",
                        lambda engine, root, **kw: "sid-1")
    monkeypatch.setattr(agent_cli, "_resume_or_create",
                        lambda engine, args: "")

    import io
    stdin = io.StringIO("")
    stdin.isatty = lambda: False           # type: ignore[method-assign]
    monkeypatch.setattr("sys.stdin", stdin)

    real_build = agent_cli._build_engine
    built: dict = {}

    def run(argv, streams, *, prepare=None):
        calls = {"n": 0}

        def _create(**kwargs):
            s = streams[min(calls["n"], len(streams) - 1)]
            calls["n"] += 1
            return s

        def _build(args):
            engine = real_build(args)
            engine.client.client.chat.completions.create = _create
            if prepare is not None:
                prepare(engine)
            built["engine"] = engine
            return engine

        monkeypatch.setattr(agent_cli, "_build_engine", _build)
        rc = agent_cli.main(["-p", "do the task", "--cwd", str(ws),
                             "--provider", "ollama", "--model", "m",
                             "--permission-mode", "default", *argv])
        return rc, built.get("engine")

    run.workspace = ws
    run.home = home
    return run


def _events(capsys) -> list[dict]:
    return [json.loads(line) for line in capsys.readouterr().out.splitlines()
            if line.strip()]


def _types(events) -> list[str]:
    return [e["type"] for e in events if e["type"] not in ("wait", "notice")]


def test_a_refused_call_is_an_event_with_a_reason_code(headless, capsys):
    """Refused, then recovered: the denial is reported, the run is done."""
    rc, _ = headless(["--output-format", "stream-json"], [
        _call("bash", {"command": "./build.sh"}),
        _call("write_file", {"path": "notes.txt", "content": "x\n"}, "c2"),
        _answer(),
    ])
    events = _events(capsys)
    assert _types(events) == [
        "tool_use", "permission_denied", "tool_use", "text", "result"]
    denied = next(e for e in events if e["type"] == "permission_denied")
    assert denied["tool"] == "bash"
    assert denied["reason_code"] == "not_allowlisted"
    assert denied["needs_approval"] is True
    assert "auto-allow" in denied["text"]
    result = events[-1]
    assert result["stop_reason"] == "done"
    assert result["denials"] == [
        {"tool": "bash", "reason_code": "not_allowlisted"}]
    assert result["files_changed"] == ["notes.txt"]
    assert result["schema_version"] == TO.SCHEMA_VERSION
    assert rc == 0


def test_a_turn_that_ends_on_an_approvable_refusal_needs_approval(
        headless, capsys):
    rc, _ = headless(["--output-format", "stream-json"], [
        _call("bash", {"command": "./build.sh"}),
        _answer("I need that command approved."),
    ])
    result = _events(capsys)[-1]
    assert result["stop_reason"] == "needs_approval"
    assert rc == TO.EXIT_CODES["needs_approval"] == 4


def test_a_turn_that_ends_on_a_hard_refusal_is_refused(headless, capsys):
    rc, _ = headless(["--output-format", "json"], [
        _call("bash", {"command": "cat ~/.ssh/id_rsa"}),
        _answer("That file is off limits."),
    ])
    out = _events(capsys)
    assert len(out) == 1, "json is one object"
    result = out[0]
    assert result["stop_reason"] == "refused"
    assert result["denials"] == [{"tool": "bash",
                                  "reason_code": "secret_path"}]
    assert rc == TO.EXIT_CODES["refused"] == 5


def test_the_tool_round_ceiling_is_max_turns(headless, capsys):
    rc, _ = headless(["--output-format", "json", "--max-turns", "1"], [
        _call("read_file", {"path": "missing.txt"}),
        _call("read_file", {"path": "missing.txt"}, "c2"),
    ])
    result = _events(capsys)[0]
    assert result["stop_reason"] == "max_turns"
    assert rc == TO.EXIT_CODES["max_turns"] == 6


def test_the_run_clock_is_timeout(headless, capsys):
    def _started_long_ago(engine):
        engine._run_started_at = time.time() - 1000.0

    rc, _ = headless(["--output-format", "json", "--max-run-seconds", "10"],
                     [_answer()], prepare=_started_long_ago)
    result = _events(capsys)[0]
    assert result["stop_reason"] == "timeout"
    assert rc == TO.EXIT_CODES["timeout"] == 8


def test_a_clean_turn_is_done_and_exits_zero(headless, capsys):
    rc, _ = headless(["--output-format", "json"], [_answer("hello")])
    result = _events(capsys)[0]
    assert result["stop_reason"] == "done"
    assert result["denials"] == []
    assert result["files_changed"] == []
    assert rc == 0


def test_a_secret_in_a_refused_command_never_reaches_the_event(
        headless, capsys):
    secret = "sk-live-abcdef0123456789"
    headless(["--output-format", "stream-json"], [
        _call("bash", {"command": f"TOKEN={secret} python3 -c 1"}),
        _answer(),
    ])
    for event in _events(capsys):
        if event["type"] in ("permission_denied", "permission_request"):
            assert secret not in json.dumps(event)
            assert "<redacted" in json.dumps(event)


def test_without_supervised_nobody_is_asked(headless, capsys):
    _, engine = headless(["--output-format", "stream-json"], [
        _call("bash", {"command": "./build.sh"}), _answer()])
    assert engine.kit_permissions.confirm_callback is None
    assert "permission_request" not in _types(_events(capsys))


def test_under_supervised_the_request_is_an_event_a_human_answers(
        headless, capsys):
    """The request reaches the stream with the approvals-directory id; the
    answer here is a refusal written the way `delfin-agent approvals` does."""
    from delfin.agent import file_confirm as FC

    secret = "sk-live-abcdef0123456789"
    stop = threading.Event()

    def _deny_when_asked():
        while not stop.is_set():
            for req in FC.pending():
                FC.answer(req["id"], False)
            time.sleep(0.05)

    answerer = threading.Thread(target=_deny_when_asked, daemon=True)
    answerer.start()
    try:
        rc, _ = headless(["--output-format", "stream-json", "--supervised"], [
            _call("bash", {"command": f"TOKEN={secret} ./build.sh"}),
            _answer("Not approved."),
        ])
    finally:
        stop.set()
        answerer.join(2)
    events = _events(capsys)
    assert _types(events) == [
        "tool_use", "permission_request", "permission_denied", "text",
        "result"]
    request = next(e for e in events if e["type"] == "permission_request")
    assert set(request) == TO.PERMISSION_REQUEST_KEYS
    assert request["tool"] == "bash"
    assert request["answered_by"] == "human"
    assert request["approvals_dir"] is True
    assert request["mode"] == "default"
    assert request["request_id"]
    assert "./build.sh" in request["shown"]
    assert secret not in json.dumps(request)
    denied = next(e for e in events if e["type"] == "permission_denied")
    assert denied["reason_code"] == "denied_by_user"
    assert events[-1]["stop_reason"] == "refused"
    assert rc == 5


# --- the schema, pinned -----------------------------------------------------

def test_the_result_keys_are_pinned():
    from delfin.agent.repl import TURN_KEYS
    assert TURN_KEYS == frozenset({
        "text", "tool_calls", "input_tokens", "output_tokens", "error",
        "schema_version", "stop_reason", "denials", "files_changed",
        "model_switches"})
    assert TO.SCHEMA_VERSION == 1


def test_the_stop_reasons_and_exit_codes_are_pinned():
    assert TO.EXIT_CODES == {
        "done": 0, "error": 1, "needs_approval": 4, "refused": 5,
        "max_turns": 6, "budget": 7, "timeout": 8}
    assert set(TO.STOP_REASONS) == set(TO.EXIT_CODES)
    # 2 (usage) and 3 (engine init failed) keep the meaning they had
    # before this schema: no run happened, so there is no result.
    assert not {2, 3} & set(TO.EXIT_CODES.values())


def test_the_event_keys_are_pinned():
    assert TO.PERMISSION_REQUEST_KEYS == frozenset({
        "type", "request_id", "tool", "shown", "preview", "answered_by",
        "mode", "approvals_dir"})
    assert TO.PERMISSION_DENIED_KEYS == frozenset({
        "type", "tool", "reason_code", "needs_approval", "text"})


def test_every_reason_code_is_a_gate_category():
    """The codes are the panel's kinds plus the refusals the panel does not
    record; nothing else may appear in a denial."""
    from delfin.agent import security_events
    assert TO.LOCAL_REASON_CODES.isdisjoint(security_events.known_kinds())
    assert TO.REASON_CODES == (security_events.known_kinds()
                               | TO.LOCAL_REASON_CODES)
    assert TO.APPROVABLE_REASON_CODES <= TO.REASON_CODES


def test_engine_stop_reasons_map_to_the_documented_set():
    cases = {
        "end_turn": "done", "": "done", "end_turn_open_tasks": "done",
        "max_tool_rounds": "max_turns", "cost_cap": "budget",
        "max_turn_output": "budget", "run_budget_usd": "budget",
        "run_budget_turns": "budget", "run_budget_s": "timeout",
        "error": "error", "no_progress_loop": "error",
        "consecutive_identical_errors": "error", "max_context": "error",
        "stopped": "error",
    }
    for engine_reason, expected in cases.items():
        assert TO.classify(error="", engine_stop_reason=engine_reason,
                           denials=[], ended_on_denial=False) == expected
    assert TO.classify(error="boom", engine_stop_reason="end_turn",
                       denials=[], ended_on_denial=False) == "error"


def test_the_interactive_turn_carries_the_same_keys():
    from delfin.agent import repl

    class _Stub:
        token_usage = {"input": 0, "output": 0}

        def stream_response(self, **kwargs):
            return "x"

    headless_keys = set(agent_cli._run_once(_Stub(), "hi"))
    interactive = set(repl.TurnResult().to_dict())
    assert headless_keys == interactive == set(repl.TURN_KEYS)
