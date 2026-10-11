"""The second opinion, driven through the real gate at the points review found.

Each test injects a stub runner on the permissions object; no test contacts
a model. Covered: the MCP side-effect gate and the gate's tail under Bypass,
the conversation bound before an MCP call, the per-turn lockout outside the
classified scope, the per-turn stall breaker, an "always allow" rule against
a chained command, git branch options that write, and a model refusal that
the user's next message lifts.
"""

import json
import time

import pytest

from delfin.agent import action_classifier as C
from delfin.agent import api_client as A
from delfin.agent import mcp_client as M


ALLOW = "DECISION: allow\nREASON: part of the task."
REFUSE = "DECISION: refuse\nREASON: nobody asked for this."
TASK = [{"role": "user", "content": "build the documentation"}]


class Runner:
    def __init__(self, *answers, delay=0.0):
        self.answers = list(answers)
        self.prompts = []
        self.delay = delay

    def __call__(self, prompt, model):
        self.prompts.append(prompt)
        if self.delay:
            time.sleep(self.delay)
        return self.answers.pop(0) if len(self.answers) > 1 else self.answers[0]


def _settings(broad=False, rules=()):
    return lambda: {"model": "m", "rules": list(rules), "broad": broad}


def _perms(tmp_path, mode, runner, messages=TASK):
    perms = A.KitToolPermissions(workspace=tmp_path, mode=mode)
    perms.classifier_runner = runner
    perms.classifier_messages = messages
    return perms


def _gate(perms, name, args):
    return A._DocToolExecutor()._run_permission_gate(name, args, perms)


# -- coverage under Bypass -------------------------------------------------

def test_a_side_effecting_mcp_tool_under_bypass_is_judged(tmp_path, monkeypatch):
    monkeypatch.setattr(A, "_classifier_settings", _settings())
    runner = Runner(REFUSE)
    perms = _perms(tmp_path, "bypassPermissions", runner)
    out = A._DocToolExecutor()._gate_mcp_tool(
        "mcp__github__create_pull_request", {"title": "x"}, perms)
    assert runner.prompts, "the classifier was never asked"
    assert out is not None and "classifier" in out


def test_run_tests_under_bypass_is_judged(tmp_path, monkeypatch):
    monkeypatch.setattr(A, "_classifier_settings", _settings())
    runner = Runner(REFUSE)
    perms = _perms(tmp_path, "bypassPermissions", runner)
    out = _gate(perms, "run_tests", {"target": "tests"})
    assert runner.prompts and out is not None


def test_run_tests_outside_bypass_in_the_default_scope_is_not(tmp_path,
                                                             monkeypatch):
    monkeypatch.setattr(A, "_classifier_settings", _settings())
    runner = Runner(REFUSE)
    perms = _perms(tmp_path, "acceptEdits", runner)
    assert _gate(perms, "run_tests", {"target": "tests"}) is None
    assert runner.prompts == []


# -- the conversation an MCP call is judged against -------------------------

class _Fn:
    def __init__(self, name=None, arguments=None):
        self.name, self.arguments = name, arguments


class _TC:
    def __init__(self, index, id=None, name=None, arguments=None):
        self.index, self.id, self.function = index, id, _Fn(name, arguments)


class _Delta:
    def __init__(self, content=None, tool_calls=None):
        self.content, self.tool_calls = content, tool_calls


class _Choice:
    def __init__(self, delta, finish=None):
        self.delta, self.finish_reason = delta, finish


class _Usage:
    prompt_tokens = 5
    completion_tokens = 3


class _Chunk:
    def __init__(self, choices, usage=None):
        self.choices, self.usage = choices, usage


class _Stream:
    def __init__(self, chunks):
        self._chunks = chunks

    def __iter__(self):
        return iter(self._chunks)

    def close(self):
        pass


class _Registry:
    def __init__(self):
        self.calls = []

    def discover_all(self): return []
    def discover_resources(self): return []
    def discover_prompts(self): return []

    def call(self, name, args):
        self.calls.append((name, args))
        return json.dumps({"status": "ok"})


def test_an_mcp_shell_call_is_judged_against_this_turns_request(
        tmp_path, monkeypatch):
    from delfin.agent import model_capabilities as mc
    caps = mc.ModelCapabilities(model="m", provider="ollama",
                                context_window=200_000, supports_tools=True)
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: caps)
    monkeypatch.setattr(A, "_classifier_settings", _settings())
    monkeypatch.setattr(M, "get_registry", lambda *a, **k: _Registry())
    client = A.create_client(backend="api", provider="ollama",
                             model="qwen2.5-coder:7b", cwd=str(tmp_path))
    runner = Runner(ALLOW)
    client._permissions = _perms(tmp_path, "bypassPermissions", runner,
                                 messages=None)
    streams = [
        _Stream([
            _Chunk([_Choice(_Delta(tool_calls=[_TC(
                0, id="c1", name="mcp__kit-coding__bash",
                arguments='{"command": "python3 make_docs.py"}')]))]),
            _Chunk([_Choice(_Delta(), finish="tool_calls")], usage=_Usage()),
        ]),
        _Stream([_Chunk([_Choice(_Delta(content="done"), finish="stop")],
                        usage=_Usage())]),
    ]
    seen = {"n": 0}

    def _create(**kwargs):
        s = streams[min(seen["n"], len(streams) - 1)]
        seen["n"] += 1
        return s

    client.client.chat.completions.create = _create
    list(client.stream_message(
        "sys", [{"role": "user", "content": "build the documentation"}],
        max_tokens=100))
    assert runner.prompts, "the classifier was never asked"
    assert "build the documentation" in runner.prompts[0]


# -- the lockout ------------------------------------------------------------

def test_after_three_refusals_an_auto_allowed_action_does_not_run(
        tmp_path, monkeypatch):
    monkeypatch.setattr(A, "_classifier_settings", _settings())
    runner = Runner(REFUSE)
    perms = _perms(tmp_path, "acceptEdits", runner)
    for n in range(3):
        assert _gate(perms, "web_fetch",
                     {"url": f"https://example.org/{n}"}) is not None
    assert perms.matches_bash_auto_allow("python3 build.py")
    calls = len(runner.prompts)
    out = _gate(perms, "bash", {"command": "python3 build.py"})
    assert out is not None and "three actions were refused" in out
    out = _gate(perms, "write_file", {"path": "a.txt", "content": "x"})
    assert out is not None and "three actions were refused" in out
    assert len(runner.prompts) == calls
    # Reading still runs.
    assert _gate(perms, "bash", {"command": "ls -la"}) is None


# -- a stalled endpoint ------------------------------------------------------

def test_a_stall_is_waited_for_once_per_turn(tmp_path, monkeypatch):
    monkeypatch.setattr(A, "_classifier_settings", _settings())
    monkeypatch.setattr(A, "_CLASSIFIER_TIMEOUT_S", 0.2)
    runner = Runner(ALLOW, delay=1.0)
    perms = _perms(tmp_path, "acceptEdits", runner)
    assert _gate(perms, "web_fetch", {"url": "https://example.org/a"}) is None
    t0 = time.monotonic()
    assert _gate(perms, "web_fetch", {"url": "https://example.org/b"}) is None
    assert time.monotonic() - t0 < 0.15
    assert len(runner.prompts) == 1
    # A new user message is a new turn: the model is asked again.
    perms.classifier_messages = TASK + [
        {"role": "assistant", "content": "ok"},
        {"role": "user", "content": "fetch the next page"}]
    _gate(perms, "web_fetch", {"url": "https://example.org/c"})
    assert len(runner.prompts) == 2


def test_a_stall_under_bypass_still_refuses_without_waiting(tmp_path,
                                                            monkeypatch):
    monkeypatch.setattr(A, "_classifier_settings", _settings())
    monkeypatch.setattr(A, "_CLASSIFIER_TIMEOUT_S", 0.2)
    runner = Runner(ALLOW, delay=1.0)
    perms = _perms(tmp_path, "bypassPermissions", runner)
    assert _gate(perms, "bash", {"command": "python3 a.py"}) is not None
    t0 = time.monotonic()
    assert _gate(perms, "bash", {"command": "python3 b.py"}) is not None
    assert time.monotonic() - t0 < 0.15
    assert len(runner.prompts) == 1


# -- user rules ---------------------------------------------------------------

@pytest.mark.parametrize("cmd", [
    "pytest -q && python3 upload.py",
    "pytest -q; python3 upload.py",
    "pytest_upload.sh",
    "pytest -q\npython3 upload.py",
])
def test_an_always_rule_does_not_cover_a_chained_or_longer_command(cmd):
    assert C.apply_user_rules(
        {"tool": "bash", "kind": "bash", "command": cmd},
        ["always allow: pytest"]) is None


def test_an_always_rule_still_covers_the_command_it_names():
    v = C.apply_user_rules(
        {"tool": "bash", "kind": "bash", "command": "pytest -q tests/"},
        ["always allow: pytest"])
    assert v is not None and v.decision == "allow"


def test_an_always_rule_reaches_the_model_through_the_gate(tmp_path,
                                                          monkeypatch):
    monkeypatch.setattr(A, "_classifier_settings",
                        _settings(rules=["always allow: pytest"]))
    runner = Runner(REFUSE)
    perms = _perms(tmp_path, "bypassPermissions", runner)
    assert _gate(perms, "bash",
                 {"command": "pytest -q && python3 upload.py"}) is not None
    assert runner.prompts


# -- the plain-read floor -------------------------------------------------------

@pytest.mark.parametrize("cmd", [
    "git branch --set-upstream-to=origin/main",
    "git branch --unset-upstream",
    "git branch --edit-description",
])
def test_a_git_branch_option_that_writes_is_not_a_plain_read(cmd):
    assert not C._bash_is_plain_read(cmd)


# -- a model refusal and the user's next word ------------------------------------

def test_the_users_next_message_lifts_a_model_refusal(tmp_path, monkeypatch):
    monkeypatch.setattr(A, "_classifier_settings", _settings())
    runner = Runner(REFUSE, ALLOW)
    perms = _perms(tmp_path, "bypassPermissions", runner)
    assert _gate(perms, "bash", {"command": "python3 deploy.py"}) is not None
    # Same turn: refused again, the model is not asked.
    assert _gate(perms, "bash", {"command": "python3 deploy.py"}) is not None
    assert len(runner.prompts) == 1
    perms.classifier_messages = TASK + [
        {"role": "assistant", "content": "May I run deploy.py?"},
        {"role": "user", "content": "yes, run deploy.py"}]
    assert _gate(perms, "bash", {"command": "python3 deploy.py"}) is None
    assert len(runner.prompts) == 2


def test_a_no_in_the_dialog_holds_past_the_users_next_message(tmp_path,
                                                             monkeypatch):
    monkeypatch.setattr(A, "_classifier_settings", _settings())
    runner = Runner("DECISION: ask\nREASON: unclear.")
    perms = _perms(tmp_path, "bypassPermissions", runner)
    perms.confirm_callback = lambda name, args, preview: False
    assert _gate(perms, "web_fetch", {"url": "https://example.org/x"}) is not None
    perms.classifier_messages = TASK + [
        {"role": "assistant", "content": "ok"},
        {"role": "user", "content": "go on"}]
    out = _gate(perms, "web_fetch", {"url": "https://example.org/x"})
    assert out is not None and "already refused" in out
    assert len(runner.prompts) == 1
