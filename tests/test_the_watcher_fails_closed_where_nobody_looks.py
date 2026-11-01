"""The model second opinion (the watcher) when its model fails.

Input: actions the permission rules already allowed, judged through the
real gate (``_DocToolExecutor._run_permission_gate`` / ``execute``) with
an injected classifier runner that raises, hangs, garbles or answers, a
recording confirm callback (attended) or none (unattended), in the modes
default, acceptEdits and bypassPermissions and both scopes. A fake clock
replaces ``time.monotonic`` for the probe interval. No test contacts a
model.

Output, pinned here:

* a failure on an action that LEAVES the workspace (network, a message to
  another agent, a push or upload, a write outside the workspace or into
  a file that instructs an agent, under Bypass a shell command outside
  the default auto-allow list) asks when someone can be asked and is
  refused when nobody can; a failure on a workspace-local action keeps the
  rule result;
* one short retry on a model error before it counts as a failure;
* after ``_CLASSIFIER_BREAKER_FAILURES`` consecutive failures the watcher
  is down: one notice, workspace work continues without model calls,
  leaving actions are held, the model is probed at most once per
  ``_CLASSIFIER_PROBE_INTERVAL_S``, and the first answer closes it again;
* a watcher refusal never counts towards the three-refusals lockout;
* in an unattended session an auto-allowed command that leaves the
  workspace (``git fetch``) is judged in the narrow scope too;
* ``DELFIN_CLASSIFIER=off`` is ignored in a process an agent session
  started (``DELFIN_AGENT_SESSION``), honoured elsewhere;
* the state (on/off, model, scope, last failure) is in the start banner
  and ``/status``;
* on the success path the model is called exactly as often as before.
"""

from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

import pytest

from delfin.agent import action_classifier as C
from delfin.agent import api_client as A


ALLOW = "DECISION: allow\nREASON: part of the task."
TASK = [{"role": "user", "content": "rebuild the report in this workspace"}]
ROOT = Path(__file__).resolve().parents[1]


class Clock:
    def __init__(self):
        self.t = 1000.0

    def __call__(self):
        return self.t


class Runner:
    """``answer`` None raises; ``fail`` > 0 raises that many times first."""

    def __init__(self, answer=None, *, fail=0, delay=0.0):
        self.answer = answer
        self.fail = fail
        self.delay = delay
        self.prompts = []

    def __call__(self, prompt, model):
        self.prompts.append(prompt)
        if self.delay:
            time.sleep(self.delay)
        if self.fail > 0:
            self.fail -= 1
            raise ConnectionError("endpoint refused")
        if self.answer is None:
            raise ConnectionError("endpoint refused")
        return self.answer


class Dialog:
    last_timed_out = False

    def __init__(self, answer=True):
        self.answer = answer
        self.shown = []

    def ask(self, name, args, preview):
        self.shown.append((name, preview))
        return self.answer


@pytest.fixture(autouse=True)
def _env(monkeypatch):
    monkeypatch.setattr(A, "_CLASSIFIER_RETRY_DELAY_S", 0.0, raising=False)
    clock = Clock()
    monkeypatch.setattr(A, "_classifier_clock", clock, raising=False)
    monkeypatch.setattr(
        A, "_classifier_settings",
        lambda: {"model": "m", "rules": [], "broad": False})
    return clock


def _scope(monkeypatch, broad):
    monkeypatch.setattr(
        A, "_classifier_settings",
        lambda: {"model": "m", "rules": [], "broad": broad})


def _perms(tmp_path, mode, runner, dialog=None):
    ws = tmp_path / "ws"
    ws.mkdir(exist_ok=True)
    perms = A.KitToolPermissions(workspace=ws.resolve(), mode=mode)
    perms.classifier_runner = runner
    perms.classifier_messages = TASK
    if dialog is not None:
        perms.confirm_callback = dialog.ask
    return perms


def _gate(perms, name, args):
    return A._DocToolExecutor()._run_permission_gate(name, args, perms)


def _second(perms, name, args):
    return A._DocToolExecutor()._classifier_second_opinion(
        perms, name, args, perms.mode)


MODES = ["default", "acceptEdits", "bypassPermissions"]


# -- the failure policy ------------------------------------------------------

@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("attended", [True, False])
@pytest.mark.parametrize("failure", ["error", "timeout", "garbled"])
def test_a_failure_on_a_network_call_asks_or_refuses(
        tmp_path, monkeypatch, mode, attended, failure):
    monkeypatch.setattr(A, "_CLASSIFIER_TIMEOUT_S", 0.1)
    runner = {"error": Runner(),
              "timeout": Runner(ALLOW, delay=0.5),
              "garbled": Runner("no idea")}[failure]
    dialog = Dialog(False) if attended else None
    perms = _perms(tmp_path, mode, runner, dialog)
    out = _gate(perms, "web_fetch", {"url": "https://example.org/a"})
    assert out is not None
    if attended:
        assert len(dialog.shown) == 1
        assert "watcher" in dialog.shown[0][1].lower()
    else:
        assert "watcher" in out.lower()


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("attended", [True, False])
def test_a_failure_on_a_message_to_another_agent_asks_or_refuses(
        tmp_path, mode, attended):
    dialog = Dialog(False) if attended else None
    perms = _perms(tmp_path, mode, Runner(), dialog)
    out = _second(perms, "session_message",
                  {"to": "ec9c6393", "message": "status: done, your turn"})
    if mode == "bypassPermissions":
        # Bypass screens a sibling message without the model
        # (test_bypass_asks_only_where_a_decision_is_needed): no model
        # call, so no failure to handle.
        assert out is None
        assert dialog is None or dialog.shown == []
        return
    assert out is not None
    if attended:
        assert len(dialog.shown) == 1


@pytest.mark.parametrize("cmd", [
    "make install-data",
    "git clone https://example.org/x/y.git reports/y",
])
@pytest.mark.parametrize("attended", [True, False])
def test_under_bypass_a_command_outside_the_allow_list_fails_closed(
        tmp_path, monkeypatch, attended, cmd):
    # A host without a bash cage; under a cage plain workspace work runs
    # free (test_bypass_caged_workspace_work_runs_free.py), a fetch never.
    monkeypatch.setattr(A, "_bash_runs_caged", lambda perms: False)
    dialog = Dialog(False) if attended else None
    runner = Runner()
    perms = _perms(tmp_path, "bypassPermissions", runner, dialog)
    out = _gate(perms, "bash", {"command": cmd})
    assert out is not None
    if attended:
        assert len(dialog.shown) == 1
    else:
        assert "watcher" in out.lower()
        assert "inside the workspace" in out


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("attended", [True, False])
def test_a_failure_on_workspace_work_keeps_the_rule_result(
        tmp_path, monkeypatch, mode, attended):
    _scope(monkeypatch, True)
    dialog = Dialog(False) if attended else None
    perms = _perms(tmp_path, mode, Runner(), dialog)
    if mode != "default":
        assert _gate(perms, "write_file",
                     {"path": "report.md", "content": "x\n"}) is None
    assert _gate(perms, "bash", {"command": "python3 build.py"}) is None
    assert _gate(perms, "bash", {"command": "pytest -q"}) is None
    assert dialog is None or dialog.shown == []


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("attended", [True, False])
def test_a_plain_read_never_calls_the_model(tmp_path, mode, attended):
    runner = Runner()
    dialog = Dialog(False) if attended else None
    perms = _perms(tmp_path, mode, runner, dialog)
    for cmd in ("ls -la", "git status", "cat notes.md"):
        assert _gate(perms, "bash", {"command": cmd}) is None
    assert _gate(perms, "read_file", {"path": "notes.md"}) is None
    assert runner.prompts == []


def test_one_error_is_retried_once_before_it_counts(tmp_path):
    runner = Runner(ALLOW, fail=1)
    perms = _perms(tmp_path, "acceptEdits", runner)
    assert _gate(perms, "web_fetch", {"url": "https://example.org/a"}) is None
    assert len(runner.prompts) == 2
    assert perms.classifier_health.get("failures", 0) == 0


def test_a_timeout_is_not_retried(tmp_path, monkeypatch):
    monkeypatch.setattr(A, "_CLASSIFIER_TIMEOUT_S", 0.1)
    runner = Runner(ALLOW, delay=0.5)
    perms = _perms(tmp_path, "acceptEdits", runner)
    assert _gate(perms, "web_fetch", {"url": "https://example.org/a"})
    assert len(runner.prompts) == 1


def test_a_watcher_refusal_does_not_lock_the_turn(tmp_path):
    perms = _perms(tmp_path, "bypassPermissions", Runner())
    for i in range(5):
        assert _gate(perms, "web_fetch", {"url": f"https://example.org/{i}"})
    assert _gate(perms, "bash", {"command": "python3 build.py"}) is None
    assert _gate(perms, "write_file",
                 {"path": "report.md", "content": "x\n"}) is None


# -- the circuit breaker -----------------------------------------------------

def test_the_breaker_opens_holds_and_closes(tmp_path, monkeypatch, _env):
    # `make` is the leaving probe here, as on a host without a bash cage.
    monkeypatch.setattr(A, "_bash_runs_caged", lambda perms: False)
    clock = _env
    runner = Runner()
    perms = _perms(tmp_path, "bypassPermissions", runner)
    for i in range(A._CLASSIFIER_BREAKER_FAILURES):
        assert _gate(perms, "web_fetch", {"url": f"https://example.org/{i}"})
        clock.t += 1
    health = perms.classifier_health
    assert health["down"] is True
    notes = A.drain_watcher_notices(perms)
    assert len(notes) == 1 and "down" in notes[0]
    assert A.drain_watcher_notices(perms) == []
    calls = len(runner.prompts)
    # Down: leaving actions are held without a model call; workspace work
    # runs.
    for i in range(10):
        out = _gate(perms, "bash", {"command": f"make target{i}"})
        assert out is not None and "watcher" in out.lower()
        assert _gate(perms, "write_file",
                     {"path": f"r{i}.md", "content": "x\n"}) is None
        assert _gate(perms, "bash", {"command": f"python3 s{i}.py"}) is None
        clock.t += 1
    assert len(runner.prompts) == calls
    # A probe after the interval; the model answers: closed again.
    runner.answer = ALLOW
    clock.t += A._CLASSIFIER_PROBE_INTERVAL_S
    assert _gate(perms, "bash", {"command": "make final"}) is None
    assert perms.classifier_health["down"] is False
    notes = A.drain_watcher_notices(perms)
    assert len(notes) == 1 and "again" in notes[0]


def test_a_probe_runs_at_most_once_a_minute(tmp_path, _env):
    clock = _env
    runner = Runner()
    perms = _perms(tmp_path, "bypassPermissions", runner)
    for i in range(A._CLASSIFIER_BREAKER_FAILURES):
        _gate(perms, "web_fetch", {"url": f"https://example.org/{i}"})
    calls = len(runner.prompts)
    for i in range(240):                 # four minutes, one action a second
        _gate(perms, "web_fetch", {"url": f"https://example.org/x{i}"})
        clock.t += 1
    probes = len(runner.prompts) - calls
    # One probe per minute, each retried once on its error.
    assert probes <= 2 * (240 // int(A._CLASSIFIER_PROBE_INTERVAL_S) + 1)


def test_a_long_unattended_run_survives_an_outage(tmp_path, _env):
    """A headless Bypass run of 120 steps, the watcher down for steps
    20-79. Every workspace step succeeds throughout; leaving steps are held
    while it is down and judged again once it answers."""
    clock = _env
    runner = Runner(ALLOW)
    perms = _perms(tmp_path, "bypassPermissions", runner)
    held, judged_after, workspace_fail = 0, 0, 0
    for step in range(120):
        runner.answer = None if 20 <= step < 80 else ALLOW
        before = len(runner.prompts)
        for name, args in (
                ("write_file", {"path": f"out/{step}.md", "content": "x\n"}),
                ("bash", {"command": f"python3 analyse.py {step}"}),
                ("bash", {"command": f"ls out/{step}.md"})):
            if _gate(perms, name, args) is not None:
                workspace_fail += 1
        if step % 10 == 5:
            out = _gate(perms, "web_fetch",
                        {"url": f"https://example.org/{step}"})
            if out is not None:
                held += 1
            elif step >= 80 and len(runner.prompts) > before:
                judged_after += 1
        clock.t += 5
    assert workspace_fail == 0
    assert held == 6                     # steps 25, 35, ..., 75
    assert judged_after == 4             # steps 85, 95, 105, 115
    assert perms.classifier_health["down"] is False


# -- scope: an unattended session ---------------------------------------------

def test_an_unattended_auto_allowed_fetch_is_judged(tmp_path):
    runner = Runner(ALLOW)
    perms = _perms(tmp_path, "default", runner)
    assert _gate(perms, "bash", {"command": "git fetch origin"}) is None
    assert len(runner.prompts) == 1


def test_an_attended_auto_allowed_fetch_is_not(tmp_path):
    runner = Runner(ALLOW)
    perms = _perms(tmp_path, "default", runner, Dialog())
    assert _gate(perms, "bash", {"command": "git fetch origin"}) is None
    assert runner.prompts == []


# -- the latency cost on the success path --------------------------------------

SUCCESS_PATH = [
    ("bash", {"command": "ls -la"}),
    ("bash", {"command": "python3 build.py"}),
    ("bash", {"command": "pytest -q"}),
    ("write_file", {"path": "report.md", "content": "x\n"}),
    ("web_fetch", {"url": "https://example.org/a"}),
]


@pytest.mark.parametrize("mode,attended,expected", [
    ("default", True, 1), ("default", False, 1),
    ("acceptEdits", True, 1), ("acceptEdits", False, 1),
    ("bypassPermissions", True, 3), ("bypassPermissions", False, 3),
])
def test_the_success_path_calls_the_model_as_often_as_before(
        tmp_path, mode, attended, expected):
    """Model calls per action, narrow scope, the model answering allow:
    default/acceptEdits judge the network call only (1 of 5); Bypass also
    the two commands (3 of 5; the workspace write is screened)."""
    runner = Runner(ALLOW)
    perms = _perms(tmp_path, mode, runner, Dialog() if attended else None)
    for name, args in SUCCESS_PATH:
        _gate(perms, name, args)
    assert len(runner.prompts) == expected


# -- switching it off ----------------------------------------------------------

def _real_settings(monkeypatch, tmp_path, payload=None):
    monkeypatch.undo()
    import delfin.user_settings as US
    monkeypatch.setattr(US, "load_settings",
                        lambda *a, **k: payload or {})


def test_off_is_honoured_outside_an_agent_session(monkeypatch, tmp_path):
    _real_settings(monkeypatch, tmp_path)
    monkeypatch.delenv(C.AGENT_SESSION_ENV, raising=False)
    monkeypatch.setenv("DELFIN_CLASSIFIER", "off")
    assert A._classifier_settings() is None


def test_off_is_ignored_inside_an_agent_session(monkeypatch, tmp_path):
    _real_settings(monkeypatch, tmp_path)
    monkeypatch.setenv(C.AGENT_SESSION_ENV, "4242")
    monkeypatch.setenv("DELFIN_CLASSIFIER", "off")
    assert A._classifier_settings() is not None


def test_the_users_settings_switch_still_works(monkeypatch, tmp_path):
    _real_settings(monkeypatch, tmp_path,
                   {"agent": {"classifier": {"enabled": False}}})
    monkeypatch.delenv("DELFIN_CLASSIFIER", raising=False)
    assert A._classifier_settings() is None


def test_every_agent_child_carries_the_session_marker():
    env = A._scrubbed_bash_env()
    assert env.get(C.AGENT_SESSION_ENV)


def test_a_child_delfin_started_by_an_agent_keeps_the_watcher(tmp_path):
    from tests.conftest import child_env
    env = child_env(tmp_path)
    env.update({k: v for k, v in A._scrubbed_bash_env().items()
                if k == C.AGENT_SESSION_ENV})
    env["DELFIN_CLASSIFIER"] = "off"
    env.pop("PYTHONPATH", None)
    code = ("from delfin.agent import api_client as A; "
            "print(A._classifier_settings() is not None)")
    done = subprocess.run([sys.executable, "-c", code], cwd=str(ROOT),
                          env=env, capture_output=True, text=True,
                          timeout=300)
    assert done.stdout.strip().splitlines()[-1] == "True", done.stderr
    env.pop(C.AGENT_SESSION_ENV)
    done = subprocess.run([sys.executable, "-c", code], cwd=str(ROOT),
                          env=env, capture_output=True, text=True,
                          timeout=300)
    assert done.stdout.strip().splitlines()[-1] == "False", done.stderr


@pytest.mark.parametrize("mode", MODES)
def test_the_settings_file_is_a_protected_path(tmp_path, mode):
    perms = _perms(tmp_path, mode, Runner(ALLOW))
    assert perms.matches_path_protected(".delfin_settings.json")
    out = _gate(perms, "write_file",
                {"path": ".delfin_settings.json",
                 "content": '{"agent": {"classifier": {"enabled": false}}}'})
    assert out is not None
    dialog = Dialog(False)
    perms = _perms(tmp_path, mode, Runner(ALLOW), dialog)
    _gate(perms, "write_file",
          {"path": ".delfin_settings.json", "content": "{}"})
    assert len(dialog.shown) == 1


# -- what the user sees --------------------------------------------------------

def test_the_status_line_names_switch_model_scope_and_failure(tmp_path):
    perms = _perms(tmp_path, "bypassPermissions", Runner())
    line = A.watcher_status_line(perms)
    assert "on" in line and "m" in line and "unattended" in line
    assert "no failure" in line
    _gate(perms, "web_fetch", {"url": "https://example.org/a"})
    line = A.watcher_status_line(perms)
    assert "ConnectionError" in line


def test_the_status_line_says_off(tmp_path, monkeypatch):
    monkeypatch.setattr(A, "_classifier_settings", lambda: None)
    perms = _perms(tmp_path, "default", Runner())
    assert "off" in A.watcher_status_line(perms)
    perms.classifier_runner = None
    assert "off" in A.watcher_status_line(perms)


class _Git:
    is_repo = False
    branch = ""
    dirty = ()
    unborn = False


class _Report:
    git = _Git()
    granted_dirs = ()
    read_dirs = ()


def test_the_start_banner_and_status_show_the_watcher(tmp_path):
    from delfin.agent import cli as agent_cli
    from delfin.agent import repl_commands as rc

    perms = _perms(tmp_path, "bypassPermissions", Runner(ALLOW))

    class _Engine:
        client = type("C", (), {"model": "m"})()
        provider = "kit"
        mode = "solo"
        session_id = "0" * 32
        kit_permissions = perms

        def get_status(self):
            return {"watcher": A.watcher_status_line(perms)}

    text = agent_cli._startup_banner(_Engine(), _Report(), Path("/tmp/p"))
    assert "watcher" in text
    out = rc.dispatch("/status", rc.ReplContext(engine=_Engine())).output
    assert "watcher" in out


@pytest.mark.parametrize("attended", [True, False])
def test_under_a_caged_bypass_a_fetch_still_fails_closed(
        tmp_path, monkeypatch, attended):
    monkeypatch.setattr(A, "_bash_runs_caged", lambda perms: True)
    dialog = Dialog(False) if attended else None
    perms = _perms(tmp_path, "bypassPermissions", Runner(), dialog)
    out = _gate(perms, "bash", {
        "command": "git clone https://example.org/x/y.git reports/y"})
    assert out is not None
