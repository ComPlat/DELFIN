"""A headless plan is approvable, and an unanswered one says so.

``delfin-agent chat -p --permission-mode plan`` is how a script or a
supervisor runs a planning turn. What that caller is entitled to know, and
what this file pins:

* with ``--supervised`` -- somebody watching the approvals directory -- the
  plan reaches that directory, and an approval there lets the plan execute.
  Before this the run was told "no approval channel is available in this
  context", so its plan could never be approved and so could never run
  (measured on the base commit: 0 of 4 supervised plan runs produced a
  question in the directory);
* a turn that ends with its plan still unanswered is not a finished turn:
  ``stop_reason`` is ``awaiting_approval``, exit code 11;
* a plan the person refused ends the turn as a refusal (exit 5), not as
  success -- the requested change will not be made;
* nothing is written before an approval: the workspace is listed at the
  moment the question appears and holds no plan artifact there.

The turns are real: the real engine, the real OpenAI-compatible tool loop and
the real permission gate, with only the model's stream scripted. The answer
goes through ``file_confirm.answer`` -- the same function
``delfin-agent approvals`` calls -- so the one-time-code and signature checks
stay on the path under test.

Controls that hold with and without this behaviour are included: an ordinary
write in default mode, and a planning turn that submits no plan. A case that
is green on the base commit measures nothing, so each red assertion here is
paired with one of those.
"""

from __future__ import annotations

import io
import json
import os
import threading
import time
from pathlib import Path

import pytest

from delfin.agent import api_client as _A
from delfin.agent import cli as agent_cli
from delfin.agent import file_confirm as FC
from delfin.agent import mcp_client as M
from delfin.agent import turn_outcome as TO

PLAN = "1. Write notes.txt with the outcome\n2. Report what was written"


# --- the scripted model ----------------------------------------------------

class _Fn:
    def __init__(self, name=None, arguments=None):
        self.name, self.arguments = name, arguments


class _TCDelta:
    def __init__(self, index, id=None, name=None, arguments=None):
        self.index, self.id = index, id
        self.function = _Fn(name, arguments)


class _Delta:
    def __init__(self, content=None, tool_calls=None):
        self.content, self.tool_calls = content, tool_calls


class _Choice:
    def __init__(self, delta, finish=None):
        self.delta, self.finish_reason = delta, finish


class _Usage:
    prompt_tokens, completion_tokens = 5, 3


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


class _EmptyReg:
    def discover_all(self):
        return []

    def discover_resources(self):
        return []

    def discover_prompts(self):
        return []


def _call(name, args, call_id="c1"):
    return _Stream([_Chunk([_Choice(
        _Delta(tool_calls=[_TCDelta(0, id=call_id, name=name,
                                    arguments=json.dumps(args))]),
        finish="tool_calls")], usage=_Usage())])


def _answer(text="done."):
    return _Stream([_Chunk([_Choice(_Delta(content=text), finish="stop")],
                           usage=_Usage())])


SUBMIT = _call("exit_plan_mode", {"plan": PLAN})
WRITE = _call("write_file", {"path": "notes.txt", "content": "planned\n"},
              "c2")


# --- the headless run, end to end ------------------------------------------

@pytest.fixture
def headless(monkeypatch, tmp_path):
    """``run(argv, streams, decide=...) -> (rc, engine, seen)``.

    ``decide(request_record) -> bool`` answers each question the run
    publishes (True approves); with ``decide=None`` nobody answers, which is
    what ``--confirm-timeout`` then expires. ``seen`` records what the
    approvals directory held: the first question seen, with the workspace
    listed at that moment, then each answer given.
    """
    home = tmp_path / "home"
    ws = home / "ws"
    ws.mkdir(parents=True)
    (ws / "readme.txt").write_text("readable\n", encoding="utf-8")
    monkeypatch.chdir(ws)                  # cmd_run moves into --cwd
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))
    # A question and its one-time code stay inside this test.
    room = tmp_path / "approvals"
    monkeypatch.setattr(FC, "requests_dir", lambda: room)
    from delfin.agent import approval_code as AC
    monkeypatch.setattr(AC, "codes_dir", lambda: tmp_path / "codes")

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
    monkeypatch.setattr(agent_cli, "_resume_or_create", lambda engine, a: "")
    # The off-list command behind some refusals must not run unasked here,
    # or the verdict would depend on this host's cage.
    monkeypatch.setattr(_A, "_caged_command_may_run_unasked",
                        lambda *a, **k: False)

    stdin = io.StringIO("")
    stdin.isatty = lambda: False           # type: ignore[method-assign]
    monkeypatch.setattr("sys.stdin", stdin)

    real_build = agent_cli._build_engine
    built: dict = {}

    def run(argv, streams, *, decide=None):
        calls = {"n": 0}

        def _create(**kwargs):
            s = streams[min(calls["n"], len(streams) - 1)]
            calls["n"] += 1
            return s

        def _build(args):
            engine = real_build(args)
            engine.client.client.chat.completions.create = _create
            built["engine"] = engine
            return engine

        monkeypatch.setattr(agent_cli, "_build_engine", _build)
        seen: list = []
        stop = threading.Event()
        if decide is not None:
            watcher = threading.Thread(target=_watch, daemon=True,
                                       args=(room, ws, decide, seen, stop))
            watcher.start()
        try:
            rc = agent_cli.main(["-p", "do the task", "--cwd", str(ws),
                                 "--provider", "ollama", "--model", "m",
                                 "--output-format", "json", *argv])
        finally:
            if decide is not None:
                stop.set()
                watcher.join(timeout=10)
        return rc, built.get("engine"), seen

    run.workspace = ws
    run.room = room
    return run


def _watch(room, ws, decide, seen, stop):
    """Answer every question this run publishes, the way a watcher would."""
    answered: set = set()
    deadline = time.time() + 60
    while not stop.is_set() and time.time() < deadline:
        for path in sorted(room.glob("*.request.json")):
            try:
                rec = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                continue
            rid = str(rec.get("id") or "")
            if not rid or rid in answered:
                continue
            if not seen:
                # The measured state of the workspace BEFORE any approval.
                seen.append({"tool": rec.get("tool", ""),
                             "ws": sorted(os.listdir(ws))})
            answered.add(rid)
            yes = bool(decide(rec))
            seen.append({"answered": yes,
                         "written": FC.answer(rid, yes, room=room,
                                              by="test")})
        time.sleep(0.05)


def _result(capsys) -> dict:
    """The run's result object. ``--output-format json`` prints exactly one,
    as a bare object (the ``type`` key is the stream-json events' shape)."""
    raw = capsys.readouterr().out
    events = []
    for line in raw.splitlines():
        if line.strip().startswith("{"):
            try:
                events.append(json.loads(line))
            except ValueError:
                pass
    for e in reversed(events):
        if "stop_reason" in e:
            return e
    raise AssertionError(f"no result object in {raw[-1500:]!r}")


# --- the plan reaches the person, and an approval lets it run ---------------

def test_a_supervised_run_puts_its_plan_where_a_person_can_answer_it(headless,
                                                                      capsys):
    """Red on the base: --supervised bound no plan channel, so the run was
    told nobody could be asked and no request file ever appeared."""
    rc, engine, seen = headless(
        ["--permission-mode", "plan", "--unattended", "--supervised",
         "--confirm-timeout", "25"],
        [SUBMIT, WRITE, _answer("Wrote it.")], decide=lambda rec: True)

    asked = [s for s in seen if s.get("tool") == "exit_plan_mode"]
    assert asked, f"no plan question reached the approvals dir: {seen}"
    assert engine.kit_permissions.plan_approval_callback is not None, (
        "--supervised bound no plan-approval channel")
    # Nothing existed in the workspace before the approval was given.
    assert "notes.txt" not in asked[0]["ws"], asked[0]["ws"]
    # And the approved plan did execute: the file exists only after it.
    assert (headless.workspace / "notes.txt").is_file()
    assert rc == 0 or rc == TO.EXIT_CODES["incomplete"], rc


def test_a_refused_plan_ends_the_turn_as_a_refusal_not_as_success(headless,
                                                                  capsys):
    """Red on the base: a person's "no" left the turn reporting done, exit
    0 -- the same answer a finished turn gives."""
    rc, engine, seen = headless(
        ["--permission-mode", "plan", "--unattended", "--supervised",
         "--confirm-timeout", "25"],
        [SUBMIT, _answer("Understood, stopping.")],
        decide=lambda rec: False)

    assert [s for s in seen if s.get("tool") == "exit_plan_mode"], seen
    assert engine.kit_permissions.mode == "plan", "a refusal must not leave plan"
    assert not (headless.workspace / "notes.txt").exists()
    assert rc == TO.EXIT_CODES["refused"] == 5, rc
    assert _result(capsys)["stop_reason"] == "refused"


def test_an_unanswered_plan_is_not_reported_as_a_finished_turn(headless,
                                                                capsys):
    """Red on the base: the turn ended done with exit 0, so a supervisor
    resuming it saw success and re-asked a question still on screen."""
    rc, engine, seen = headless(
        ["--permission-mode", "plan", "--unattended", "--supervised",
         "--confirm-timeout", "3"],
        [SUBMIT, _answer("Waiting."), _answer("Still waiting.")])

    assert engine.kit_permissions.mode == "plan"
    assert not (headless.workspace / "notes.txt").exists()
    assert rc == TO.EXIT_CODES["awaiting_approval"] == 11, rc
    assert _result(capsys)["stop_reason"] == "awaiting_approval"


def test_an_unsupervised_run_reports_its_plan_as_unanswered(headless, capsys):
    """No watcher at all: the plan genuinely cannot be approved here, and the
    outcome has to say so rather than claim the work is finished."""
    rc, engine, _ = headless(["--permission-mode", "plan", "--unattended"],
                             [SUBMIT, _answer()])
    assert engine.kit_permissions.mode == "plan"
    assert rc == TO.EXIT_CODES["awaiting_approval"] == 11, rc


# --- controls: these hold on the base commit too ----------------------------

def test_a_plain_write_in_default_mode_is_still_done(headless, capsys):
    """Control. Nothing about a plan, nothing about approval: a normal
    permitted turn still ends done with exit 0."""
    rc, _, _ = headless(["--permission-mode", "default", "--unattended"],
                        [WRITE, _answer("Wrote it.")])
    assert rc == 0, rc
    assert _result(capsys)["stop_reason"] == "done"


def test_a_planning_turn_that_submits_nothing_is_not_awaiting(headless,
                                                               capsys):
    """Control. Reading and answering in plan mode is the ordinary case and
    must not pick up a waiting-for-you outcome."""
    read = _call("read_file", {"path": "readme.txt"})
    rc, _, _ = headless(["--permission-mode", "plan", "--unattended"],
                        [read, _answer("Read.")])
    assert rc == 0, rc
    assert _result(capsys)["stop_reason"] == "done"


def test_an_empty_plan_is_not_a_submission(headless, capsys):
    """Control with a failing edge: a submission with no plan body was
    refused by the tool, so the turn is not waiting on anyone."""
    empty = _call("exit_plan_mode", {"plan": "   "})
    rc, _, _ = headless(["--permission-mode", "plan", "--unattended"],
                        [empty, _answer()])
    assert rc == 0, rc
    assert _result(capsys)["stop_reason"] == "done"
