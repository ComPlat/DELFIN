"""The terminal session runs auto iterations at its idle prompt.

Input: TerminalAgent.run() over a scripted engine whose "turn" writes a
file into a real git workspace (no model), a scripted keyboard
(``_line_ready`` / ``read_line``) and an AutoController with a fake clock.
Output: the turns the engine saw, the commits on the auto branch, what
the veto said, and the transcript.

Judged: ``/auto on`` is a command, not a turn; an idle prompt with work in
scope starts an iteration as a harness turn (no push grant: HarnessText);
each iteration is committed; the budget ends it with a report on screen; a
line typed during an iteration refuses its next tool call and asks the
engine to stop, the user's line runs next, and auto mode stays paused.
"""

from __future__ import annotations

import io
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest

from delfin.agent import agent_tasks, auto_mode
from delfin.agent import repl as R
from delfin.agent.harness_rows import HarnessText


def _git(ws, *a):
    return subprocess.run(["git", "-C", str(ws), *a], capture_output=True,
                          text=True, check=True).stdout.strip()


@pytest.fixture()
def ws(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))
    w = tmp_path / "ws"
    w.mkdir()
    _git(w, "init", "-q", "-b", "main")
    _git(w, "config", "user.name", "t")
    _git(w, "config", "user.email", "t@example.org")
    (w / "notes.txt").write_text("0\n", encoding="utf-8")
    _git(w, "add", "-A")
    _git(w, "commit", "-q", "-m", "init")
    return w


class Engine:
    session_id = "s-auto"
    last_turn_stop_reason = "end_turn"

    def __init__(self, ws, on_turn=None):
        self.ws = ws
        self.kit_permissions = SimpleNamespace(mode="acceptEdits",
                                               auto_veto=None)
        self.token_usage = {"input": 0, "output": 0}
        self.messages = []
        self.prompts = []
        self.stops = 0
        self.on_turn = on_turn
        self.pending_turn_continuation = ""

    def stream_response(self, user_message="", max_tokens=0, **kw):
        self.prompts.append(user_message)
        self.token_usage["input"] += 100
        answer = "done"
        if self.on_turn:
            answer = self.on_turn(self, user_message) or answer
        self.messages.append({"role": "user", "content": str(user_message)})
        self.messages.append({"role": "assistant", "content": answer})
        return answer

    def request_stop(self):
        self.stops += 1

    def clear_stop(self):
        pass

    def get_status(self):
        return {}

    def export_state(self):
        return {"engine_messages": list(self.messages),
                "token_usage": dict(self.token_usage)}


class Clock:
    t = 5_000_000.0

    def __call__(self):
        Clock.t += 1
        return Clock.t


def _agent(ws, engine, keys, *, max_iterations=2):
    """keys: a list of typed lines; None means 'nothing typed, idle'."""
    out, err = io.StringIO(), io.StringIO()
    out.isatty = lambda: False
    err.isatty = lambda: False
    script = list(keys)

    def read_line(_p=""):
        if not script:
            raise EOFError
        k = script.pop(0)
        if k is None:
            raise EOFError
        return k

    agent = R.TerminalAgent(engine, opts=R.ReplOptions(cwd=ws, max_tokens=0),
                            out=out, err=err, read_line=read_line)
    agent._stdin = type("S", (), {"isatty": lambda self: False})()
    agent._WAKE_EVERY_S = 0.0
    idle_polls = {"n": 0}

    def ready():
        if script and script[0] is None:
            # idle: let the wake run a few times, then the next key
            idle_polls["n"] += 1
            if idle_polls["n"] > 3:
                idle_polls["n"] = 0
                script.pop(0)
                return True
            return False
        return True

    agent._line_ready = ready
    agent._auto = auto_mode.AutoController(
        ws, session_id=engine.session_id, clock=Clock(),
        run_check=lambda e, c: (1, ""), notify=lambda *a: None,
        settings={"idle_seconds": 0, "max_iterations": max_iterations,
                  "plateau": 9, "repeat_failures": 9})
    return agent, err


def _task(ws, engine, subject="write the notes"):
    store = agent_tasks.get_store(ws)
    t = store.create(subject, session_id=engine.session_id)
    store.update(t["id"], status="in_progress")


def test_auto_on_is_a_command_and_iterations_follow_at_idle(ws):
    n = {"i": 0}

    def turn(engine, prompt):
        n["i"] += 1
        (ws / "notes.txt").write_text(f"{n['i']}\n", encoding="utf-8")
        return f"wrote {n['i']}\nLOOP: progress notes {n['i']}"

    eng = Engine(ws, on_turn=turn)
    _task(ws, eng)
    agent, err = _agent(ws, eng, ["/auto on", None, None, None])
    agent.run()
    # Two iterations, both harness turns, no turn for the command itself.
    assert len(eng.prompts) == 2, eng.prompts
    assert all(isinstance(p, HarnessText) for p in eng.prompts)
    assert all("write the notes" in p for p in eng.prompts)
    branch = _git(ws, "branch", "--show-current")
    assert branch.startswith("auto/")
    subjects = _git(ws, "log", "--format=%s", "-3").splitlines()
    assert subjects[0].startswith("auto iteration 2")
    assert subjects[1].startswith("auto iteration 1")
    text = err.getvalue()
    assert "auto mode on" in text
    assert "auto stopped" in text and "2 of 2 iterations" in text
    assert eng.kit_permissions.auto_veto is None


def test_a_line_typed_during_an_iteration_comes_first(ws):
    seen = {}

    def turn(engine, prompt):
        veto = engine.kit_permissions.auto_veto
        if not isinstance(prompt, HarnessText):
            seen["user_turn_veto"] = veto
            seen["paused_in_user_turn"] = agent._auto.paused
            return "answered the user"
        seen["before"] = veto("read_file", {"path": "notes.txt"})
        agent.queued.append("what are you doing?")   # typed meanwhile
        seen["after"] = veto("read_file", {"path": "notes.txt"})
        return "stopping\nLOOP: progress none"

    eng = Engine(ws, on_turn=turn)
    _task(ws, eng)
    agent, err = _agent(ws, eng, ["/auto on", None, None])
    agent.run()
    assert seen["before"] is None
    assert seen["after"] and "user" in seen["after"]
    assert eng.stops >= 1
    # The user's line ran next, outside the frame, and auto stayed paused.
    assert "what are you doing?" in eng.prompts[1]
    assert not isinstance(eng.prompts[1], HarnessText)
    assert seen["user_turn_veto"] is None
    assert len(eng.prompts) == 2
    assert seen["paused_in_user_turn"]
    # Leaving the session ends auto mode with its report.
    assert "session ended" in agent._auto.stop_report


def test_a_push_inside_an_iteration_is_queued(ws):
    seen = {}

    def turn(engine, prompt):
        seen["push"] = engine.kit_permissions.auto_veto(
            "bash", {"command": "git push origin HEAD"})
        return "LOOP: needs-human may I push?"

    eng = Engine(ws, on_turn=turn)
    _task(ws, eng)
    agent, err = _agent(ws, eng, ["/auto on", None])
    agent.run()
    assert "queued" in seen["push"]
    kinds = [n["kind"] for n in agent._auto.needs_you]
    assert "outward" in kinds and "question" in kinds


def test_the_window_and_terminal_share_the_controller_in_commands(ws):
    eng = Engine(ws)
    agent, _ = _agent(ws, eng, [])
    assert agent._ctx().auto is agent._auto
