"""A turn-scoped undo reverts the files of ITS turn, not of the turn before.

What the user asks ``undo_changes(scope="turn")`` is "take back what this
turn did". The window it reads is the per-session ``ToolContext``: a client
empties the window at the start of every turn (``ToolContext.begin_turn``)
and every captured file change adds its journal seq to it.

``begin_turn`` reset only the flat list. The session-keyed view of the same
window, ``_turn_seqs_by_session``, was cleared lazily -- and only when the
next change was recorded -- so a turn that writes nothing at all kept
reporting the PREVIOUS turn's seqs. The product then did the one thing the
scope promises never to do: it reverted a file the undoing turn never
touched, and reported it as a successful undo of this turn. The terminal's
``/undo-file turn`` reads the same window and did the same thing.

Three things, all driven through the terminal's own turn path
(``cli._build_engine`` then ``repl.run_turn``) with only the model's streamed
reply scripted:

* a turn that wrote nothing and asks to undo "this turn" -> reverts NOTHING
  (red on the base: it deleted the file the earlier turn wrote);
* CONTROL -- the turn that wrote a file undoes its own turn -> that file
  comes back, and an earlier turn's file is left alone;
* CONTROL -- ``scope="session"`` is meant to reach back over turns -> it
  still reverts the earlier file.

The first test fails on the base for exactly the reason it names, and the
two controls pass on both sides: that is what stops the first from being an
artefact of how the turn was scripted.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent import cli as agent_cli
from delfin.agent import repl as _repl


# --- the scripted streamed reply -------------------------------------------

class _Part:
    def __init__(self, **kw):
        self.__dict__.update(kw)


class _Chunk:
    def __init__(self, delta, finish):
        self.choices = [_Part(delta=delta, finish_reason=finish)]
        self.usage = _Part(prompt_tokens=20, completion_tokens=5,
                           total_tokens=25)


class _Stream:
    def __init__(self, chunks):
        self._chunks = chunks

    def __iter__(self):
        return iter(self._chunks)

    def close(self):
        pass


def _reply(steps: list[dict]):
    """Answer the next request with one step per round of THIS turn.

    The round is counted from the last user message, so a second turn of one
    session gets its own plan instead of starting past the end of the first.
    """

    def _create(**kwargs):
        messages = kwargs.get("messages") or []
        start = 0
        for i, msg in enumerate(messages):
            if msg.get("role") == "user":
                start = i + 1
        done = sum(1 for m in messages[start:] if m.get("role") == "assistant")
        if done >= len(steps):
            return _Stream([_Chunk(_Part(content="Nothing more.",
                                         tool_calls=None), "stop")])
        step = steps[done]
        if not step.get("tool"):
            return _Stream([_Chunk(_Part(content=step.get("text", "done"),
                                         tool_calls=None), "stop")])
        call = _Part(index=0, id=f"call_{done + 1}", function=_Part(
            name=step["tool"], arguments=json.dumps(step.get("args") or {})))
        return _Stream([_Chunk(_Part(content=step.get("say"),
                                     tool_calls=[call]), "tool_calls")])

    return _create


# --- one real session, two real turns --------------------------------------

@pytest.fixture
def session(tmp_path, monkeypatch):
    """A real agent session in a scratch workspace; two turns on one engine.

    Only the model's HTTP reply is scripted. The engine, the tool loop, the
    permission gate, the change journal and the undo tool are production's.
    """
    home = tmp_path / "home"
    (home / ".delfin").mkdir(parents=True)
    ws = tmp_path / "ws"
    ws.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("DELFIN_STATE", str(home / ".delfin"))
    monkeypatch.setenv("DELFIN_CLASSIFIER", "off")
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))
    monkeypatch.chdir(ws)

    args = agent_cli.build_parser().parse_args([
        "chat", "-p", "unused", "--cwd", str(ws), "--provider", "ollama",
        "--model", "scripted-model", "--permission-mode", "acceptEdits"])
    engine = agent_cli._build_engine(args)

    def _turn(steps):
        engine.client.client.chat.completions.create = _reply(steps)
        return _repl.run_turn(engine, "do it", sink=lambda item: None,
                              memory_context="")

    _turn.workspace = ws
    _turn.engine = engine
    return _turn


def _write(path, content):
    return {"tool": "write_file", "args": {"path": path, "content": content}}


def _undo(scope):
    return {"tool": "undo_changes", "args": {"scope": scope}}


def _undo_answer(engine) -> dict:
    """What the product answered the last undo_changes call."""
    names: dict = {}
    answer = ""
    for msg in (getattr(engine, "messages", None) or []):
        for call in (msg.get("tool_calls") or []):
            names[str(call.get("id", ""))] = (
                (call.get("function") or {}).get("name") or "")
        if msg.get("role") == "tool" and names.get(
                str(msg.get("tool_call_id", ""))) == "undo_changes":
            answer = str(msg.get("content") or "")
    assert answer, "the undo call never reached the tool executor"
    got = json.loads(answer)
    assert isinstance(got, dict), f"undo_changes answered {answer[:120]!r}"
    return got


# --- the defect -------------------------------------------------------------

def test_a_turn_that_wrote_nothing_undoes_nothing(session):
    """The scope is THIS turn; a turn with no writes has nothing to undo.

    Red on the base: the per-session seq entry of the previous turn survives
    begin_turn, so the undo reverts the file that earlier turn wrote -- a
    file this turn never touched, deleted and reported as this turn's undo.
    """
    ws, engine = session.workspace, session.engine
    (ws / "keep.txt").write_text("written by the person\n")

    session([_write("note.txt", "turn one\n"), {"text": "wrote it"}])
    assert (ws / "note.txt").read_text() == "turn one\n", (
        "turn 1 wrote nothing: the harness, not the undo, is at fault")

    session([_undo("turn"), {"text": "undid it"}])

    answer = _undo_answer(engine)
    assert answer.get("reverted") == [], (
        f"scope=turn reverted {answer.get('reverted')!r}, but this turn "
        "changed no file -- the window still held the previous turn's seqs")
    assert (ws / "note.txt").exists() and (
        ws / "note.txt").read_text() == "turn one\n", (
        "the undo took back a file this turn never touched")
    assert (ws / "keep.txt").read_text() == "written by the person\n"


def test_the_terminal_builtin_does_not_revert_the_previous_turn_either(
        session):
    """/undo-file turn reads the same window and must answer the same way."""
    from delfin.agent import repl_commands as RC

    ws = session.workspace
    session([_write("note.txt", "turn one\n"), {"text": "wrote it"}])
    session([{"text": "nothing to do here"}])

    ctx = _ReplContext(engine=session.engine, workspace=ws,
                       session_id=session.engine.session_id)
    out = str(getattr(RC.dispatch("/undo-file turn", ctx), "output", "") or "")
    assert "restored" not in out, f"/undo-file turn reverted a file: {out!r}"
    assert (ws / "note.txt").exists() and (
        ws / "note.txt").read_text() == "turn one\n", (
        f"/undo-file turn undid the previous turn: {out!r}")


class _ReplContext:
    """The handle /undo-file reads: an engine, its workspace and session id."""

    def __init__(self, engine, workspace: Path, session_id: str) -> None:
        self.engine = engine
        self.workspace = workspace
        self.session_id = session_id
        self.ask = None
        self.say = None
        self.auto = None


# --- controls: these must pass on the base too ------------------------------

def test_the_turn_that_wrote_undoes_its_own_write(session):
    """CONTROL: the scope does work when the turn really did write."""
    ws, engine = session.workspace, session.engine
    session([_write("note.txt", "turn one\n"), {"text": "wrote it"}])
    session([_write("t2.txt", "turn two\n"), _undo("turn"),
             {"text": "undid my own write"}])

    answer = _undo_answer(engine)
    reverted = [Path(str(p)).name for p in (answer.get("reverted") or [])]
    assert reverted == ["t2.txt"], answer
    assert not (ws / "t2.txt").exists(), "the undo left its own file behind"
    assert (ws / "note.txt").read_text() == "turn one\n", (
        "scope=turn reached into the previous turn")


def test_session_scope_still_reaches_back_over_turns(session):
    """CONTROL: session scope is the one meant to unwind earlier changes."""
    ws, engine = session.workspace, session.engine
    session([_write("note.txt", "turn one\n"), {"text": "wrote it"}])
    session([_undo("session"), {"text": "undid the session"}])

    answer = _undo_answer(engine)
    reverted = [Path(str(p)).name for p in (answer.get("reverted") or [])]
    assert reverted == ["note.txt"], answer
    assert not (ws / "note.txt").exists()
