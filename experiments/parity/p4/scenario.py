"""Measure package p4: a subagent in its own worktree, and undoing one turn.

Run it; do not read it as a description of the result::

    python3 experiments/parity/p4/scenario.py > experiments/parity/p4/before.json

It prints ONE JSON object: the command, the tree each measure loaded, and
per variant whether the capability was REACHED, USED and CORRECT.

The two claims under test
------------------------
UNDO. ``undo_changes(scope="turn")`` reverts exactly the files THIS turn
touched and nothing else. The turn window is the per-session
``ToolContext``: a client empties it at the start of every turn
(``ToolContext.begin_turn``) and every captured file change appends its
journal seq to it. The question asked here is the one a person asks: after a
turn that wrote a file, does a LATER turn told to undo "this turn" take that
file back, or say it has nothing of its own to undo?

ISOLATION. A subagent asked to change files can run in its own git worktree,
so two writers never share a tree. The question asked here: does the
product's own fan-out point a writer at a second checkout, and is the
parent's tree still what it was afterwards?

How it drives the product
-------------------------
Never a helper called directly for the thing under test. A session is built
the way the CLI builds one -- ``cli.build_parser()`` then
``cli._build_engine(args)``, the real ``AgentEngine`` with the real tool
loop, permission gate, change journal and undo tool -- and each turn goes
through ``repl.run_turn``, the terminal's own per-turn entry point. Two turns
run on the SAME engine, which is what a turn-scoped undo needs to be
measured at all.

Only the model's HTTP response is scripted:
``engine.client.client.chat.completions.create`` is replaced by a stub that
returns the streamed chunks a reply would have carried, the technique of
``tests/test_a_headless_run_reports_a_state_a_program_can_read.py``. Nothing
here opens a network connection, starts a live model, or touches anyone's
calc/ or archive/ folder: HOME and the DELFIN state root move into a scratch
directory BEFORE DELFIN is imported, so even a store that fixes its path at
import time lands there.

One variant per process: the dashboard and the widget machinery bind
process-global state, and a second session in the same process interferes
with the first -- a harness artefact this measurement must not contain.

What each verdict means
-----------------------
REACHED  the capability's own code ran: an undo reached
         ``change_journal.revert`` (proven by the ``tool="undo"`` record it
         writes), isolation reached ``worktree.enter_worktree`` (proven by a
         second checkout appearing while the delegate ran).
USED     the product path invoked it and answered with its own result shape
         (a real ``undo_changes`` verdict object, a real worktree handed to
         the delegate) -- not "Unknown tool", not an exception. Whether
         acting was RIGHT is CORRECT's question, so USED cannot drift with
         the defect it measures.
CORRECT  the outcome is right: exactly the files of THAT turn came back and
         no others (undo); the delegate wrote in the other checkout while the
         parent's tree stayed untouched (isolation).
"""

from __future__ import annotations

import contextlib
import io
import json
import os
import subprocess
import sys
import tempfile
import types
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]

# Scratch first, imports second: a store that resolves its path at import
# time would otherwise bind the real ~/.delfin (the bash_jobs lesson).
SCRATCH = Path(os.environ.get("DELFIN_SCRATCH")
               or Path(tempfile.mkdtemp(prefix="p4-parity-"))).resolve()
SCRATCH.mkdir(parents=True, exist_ok=True)
HOME = SCRATCH / "home"
(HOME / ".delfin").mkdir(parents=True, exist_ok=True)
os.environ["HOME"] = str(HOME)
os.environ["DELFIN_STATE"] = str(HOME / ".delfin")
os.environ["DELFIN_CLASSIFIER"] = "off"
os.environ["DELFIN_SETTINGS_DEFAULTS"] = ""

if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

MARK = "p4-scenario"
PERSONS_FILE = "keep.txt"
PERSONS_TEXT = "written by the person, not the agent\n"


def _write(path, content):
    return {"tool": "write_file", "args": {"path": path, "content": content}}


def _undo(scope):
    return {"tool": "undo_changes", "args": {"scope": scope}}


# --------------------------------------------------------------------------
# the scripted model reply
# --------------------------------------------------------------------------

class _Part:
    def __init__(self, **kw):
        self.__dict__.update(kw)


class _Chunk:
    """One streamed chunk, in the shape the OpenAI SDK hands the client."""

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
    """The scripted answer to the NEXT request: one step per round.

    The round is the number of assistant messages AFTER the last user
    message -- the rounds of the CURRENT turn. Counting assistant messages
    over the whole conversation instead would make a second turn of one
    session start already past the end of its own plan and answer "Nothing
    more" without ever issuing the call under test.
    """

    def _create(**kwargs):
        messages = kwargs.get("messages") or []
        start = 0
        for i, m in enumerate(messages):
            if m.get("role") == "user":
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
        return _Stream([_Chunk(
            _Part(content=step.get("say"), tool_calls=[call]), "tool_calls")])

    return _create


# --------------------------------------------------------------------------
# one session: the CLI's engine, the terminal's turn
# --------------------------------------------------------------------------

class Session:
    """A real agent session in its own workspace; ``turn`` runs one turn."""

    def __init__(self, root: Path) -> None:
        self.root = root
        self.workspace = root / "ws"
        for name in ("calc", "archive", "office"):
            (root / name).mkdir(parents=True, exist_ok=True)
        self.workspace.mkdir(parents=True, exist_ok=True)
        (self.workspace / PERSONS_FILE).write_text(PERSONS_TEXT)
        self._engine = None

    def open(self) -> "Session":
        from delfin.agent import cli as agent_cli

        parser = agent_cli.build_parser()
        args = parser.parse_args([
            "chat", "-p", "unused", "--cwd", str(self.workspace),
            "--provider", "ollama", "--model", "scripted-model",
            "--permission-mode", "acceptEdits"])
        self._engine = agent_cli._build_engine(args)
        return self

    @property
    def engine(self):
        return self._engine

    @property
    def sid(self) -> str:
        return str(getattr(self._engine, "session_id", "") or "")

    def turn(self, steps: list[dict]) -> dict:
        """One turn: the terminal's run_turn over the scripted reply."""
        from delfin.agent import repl as _repl
        from delfin.agent.api_client import OpenAIClient

        if not isinstance(self._engine.client, OpenAIClient):
            raise AssertionError(
                f"the session runs on a {type(self._engine.client).__name__}, "
                "not the scripted OpenAI-compatible client")
        self._engine.client.client.chat.completions.create = _reply(steps)
        # run_turn answers with a TurnResult; the driver only needs the turn
        # to have run and the transcript to hold its tool results.
        return _repl.run_turn(self._engine, f"{MARK} turn",
                              sink=lambda item: None, memory_context="")

    def window(self) -> list[int]:
        """The undo window the product's own turn-scoped undo reads."""
        from delfin.agent.api_client import _doc_executor
        return list(_doc_executor._turn_seqs_for(self.sid))


# --------------------------------------------------------------------------
# UNDO
# --------------------------------------------------------------------------

def undo_variant(root: Path, spec: dict) -> dict:
    """Write in turn 1, undo in turn 2, report who/what came back."""
    from delfin.agent import change_journal as CJ
    from delfin.agent import repl_commands as RC

    files = spec["files"]
    names = sorted(files)
    s = Session(root).open()
    s.turn([_write(n, files[n]) for n in names] + [{"text": "wrote them"}])
    seqs_t1 = s.window()
    wrote = [n for n in names if (s.workspace / n).read_text() == files[n]]
    if len(wrote) != len(names):
        return {"reached": False, "used": False, "correct": False,
                "note": f"turn 1 wrote {wrote}, wanted {names}"}

    steps: list[dict] = []
    if spec.get("second_turn_writes"):
        steps.append(_write(spec["second_turn_writes"], "turn two\n"))
    if spec.get("surface", "tool") == "tool":
        steps.append(_undo(spec.get("scope", "turn")))
    s.turn(steps + [{"text": "done"}])
    seqs_t2 = s.window()

    # The terminal builtin, driven over the same session state -- the second
    # surface that reads the same window.
    repl_says = ""
    if spec.get("surface") == "repl":
        ctx = types.SimpleNamespace(
            engine=s.engine, workspace=s.workspace, session_id=s.sid,
            ask=None, say=None, auto=None)
        repl_says = str(getattr(
            RC.dispatch(f"/undo-file {spec.get('scope', 'turn')}", ctx),
            "output", "") or "")

    journal = CJ.list_changes(s.sid, last_n=40)
    undo_recs = [r for r in journal if str(r.get("tool", "")) == "undo"]
    reverted = sorted({Path(str(r.get("path", ""))).name for r in undo_recs})
    answer = _last_tool_answer(s.engine, "undo_changes")

    expect = sorted(spec.get("expect_reverted", []))
    holds = {}
    for n in names:
        p = s.workspace / n
        holds[n] = bool(p.exists() and p.read_text() == files[n])
    reverted_ok = reverted == expect
    # What the undo must NOT have touched still holds turn 1's content, and
    # what it must have taken back no longer does.
    holds_ok = all(holds.get(n) is True for n in spec.get("must_hold", []))
    moved_ok = all((not (s.workspace / n).exists())
                   or (s.workspace / n).read_text() != files[n]
                   for n in expect)
    persons_ok = (PERSONS_FILE in expect) or (
        (s.workspace / PERSONS_FILE).read_text() == PERSONS_TEXT)
    return {
        # REACHED: the capability's code ran.
        "reached": bool(answer) or bool(undo_recs) or bool(repl_says),
        # USED: the product's undo surface answered with its own shape.
        "used": _is_undo_verdict(answer) or _is_repl_verdict(repl_says),
        "correct": bool(reverted_ok and holds_ok and moved_ok and persons_ok),
        # What the undo took back although its scope never touched it: the
        # defect itself, counted apart so USED stays a wiring measure.
        "wrongful_revert": sorted(set(reverted) - set(expect)),
        "surface": spec.get("surface", "tool"),
        "scope": spec.get("scope", "turn"),
        "seqs_turn1": seqs_t1,
        "seqs_turn2": seqs_t2,
        "reverted": reverted,
        "expected_reverted": expect,
        "turn1_content_still_there": holds,
        "persons_file_untouched": persons_ok,
        "undo_answer": answer[:240],
        "repl_says": repl_says.strip()[:160],
    }


def _last_tool_answer(engine, tool: str) -> str:
    """The output the last call of ``tool`` got, read off the transcript."""
    calls: dict = {}
    answer = ""
    for msg in (getattr(engine, "messages", None) or []):
        if not isinstance(msg, dict):
            continue
        for call in (msg.get("tool_calls") or []):
            name = ((call.get("function") or {}).get("name") or "")
            calls[str(call.get("id", ""))] = name
        if msg.get("role") == "tool" and calls.get(
                str(msg.get("tool_call_id", ""))) == tool:
            answer = str(msg.get("content") or "")
    return answer


def _is_undo_verdict(answer: str) -> bool:
    """True when the answer is undo_changes' own result object.

    A refusal, an "Unknown tool" string or an exception is not a verdict:
    naming a capability in a catalogue is not a capability.
    """
    try:
        got = json.loads(answer)
    except (json.JSONDecodeError, TypeError):
        return False
    return isinstance(got, dict) and bool(
        {"reverted", "note", "error"} & set(got))


def _is_repl_verdict(text: str) -> bool:
    """True when /undo-file answered with its own report header."""
    return text.strip().startswith("undo (")


# --------------------------------------------------------------------------
# ISOLATION
# --------------------------------------------------------------------------

class _ScriptedParent:
    """The parent-client surface ``run_subagent`` uses for the delegate.

    The delegate's own stream is scripted here rather than over HTTP: what
    this measures is which WORKSPACE the delegate was pointed at, and the
    write it performs goes through that workspace path, so the isolation
    decision -- not the model -- is what is under test.
    """

    model = "scripted-model"
    _turn_thinking_budget = 0

    def __init__(self, steps: list[dict], perms) -> None:
        self._steps = steps
        self._perms = perms
        self.wrote_in: list[str] = []

    def set_permissions(self, perms):
        self._perms = perms

    def switch_model(self, name):
        self.model = name

    def should_stop(self):
        return False

    def stream_message(self, *, messages, system, max_tokens, **kwargs):
        from delfin.agent.api_client import StreamEvent

        workspace = Path(getattr(self._perms, "workspace", "") or ".")
        for step in self._steps:
            if step.get("tool") != "write_file":
                continue
            args = step["args"]
            path = workspace / args["path"]
            try:
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(args.get("content", ""))
                self.wrote_in.append(str(path))
                out = f"wrote {path}"
            except OSError as exc:
                out = f"error: {exc}"
            yield StreamEvent(type="tool_use", tool_name="write_file",
                              tool_input=json.dumps(args))
            yield StreamEvent(type="tool_result", tool_name="write_file",
                              tool_output=out)
        yield StreamEvent(type="text_delta", text="done")


def _git(repo: Path, *args) -> str:
    return subprocess.run(["git", "-C", str(repo), *args],
                          capture_output=True, text=True).stdout


@contextlib.contextmanager
def _silence():
    old = sys.stdout, sys.stderr
    sys.stdout = sys.stderr = io.StringIO()
    try:
        yield
    finally:
        sys.stdout, sys.stderr = old


def _writer_repo(root: Path):
    """A committed git workspace plus permissions that may write in it."""
    from delfin.agent import api_client as A

    ws = root / "repo"
    ws.mkdir(parents=True, exist_ok=True)
    _git(ws, "init", "-q")
    _git(ws, "config", "user.email", "test@example.invalid")
    _git(ws, "config", "user.name", "Test Person")
    (ws / "shared.py").write_text("A = 1\n")
    _git(ws, "add", "-A")
    _git(ws, "commit", "-qm", "start")
    perms = A.KitToolPermissions(workspace=ws, mode="acceptEdits")
    perms.task_session_id = "p4-isolation"
    perms.confirm_callback = lambda name, args, preview: True
    return ws, perms


def fanout_variant(root: Path, spec: dict) -> dict:
    """The engine's own parallel fan-out: what isolation does it choose?"""
    from delfin.agent import api_client as A

    _ws, perms = _writer_repo(root)
    n = int(spec.get("n", 2))
    calls = [{"id": f"c{i}", "type": "function", "function": {
        "name": "subagent",
        "arguments": json.dumps({
            "subagent_type": spec.get("preset", "general-purpose"),
            "description": f"arm {i}", "prompt": "edit shared.py"})}}
        for i in range(n)]
    seen: list[dict] = []
    # The fan-out submits to the MODULE-LEVEL executor singleton
    # (``_doc_executor.execute``), so that is the object to observe; spying on
    # a fresh _DocToolExecutor() would see nothing and report the capability
    # as unreachable when it is only watched in the wrong place.
    ex = A._doc_executor
    real = ex.execute

    def _spy(name, arguments, permissions=None, tool_context=None):
        if name == "subagent":
            seen.append(dict(arguments or {}))
            return json.dumps({"note": "measured, not run"})
        return real(name, arguments, permissions=permissions,
                    tool_context=tool_context)

    ex.execute = _spy
    futures = executor = None
    try:
        with _silence():
            futures, executor = A._fan_out_subagents(calls, perms)
        # The calls reach execute() on the pool's WORKER threads, so the spy
        # must stay in place until they have all run. Restoring it as soon as
        # the fan-out returned raced the workers: a late call went to the real
        # executor, the count came back short, and the variant measured the
        # harness's thread timing instead of the product's decision.
        if futures:
            import concurrent.futures as _cf
            _cf.wait(list(futures.values()), timeout=120)
    except Exception as exc:                        # noqa: BLE001
        return {"reached": False, "used": False, "correct": False,
                "error": f"{type(exc).__name__}: {exc}"[:200]}
    finally:
        if executor is not None:
            executor.shutdown(wait=True)
        ex.execute = real
    want = int(spec.get("expect_isolated", n))
    isolated = [c.get("isolation") == "worktree" for c in seen]
    # reached/used both require the fan-out to have handed every call to the
    # subagent path: a count that compares equal against 0 is not a verdict
    # (an assertion that cannot fail is not a test).
    got = sum(isolated)
    return {
        "reached": bool(futures) and len(seen) == n,
        "used": len(seen) == n and got == want,
        "correct": len(seen) == n and got == want,
        "submitted_isolation": isolated,
        "writers_isolated": got,
        "calls_seen": len(seen),
        "wanted_isolated": want,
        "note": "the fan-out's own decision; no delegate ran",
    }


def isolation_variant(root: Path, spec: dict) -> dict:
    """Writers really run, each asked for its own worktree."""
    from delfin.agent import subagents as SA
    from delfin.agent import worktree as WT

    ws, perms = _writer_repo(root)
    n = int(spec.get("n", 2))
    created: list[str] = []
    real_enter = WT.enter_worktree

    def _spy_enter(*args, **kwargs):
        info = real_enter(*args, **kwargs)
        created.append(str(info.path))
        return info

    WT.enter_worktree = _spy_enter
    results = []
    try:
        for i in range(n):
            parent = _ScriptedParent([_write(f"arm{i}.txt", f"arm {i}\n")],
                                     perms)
            with _silence():
                res = SA.run_subagent(
                    subagent_type=spec.get("preset", "general-purpose"),
                    description=f"arm {i}", prompt="edit shared.py",
                    parent_client=parent, parent_perms=perms,
                    isolation=spec.get("isolation", "worktree"),
                    inherit_context=False)
            results.append((parent, res))
    finally:
        WT.enter_worktree = real_enter

    errors = [r.error for _p, r in results if r.error]
    parent_dirty = _git(ws, "status", "--porcelain").strip()
    wrote = [w for p, _r in results for w in p.wrote_in]
    inside = [w for w in wrote
              if any(Path(w).is_relative_to(Path(t)) for t in created)]
    # A writer that ran WITHOUT isolation would have written into ws.
    in_parent = [w for w in wrote if Path(w).is_relative_to(ws.resolve())]
    distinct = len(set(created)) == n
    return {
        "reached": bool(created),
        "used": bool(created) and distinct and not in_parent,
        "correct": (bool(created) and distinct and not in_parent
                    and len(inside) == n and not parent_dirty and not errors),
        "worktrees": created,
        "trees_distinct": distinct,
        "delegate_wrote": [w.replace(str(ws), "<repo>") for w in wrote],
        "wrote_inside_a_worktree": len(inside),
        "wrote_in_parent": in_parent,
        "parent_dirty": parent_dirty[:200],
        "errors": [e[:200] for e in errors],
        "worktree_summary": [r.worktree for _p, r in results],
    }


# --------------------------------------------------------------------------
# the matrix
# --------------------------------------------------------------------------

def matrix() -> list[dict]:
    one = {"note.txt": "turn one\n"}
    return [
        # UNDO: a later turn undoes "this turn", which wrote nothing itself.
        {"id": "U1", "kind": "undo", "why": "one write, undo a turn later",
         "files": one, "expect_reverted": [], "must_hold": ["note.txt"]},
        {"id": "U2", "kind": "undo", "why": "two writes",
         "files": {"note.txt": "turn one\n", "second.md": "# two\n"},
         "expect_reverted": [], "must_hold": ["note.txt", "second.md"]},
        {"id": "U3", "kind": "undo", "why": "python source file",
         "files": {"prog.py": "x = 1\n"}, "expect_reverted": [],
         "must_hold": ["prog.py"]},
        {"id": "U4", "kind": "undo",
         "why": "CONTROL: the turn that wrote undoes its own turn -- only "
         "its own file may go",
         "files": one, "second_turn_writes": "t2.txt",
         "expect_reverted": ["t2.txt"], "must_hold": ["note.txt"]},
        {"id": "U5", "kind": "undo", "why": "long multi-line content",
         "files": {"long.txt": "\n".join(f"line {i}" for i in range(200))},
         "expect_reverted": [], "must_hold": ["long.txt"]},
        {"id": "U6", "kind": "undo", "why": "file in a subfolder",
         "files": {"sub/deep/note.txt": "deep\n"}, "expect_reverted": [],
         "must_hold": ["sub/deep/note.txt"]},
        {"id": "U7", "kind": "undo",
         "why": "the terminal builtin on the same state, not the model tool",
         "files": one, "surface": "repl", "expect_reverted": [],
         "must_hold": ["note.txt"]},
        {"id": "U8", "kind": "undo",
         "why": "CONTROL: session scope reaches back over turns, so it MUST "
         "take turn 1's file back",
         "files": one, "scope": "session", "expect_reverted": ["note.txt"]},
        # ISOLATION: writers must not share a tree.
        {"id": "I1", "kind": "fanout", "why": "fan-out of 2 writers", "n": 2},
        {"id": "I2", "kind": "fanout", "why": "fan-out of 3 writers", "n": 3},
        {"id": "I3", "kind": "fanout", "why": "fan-out of 6 writers (above "
         "the 4-worker pool)", "n": 6},
        {"id": "I4", "kind": "fanout", "why": "fan-out of read-only presets: "
         "no worktree wanted", "n": 2, "preset": "explore",
         "expect_isolated": 0},
        {"id": "I5", "kind": "isolation", "why": "two writers really run",
         "n": 2},
    ]


RUNNERS = {"undo": undo_variant, "fanout": fanout_variant,
           "isolation": isolation_variant}


def _run_one(variant_id: str) -> int:
    """This process runs ONE variant and prints its row on stdout.

    Everything a variant or a widget prints goes to stderr, so stdout holds
    exactly one JSON object for the parent to read.
    """
    spec = next((s for s in matrix() if s["id"] == variant_id), None)
    if spec is None:
        print(json.dumps({"id": variant_id, "error": "no such variant"}))
        return 2
    home = SCRATCH / variant_id
    home.mkdir(parents=True, exist_ok=True)
    row = {"id": spec["id"], "kind": spec["kind"], "why": spec["why"]}
    real_stdout = sys.stdout
    try:
        with contextlib.redirect_stdout(sys.stderr):
            row.update(RUNNERS[spec["kind"]](home, spec))
    except Exception as exc:                        # noqa: BLE001
        row.update({"reached": False, "used": False, "correct": False,
                    "error": f"{type(exc).__name__}: {exc}"[:300]})
    print(json.dumps(row), file=real_stdout, flush=True)
    return 0


def main() -> int:
    import delfin
    import delfin.agent.subagents as SA
    import delfin.agent.tool_context as TC
    import delfin.agent.worktree as WT

    args = list(sys.argv[1:])
    if args[:1] == ["--one"]:
        return _run_one(args[1].upper())

    # Parent: one clean measurement per variant, each in its own process.
    wanted = {part.strip().upper() for arg in args
              for part in arg.split(",") if part.strip()}
    rows = []
    for spec in matrix():
        if wanted and spec["id"] not in wanted:
            continue
        env = dict(os.environ, DELFIN_SCRATCH=str(SCRATCH),
                   PYTHONPATH=str(REPO))
        proc = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "--one",
             spec["id"]],
            capture_output=True, text=True, env=env, cwd=str(REPO),
            timeout=900)
        try:
            rows.append(json.loads(proc.stdout.strip().splitlines()[-1]))
        except (IndexError, json.JSONDecodeError):
            rows.append({"id": spec["id"], "kind": spec["kind"],
                         "why": spec["why"], "reached": False, "used": False,
                         "correct": False,
                         "error": f"no row from the child (rc "
                                  f"{proc.returncode}): "
                                  f"{proc.stderr.strip()[-200:]}"})

    def n(pred):
        return sum(1 for r in rows if r.get(pred) is True)

    text = json.dumps({
        "package": "p4",
        "command": "python3 experiments/parity/p4/scenario.py "
                   "> experiments/parity/p4/<before|after>.json",
        "scratch": str(SCRATCH),
        "n_variants": len(rows),
        "population": "synthetic scenarios on the terminal's real turn path "
                      "(cli._build_engine + repl.run_turn) in a temporary "
                      "workspace, with only the model reply scripted; no "
                      "live model, no network, no user calc/ or archive/ "
                      "folder",
        "loaded": {
            "delfin": str(getattr(delfin, "__file__", "")),
            "tool_context": str(TC.__file__),
            "subagents": str(SA.__file__),
            "worktree": str(WT.__file__),
            "head": _head(REPO),
        },
        "reached": n("reached"), "used": n("used"), "correct": n("correct"),
        "rows": rows,
    }, indent=2)
    # No home path or account name in a committed measurement: every path
    # here lives under the scratch root or this checkout, so both are
    # replaced by their placeholder before the text is written.
    for real, alias in ((str(SCRATCH), "<SCRATCH>"), (str(REPO), "<REPO>")):
        text = text.replace(real, alias)
    # Fail loudly rather than commit a personal path: every path this
    # measurement can name lives under one of the two roots above, so no
    # account directory may survive the aliases.
    if "/home/" in text:
        raise AssertionError("an unaliased personal path reached the output")
    print(text)
    return 0


def _head(repo: Path) -> str:
    try:
        return subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"],
                              capture_output=True, text=True).stdout.strip()
    except OSError:
        return ""


if __name__ == "__main__":
    raise SystemExit(main())
