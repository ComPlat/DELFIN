"""One agent session, driven over the TUI protocol on stdin/stdout.

    python -m delfin.agent.service.worker [--cwd DIR] [--resume ID] [--fake]

The worker runs the session loop the terminal REPL runs -- ``TerminalAgent``
itself, not a copy of it -- with the terminal swapped out:

* typed lines arrive as ``session.prompt`` requests instead of keystrokes,
  and are handled exactly as the REPL handles a line (messages, /commands,
  !shell, # notes, queueing during a turn, steering into it);
* what the REPL would print arrives at the window as events
  (``transcript.ProtocolTranscript``);
* approval questions go to the window (``approvals.ProtocolBroker``) with
  the options the REPL offers, and wait for its answer.

Everything that makes a session behave -- continuing a turn cut off at the
token limit, waiting out an endpoint outage, waking for a finished job or a
message from another session, presence, the per-turn checkpoint, the
session store, the end-of-session report and learning -- is the REPL's own
code, so the TUI session and ``delfin-agent chat`` cannot drift apart.

stdout carries protocol messages and nothing else: file descriptors 0 and 1
are taken over at start-up and pointed at /dev/null and stderr, so neither
a stray print nor a child process can write into the protocol stream.
"""

from __future__ import annotations

import argparse
import os
import queue
import sys
import threading
import time
from pathlib import Path
from typing import Any

from . import protocol
from .approvals import ProtocolBroker
from .eventlog import EventLog, log_path
from .transcript import ProtocolTranscript

_EOF = object()

#: How often the usage counters are sampled while a turn runs.
_USAGE_TICK_S = 0.5


class SessionRefused(Exception):
    """The session cannot be opened; the message says why and what to do."""


# ---------------------------------------------------------------------------
# stdio
# ---------------------------------------------------------------------------

def claim_stdio():
    """Take fds 0/1 for the protocol; leave the process a harmless 0/1/stdout.

    Returns ``(reader, writer)`` binary streams on duplicates of the
    original stdin/stdout. fd 0 becomes /dev/null and fd 1 a duplicate of
    stderr, and ``sys.stdout`` follows, so whatever the engine, its tools or
    their children print ends up in the worker's log instead of in the
    middle of a JSON line.
    """
    proto_in = os.fdopen(os.dup(0), "rb", buffering=0)
    proto_out = os.fdopen(os.dup(1), "wb", buffering=0)
    devnull = os.open(os.devnull, os.O_RDONLY)
    os.dup2(devnull, 0)
    os.close(devnull)
    os.dup2(2, 1)
    sys.stdout = sys.stderr
    return proto_in, proto_out


class _Sink:
    """A write-only stream that is not a terminal (TerminalAgent's out/err)."""

    def write(self, _text: str) -> int:
        return 0

    def flush(self) -> None:
        pass

    def isatty(self) -> bool:
        return False


class _NoTerminal:
    """What TerminalAgent sees as stdin: not a terminal, nothing to read."""

    def isatty(self) -> bool:
        return False

    def fileno(self) -> int:
        raise OSError("the worker has no terminal")

    def read(self, *_args) -> str:
        return ""

    def readline(self, *_args) -> str:
        return ""


# ---------------------------------------------------------------------------
# Events
# ---------------------------------------------------------------------------

class Emitter:
    """Stamps every event (session, sequence, time), logs, sends -- in order."""

    def __init__(self, writer: protocol.LineWriter, log: EventLog,
                 session_id: str):
        self.writer = writer
        self.log = log
        self.session_id = session_id
        self._lock = threading.Lock()

    def emit(self, type_: str, **fields: Any) -> None:
        event: dict[str, Any] = {"type": type_, "session": self.session_id,
                                 "ts": round(time.time(), 3)}
        event.update({k: v for k, v in fields.items() if v is not None})
        with self._lock:
            self.log.stamp(event)
            self.writer.send(protocol.event(event))


# ---------------------------------------------------------------------------
# The REPL loop with the terminal swapped out
# ---------------------------------------------------------------------------

def _whole_message(_read_line, first: str, *, cont: str = "") -> str:
    """A prompt from a window is one complete message, whatever it contains."""
    return first


def _agent_class():
    """``TerminalAgent`` subclass, built on first use (the import is heavy)."""
    from .. import repl as _repl
    from ..repl import HarnessText, TerminalAgent

    # The REPL reads a line that is just three quotes as the start of a
    # block and a trailing backslash as "more follows": conventions for
    # typing several lines into a one-line prompt. A window sends each
    # message whole, so in this process -- one session, no terminal -- the
    # continuation grammar is switched off; with it on, a message ending in
    # a backslash would leave the session waiting for input nobody sends.
    _repl._read_block_tail = _whole_message

    class ProtocolAgent(TerminalAgent):
        """TerminalAgent whose keyboard is a queue and whose screen is events."""

        def __init__(self, engine, *, broker, opts, emit, inbox: queue.Queue):
            self._inbox = inbox
            self._emit = emit
            super().__init__(engine, opts, out=_Sink(), err=_Sink(),
                             read_line=self._read_next, broker=broker)
            self.transcript = ProtocolTranscript(emit)
            self._stdin = _NoTerminal()
            self._turn_seq = 0
            self._next_source = "user"
            self._state = ""
            self.closing = threading.Event()

        # -- input: the protocol's queue instead of the keyboard -----------
        def _read_next(self, _prompt: str = "") -> str:
            while True:
                try:
                    item = self._inbox.get(timeout=0.25)
                except queue.Empty:
                    if self.closing.is_set():
                        raise EOFError
                    continue
                if item is _EOF:
                    raise EOFError
                return item

        def _line_ready(self) -> bool:
            return self.closing.is_set() or not self._inbox.empty()

        def _wake_text(self, typed: str) -> str:
            woke = super()._wake_text(typed)
            if woke:
                self._next_source = "harness"
            return woke

        def _show_user_input(self, text: str, *, queued: bool = False) -> None:
            source = "harness" if (queued or self._next_source == "harness") \
                else "user"
            self._next_source = "user"
            cleaned = str(text or "").strip()
            if cleaned:
                self._emit("user.message", text=cleaned, source=source)

        # -- a turn, bracketed by events ------------------------------------
        def turn(self, prompt: str):
            self._turn_seq += 1
            turn_id = f"t{self._turn_seq}"
            source = "harness" if isinstance(prompt, HarnessText) else "user"
            self.transcript.begin_turn(turn_id)
            self._emit("turn.started", turnId=turn_id, source=source)
            self.set_state("running")
            started = time.monotonic()
            before = self._status()
            error = ""
            try:
                result = super().turn(prompt)
                error = str(getattr(result, "error", "") or "")
                return result
            except BaseException as exc:            # noqa: BLE001
                error = str(exc) or type(exc).__name__
                raise
            finally:
                self.transcript.finish()
                after = self._status()
                self._emit(
                    "turn.ended", turnId=turn_id,
                    durationMs=int((time.monotonic() - started) * 1000),
                    inputTokens=max(0, after[0] - before[0]),
                    outputTokens=max(0, after[1] - before[1]),
                    costUsd=round(max(0.0, after[2] - before[2]), 6),
                    toolCalls=self.transcript.tool_calls,
                    error=error or None)
                self.set_state("needs_input" if self.broker and
                               self.broker.pending() else "idle")

        def _status(self) -> tuple[int, int, float]:
            try:
                st = self.engine.get_status() or {}
                return (int(st.get("input_tokens", 0) or 0),
                        int(st.get("output_tokens", 0) or 0),
                        float(st.get("total_cost_usd",
                                     st.get("cost_usd", 0.0)) or 0.0))
            except Exception:
                return (0, 0, 0.0)

        def set_state(self, state: str) -> None:
            if state != self._state:
                self._state = state
                self._emit("session.state", state=state)

        @property
        def turn_running(self) -> bool:
            return self._turn_active.is_set()

        # -- approvals: answered by a window, waited for here ----------------
        def _answer(self, req, raw) -> None:
            """The pump hands us each pending request; the reader thread
            answers it. Wait here so the pump does not spin on it."""
            while not req.event.wait(0.1):
                if self.closing.is_set():
                    return

        # -- steering from the reader thread ---------------------------------
        def steer_text(self, text: str) -> None:
            self._emit("user.message", text=text, source="steer")
            self._steer(text)

        # -- nothing of the terminal remains -----------------------------------
        def _load_history(self) -> None:
            pass

        def _save_history(self) -> None:
            pass

        def _install_completer(self) -> None:
            pass

        def _recall_later(self, text: str) -> None:
            pass

        def _install_sigwinch(self) -> None:
            pass

        def _restore_sigwinch(self) -> None:
            pass

        def _teardown_screen(self) -> None:
            pass

    return ProtocolAgent


# ---------------------------------------------------------------------------
# Opening a session the way `delfin-agent chat` does
# ---------------------------------------------------------------------------

def chat_namespace(args: argparse.Namespace, workspace: Path) -> argparse.Namespace:
    """The namespace `delfin-agent chat` would have, from the worker's args.

    Parsed by the CLI's own parser, so every default the chat command has
    -- and every one upstream adds later -- applies here unchanged.
    """
    from .. import cli as _cli

    argv = ["chat", "--cwd", str(workspace)]
    for flag, value in (("--model", args.model), ("--provider", args.provider),
                        ("--effort", args.effort), ("--mode", args.mode),
                        ("--permission-mode", args.permission_mode),
                        ("-n", args.name), ("--session", args.resume)):
        if value:
            argv += [flag, value]
    if args.continue_session:
        argv.append("--continue")
    if args.fork_session:
        argv.append("--fork-session")
    for path in args.add_dir:
        argv += ["--add-dir", path]
    for path in args.read_dir:
        argv += ["--read-dir", path]
    if args.no_isolate:
        argv.append("--no-isolate")
    elif args.isolate:
        argv.append("--isolate")
    return _cli.build_parser().parse_args(argv)


def open_real_session(args: argparse.Namespace, workspace: Path):
    """Engine, banner lines and launch notices -- cmd_chat's steps, in order.

    Kept in step with ``cli.cmd_chat``: launch checks, isolation override,
    engine, resume/claim, approval posture, banner. A launch finding that
    needs a yes/no is refused here (a worker has nobody to ask yet).
    """
    from .. import cli as _cli
    from .. import launch_guard

    ns = chat_namespace(args, workspace)
    report = launch_guard.inspect_launch_dir(
        workspace, add_dirs=tuple(ns.add_dirs or ()),
        read_dirs=tuple(ns.read_dirs or ()))
    if report.refused:
        raise SessionRefused(report.render())
    if tuple(getattr(report, "questions", ()) or ()):
        raise SessionRefused(
            "This folder needs a yes/no before an agent may start here:\n"
            + report.render())
    os.chdir(workspace)
    ns.extra_dirs = [str(p) for p in report.granted_dirs]
    ns.read_only_dirs = [str(p) for p in report.read_dirs]
    isolation_note = ""
    if getattr(ns, "no_isolate", False):
        from ..api_client import set_bash_isolation_override
        set_bash_isolation_override("off")
        isolation_note = ("isolation  OFF by request (--no-isolate) — a command "
                          "the agent runs can write outside the workspace")
    elif getattr(ns, "isolate", False):
        from ..api_client import set_bash_isolation_override
        set_bash_isolation_override("bwrap")
    engine = _cli._build_engine(ns)
    if not _cli._open_session(engine, ns, workspace):
        raise SessionRefused("This session is being written by another process; "
                             "open a copy with --fork-session or start a new one.")
    persisted, persisted_from = _cli._persisted_default_mode(workspace)
    posture, why = launch_guard.resolve_posture(
        flag_mode=ns.permission_mode or "", persisted_mode=persisted,
        unattended_opt_in=False,
        settings_path=persisted_from or "~/.delfin/settings.json")
    try:
        engine.set_kit_permission_mode(posture)
    except Exception:
        pass
    banner = _cli._startup_banner(engine, report, workspace, why,
                                  isolation_note,
                                  tuple(_cli._bounding_notices(ns, engine)))
    notices = report.render()
    return engine, ns, banner, [notices] if notices else []


def end_real_session(engine, workspace: Path, title: str) -> None:
    """cmd_chat's finally: save, report, learn, index. Never raises."""
    from .. import cli as _cli

    _cli._save_session(engine, workspace, title=title)
    sid = str(getattr(engine, "session_id", "") or "")
    for step in ("report", "learn", "index"):
        try:
            if step == "report":
                from ..session_report import write_session_report
                write_session_report(sid)
            elif step == "learn":
                from ..session_end import learn_at_session_end
                learn_at_session_end(_cli._display_messages(engine), session_id=sid)
            else:
                from ..session_end import index_at_session_end
                index_at_session_end(sid)
        except Exception:
            pass


# ---------------------------------------------------------------------------
# The worker
# ---------------------------------------------------------------------------

class Worker:
    """Requests in, events out, one session in between."""

    def __init__(self, engine, *, writer: protocol.LineWriter, reader,
                 workspace: Path, banner: list[str], notices: list[str],
                 name: str = "", fake: bool = False):
        from ..repl import ReplOptions

        self.engine = engine
        self.writer = writer
        self.reader = reader
        self.workspace = workspace
        self.name = name
        self.fake = fake
        sid = str(getattr(engine, "session_id", "") or "")
        self.log = EventLog(log_path(sid))
        self.emitter = Emitter(writer, self.log, sid)
        emit = self.emitter.emit
        self.inbox: queue.Queue = queue.Queue()
        try:
            from ..cli import _presence_key_for
            key = _presence_key_for(name, sid)
        except Exception:
            key = name or sid
        self.broker = ProtocolBroker(
            emit, stop_turn=self._stop_turn, on_waiting=self._waiting,
            persist=lambda pat: engine.persist_kit_pattern(pat, kind="allow"),
            persist_read=lambda d: engine.add_kit_read_dir(d, persist=True),
            set_mode=engine.set_kit_permission_mode,
            session_id=sid, session_key=key)
        self._bind_broker(key)
        self.agent = _agent_class()(
            engine, broker=self.broker, emit=emit, inbox=self.inbox,
            opts=ReplOptions(cwd=workspace, session_name=name, banner=""))
        self.banner = banner
        self.notices = notices
        self._initialized = False

    def _bind_broker(self, key: str) -> None:
        bind = getattr(self.engine, "set_kit_confirm_callback", None)
        if callable(bind) and bind(self.broker.callback):
            perms = self.engine.kit_permissions
            perms.ask_user_callback = self.broker.ask_user
            perms.plan_approval_callback = self.broker.approve_plan
            perms.presence_key = key

    # -- state -------------------------------------------------------------------
    def _approval(self) -> str:
        try:
            return str(getattr(self.engine.kit_permissions, "mode", "") or "")
        except Exception:
            return ""

    def _waiting(self, waiting: bool) -> None:
        if waiting:
            self.agent.set_state("needs_input")
        else:
            self.agent.set_state("running" if self.agent.turn_running else "idle")

    def _stop_turn(self) -> None:
        self.agent._stop_engine()

    def info(self) -> dict[str, Any]:
        client = getattr(self.engine, "client", None)
        return {"sessionId": self.emitter.session_id,
                "workspace": str(self.workspace),
                "provider": str(getattr(self.engine, "provider", "") or ""),
                "model": str(getattr(client, "model", "") or ""),
                "mode": str(getattr(self.engine, "mode", "") or ""),
                "approval": self._approval(),
                "title": self.name or None,
                "banner": [l for l in self.banner if l is not None],
                "notices": list(self.notices)}

    # -- requests ------------------------------------------------------------------
    def handle(self, message: dict) -> None:
        msg_id = message.get("id")
        method = str(message.get("method") or "")
        params = message.get("params")
        params = {} if params is None else params
        if method not in protocol.REQUESTS:
            if msg_id is not None:
                self.writer.send(protocol.error(
                    msg_id, protocol.METHOD_NOT_FOUND, f"no method {method!r}"))
            return
        problem = protocol.invalid_params(method, params)
        if problem:
            self.writer.send(protocol.error(msg_id, protocol.INVALID_PARAMS,
                                            problem))
            return
        if method != "initialize" and not self._initialized:
            self.writer.send(protocol.error(msg_id, protocol.NOT_ALLOWED,
                                            "initialize first"))
            return
        try:
            value = getattr(self, "_m_" + method.replace(".", "_"))(params)
        except _Refusal as exc:
            self.writer.send(protocol.error(msg_id, exc.code, str(exc)))
            return
        except Exception as exc:                    # noqa: BLE001
            self.writer.send(protocol.error(msg_id, protocol.INTERNAL_ERROR,
                                            f"{type(exc).__name__}: {exc}"))
            return
        if msg_id is not None:
            self.writer.send(protocol.result(msg_id, value))

    def _m_initialize(self, params: dict) -> dict:
        version = int(params.get("protocolVersion", 0) or 0)
        if version != protocol.PROTOCOL_VERSION:
            raise _Refusal(protocol.VERSION_MISMATCH,
                           f"worker speaks protocol {protocol.PROTOCOL_VERSION},"
                           f" client {version}")
        self._initialized = True
        return {"protocolVersion": protocol.PROTOCOL_VERSION, "pid": os.getpid(),
                "session": self.info(),
                "state": self.agent._state or "idle",
                "pending": [self.broker.describe(r) for r in self.broker.pending()]}

    def _m_session_prompt(self, params: dict) -> dict:
        text = str(params.get("text") or "")
        if not text.strip():
            raise _Refusal(protocol.INVALID_PARAMS, "empty prompt")
        how = str(params.get("how") or "queue")
        if self.agent.turn_running:
            if how == "steer":
                self.agent.steer_text(text.strip())
                return {"accepted": True, "steered": True}
            self.agent.queued.append(text)
            self.emitter.emit("user.message", text=text.strip(), source="queued")
            return {"accepted": True, "queued": True}
        self.inbox.put(text)
        return {"accepted": True}

    def _m_session_cancel(self, _params: dict) -> dict:
        refused = self.broker.abort_all() if self.broker.pending() else []
        self.agent._stop_engine()
        return {"stopping": self.agent.turn_running, "refused": len(refused)}

    def _m_session_setMode(self, params: dict) -> dict:
        mode = str(params.get("mode") or "")
        if mode not in ("plan", "default", "acceptEdits"):
            raise _Refusal(protocol.INVALID_PARAMS,
                           "mode must be plan, default or acceptEdits "
                           "(bypassPermissions needs a new session with "
                           "--unattended)")
        self.engine.set_kit_permission_mode(mode)
        self.emitter.emit("session.updated", approval=self._approval())
        return {"approval": self._approval()}

    def _m_session_status(self, _params: dict) -> dict:
        try:
            status = dict(self.engine.get_status() or {})
        except Exception:
            status = {}
        status.update({"state": self.agent._state or "idle",
                       "approval": self._approval(),
                       "queued": len(self.agent.queued)})
        return protocol_safe(status)

    def _m_session_replay(self, params: dict) -> dict:
        turns = int(params.get("turns", 0) or 0)
        return {"events": self.log.replay(turns)}

    def _m_session_shutdown(self, _params: dict) -> dict:
        self.shutdown()
        return {"closing": True}

    def _m_permission_reply(self, params: dict) -> dict:
        from .approvals import InvalidAnswer, UnknownRequest

        try:
            landed, message = self.broker.answer(
                str(params.get("requestId") or ""),
                choice=str(params.get("choice") or ""),
                answers=list(params.get("answers") or []),
                reason=str(params.get("reason") or ""))
        except UnknownRequest as exc:
            raise _Refusal(protocol.UNKNOWN_REQUEST, str(exc)) from None
        except InvalidAnswer as exc:
            raise _Refusal(protocol.INVALID_PARAMS, str(exc)) from None
        return {"applied": landed, "message": message}

    # -- lifecycle -------------------------------------------------------------------
    def shutdown(self) -> None:
        self.agent.closing.set()
        try:
            self.broker.abort_all()
        except Exception:
            pass
        self.agent._stop_engine()
        self.inbox.put(_EOF)

    def _read_requests(self) -> None:
        def bad(problem: str) -> None:
            self.writer.send(protocol.error(None, protocol.PARSE_ERROR, problem))

        try:
            for message in protocol.read_messages(self.reader, on_bad_line=bad):
                self.handle(message)
        finally:
            # The window side went away (EOF): end the session cleanly.
            self.shutdown()

    def _sample_usage(self) -> None:
        last = None
        while not self.agent.closing.is_set():
            time.sleep(_USAGE_TICK_S)
            if not self.agent.turn_running:
                continue
            try:
                st = self.engine.get_status() or {}
            except Exception:
                continue
            sample = (int(st.get("input_tokens", 0) or 0),
                      int(st.get("output_tokens_live",
                                 st.get("output_tokens", 0)) or 0),
                      round(float(st.get("total_cost_usd",
                                         st.get("cost_usd", 0.0)) or 0.0), 6))
            if sample != last:
                last = sample
                self.emitter.emit("usage", inputTokens=sample[0],
                                  outputTokens=sample[1], costUsd=sample[2])

    def run(self) -> int:
        self.emitter.emit("worker.ready", protocolVersion=protocol.PROTOCOL_VERSION,
                          pid=os.getpid())
        self.emitter.emit("session.info", **self.info())
        self.agent.set_state("idle")
        threading.Thread(target=self._read_requests, name="protocol-reader",
                         daemon=True).start()
        threading.Thread(target=self._sample_usage, name="usage-sampler",
                         daemon=True).start()
        code = 0
        try:
            code = int(self.agent.run() or 0)
        finally:
            self.agent.closing.set()
            if not self.fake:
                end_real_session(self.engine, self.workspace, self.name)
            self.agent.set_state("ended")
            self.emitter.emit("session.ended", reason="closed", exitCode=code)
        return code


class _Refusal(Exception):
    def __init__(self, code: int, message: str):
        super().__init__(message)
        self.code = code


def protocol_safe(value: Any) -> Any:
    import json
    return json.loads(json.dumps(value, default=str))


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(prog="python -m delfin.agent.service.worker",
                                description="One agent session over the TUI protocol")
    p.add_argument("--cwd", default="")
    p.add_argument("--resume", default="", help="session id to resume")
    p.add_argument("--continue", dest="continue_session", action="store_true")
    p.add_argument("--fork-session", action="store_true")
    p.add_argument("-n", "--name", default="")
    p.add_argument("--model", default="")
    p.add_argument("--provider", default="")
    p.add_argument("--effort", default="")
    p.add_argument("--mode", default="")
    p.add_argument("--permission-mode", default="",
                   choices=["", "plan", "default", "acceptEdits"])
    p.add_argument("--add-dir", action="append", default=[])
    p.add_argument("--read-dir", action="append", default=[])
    p.add_argument("--isolate", action="store_true")
    p.add_argument("--no-isolate", action="store_true")
    p.add_argument("--fake", action="store_true",
                   help="scripted turns, no model (tests and UI work)")
    p.add_argument("--fake-delay", type=float, default=0.02)
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    from .. import state_paths

    # Before anything imports a state sink, as cli.main does.
    state_paths.scratch_state_from_environment()
    args = parse_args(argv)
    reader, raw_out = claim_stdio()
    writer = protocol.LineWriter(raw_out)
    workspace = Path(args.cwd or os.getcwd()).expanduser().resolve()
    try:
        if args.fake:
            from .fake_engine import FakeEngine
            engine = FakeEngine(str(workspace), delay=args.fake_delay)
            if args.permission_mode:
                engine.set_kit_permission_mode(args.permission_mode)
            banner = [f"delfin-agent · fake/fake-model · solo",
                      f"workspace  {workspace}",
                      f"approval   {engine.kit_permissions.mode}"]
            notices: list[str] = []
        else:
            engine, _ns, banner_text, notices = open_real_session(args, workspace)
            banner = banner_text.splitlines()
    except SessionRefused as exc:
        writer.send(protocol.event({"type": "session.ended", "session": "",
                                    "reason": str(exc), "exitCode": 2}))
        return 2
    worker = Worker(engine, writer=writer, reader=reader, workspace=workspace,
                    banner=banner, notices=notices, name=args.name,
                    fake=args.fake)
    return worker.run()


if __name__ == "__main__":
    sys.exit(main())
