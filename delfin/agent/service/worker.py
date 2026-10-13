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

Three threads share the session:

* the main thread runs the REPL loop; during a turn it is the REPL's pump,
  turning the engine's output into events in the order it was produced;
* the reader thread handles requests as they arrive -- a prompt is routed
  by what the loop is doing at that moment (``ProtocolAgent.lock``);
* the look-around thread notices what changes on its own: usage, the
  approval posture, mail for a running turn, a cancel that must hold.

Anything said from another thread that belongs in the stream -- an
approval asked or settled, a message queued -- goes through
``ProtocolAgent.post``, which puts it behind what the turn has already
produced, so a window sees things in the order they happened.

Ending a session goes the REPL's way out: ``session.shutdown`` ends the
loop as ctrl+d at the prompt does; a window that goes away (EOF) or a
SIGTERM ends it the way a closed terminal does -- no further turn starts,
every open question is refused, background shells and MCP servers end,
everything within the REPL's leave deadline -- and the session is saved.

stdout carries protocol messages and nothing else: file descriptors 0 and 1
are taken over at start-up, so neither a stray print nor a child process
can write into the protocol stream.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import os
import sys
import threading
import time
from collections import deque
from pathlib import Path
from typing import Any, Callable

from . import choices, protocol
from .approvals import ProtocolBroker
from .eventlog import REPLAY_BUDGET_BYTES, EventLog, log_path, service_dir
from .transcript import ProtocolTranscript

#: How often the look-around thread looks. Also how long a cancel that
#: came before the engine's turn began can go unnoticed: the engine clears
#: its stop flag when a turn starts, so the cancel is repeated each tick.
_TICK_S = 0.25

#: How long the idle loop waits for a message before it looks for mail
#: and finished jobs (the REPL throttles those to every few seconds).
_IDLE_WAIT_S = 0.5

#: Requests answered on a thread of their own: a provider's model list
#: can take seconds, and the reader must stay free for the next request --
#: an Esc never waits behind a picker.
_SLOW = frozenset({"session.choices"})

#: Live-only events that are not redacted: the model's own stream, which
#: the REPL shows as it arrives too, and numbers. Everything else -- what
#: is stored and what a window shows from a tool -- goes through the
#: secret redactor.
_UNREDACTED = frozenset({"text.delta", "reasoning.delta", "usage",
                         "session.state"})

#: DELFIN_TUI_TRACE=1: every request received, on stderr (the log).
_TRACE = os.environ.get("DELFIN_TUI_TRACE", "").strip() not in ("", "0")

#: Fields that identify things; never rewritten by the redactor.
_IDENTIFIERS = frozenset({"type", "session", "seq", "ts", "turnId", "partId",
                          "callId", "requestId", "sessionId", "source",
                          "kind", "state", "decision", "by", "flags"})


class SessionRefused(Exception):
    """The session cannot be opened; the message says why and what to do."""


# ---------------------------------------------------------------------------
# stdio
# ---------------------------------------------------------------------------

def _same_file(fd_a: int, fd_b: int) -> bool:
    try:
        a, b = os.fstat(fd_a), os.fstat(fd_b)
    except OSError:
        return False
    return (a.st_dev, a.st_ino) == (b.st_dev, b.st_ino)


def claim_stdio():
    """Take fds 0/1 for the protocol; leave the process a harmless 0/1/2.

    Returns ``(reader, writer)`` binary streams on duplicates of the
    original stdin/stdout. fd 0 becomes /dev/null; fd 1 points where fd 2
    points, and ``sys.stdout`` follows, so whatever the engine, its tools or
    their children print ends up in the worker's log instead of in the
    middle of a JSON line. When the caller handed the protocol pipe in as
    stderr too, stderr goes to a log file instead -- never into the stream.
    """
    # Buffered for reading: a long request line is read in blocks, not
    # byte by byte. Unbuffered for writing: LineWriter writes whole lines.
    proto_in = os.fdopen(os.dup(0), "rb")
    proto_out = os.fdopen(os.dup(1), "wb", buffering=0)
    devnull = os.open(os.devnull, os.O_RDONLY)
    os.dup2(devnull, 0)
    os.close(devnull)
    if _same_file(1, 2):
        from .. import state_paths
        logs = service_dir() / "logs"
        state_paths.ensure_dir(logs)
        log = os.open(logs / f"worker-{os.getpid()}.log",
                      os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o600)
        os.dup2(log, 2)
        os.close(log)
        sys.stderr = os.fdopen(os.dup(2), "w", buffering=1)
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

class Redactor:
    """Masks credentials in what a session shows and keeps.

    Two passes, as ``output_guard.scrub_secrets`` does: the exact values of
    the credentials DELFIN holds (a KIT key has no recognisable shape), then
    everything shaped like a credential. The known values are read once a
    minute, not once per string.
    """

    _REFRESH_S = 60.0

    def __init__(self):
        self._known: list[str] = []
        self._loaded = -1e9
        try:
            from ..output_guard import _redact_secrets
            self._shapes: Callable[[str, list], str] | None = _redact_secrets
        except Exception:
            self._shapes = None

    def _known_values(self) -> list[str]:
        now = time.monotonic()
        if now - self._loaded > self._REFRESH_S:
            self._loaded = now
            try:
                from ..web_tools import _known_secret_values
                self._known = sorted(set(_known_secret_values()), key=len,
                                     reverse=True)
            except Exception:
                self._known = []
        return self._known

    def text(self, text: str) -> str:
        for secret in self._known_values():
            if secret in text:
                text = text.replace(secret, "[redacted:credential]")
        if self._shapes is not None:
            try:
                text = self._shapes(text, [])
            except Exception:
                pass
        return text

    def value(self, value: Any) -> Any:
        if isinstance(value, str):
            return self.text(value)
        if isinstance(value, dict):
            return {k: self.value(v) for k, v in value.items()}
        if isinstance(value, list):
            return [self.value(v) for v in value]
        return value

    def event(self, event: dict) -> dict:
        return {k: (v if k in _IDENTIFIERS else self.value(v))
                for k, v in event.items()}


class Emitter:
    """Stamps every event (session, time, and a sequence number for the
    durable ones), masks secrets, logs, sends -- in order. Never raises:
    an event that cannot be sent must not end the turn that produced it."""

    def __init__(self, writer: protocol.LineWriter, log: EventLog,
                 session_id: str):
        self.writer = writer
        self.log = log
        self.session_id = session_id
        self.redactor = Redactor()
        self._lock = threading.Lock()

    def emit(self, type_: str, **fields: Any) -> None:
        try:
            event: dict[str, Any] = {"type": type_, "session": self.session_id,
                                     "ts": round(time.time(), 3)}
            event.update({k: v for k, v in fields.items() if v is not None})
            if type_ not in _UNREDACTED:
                event = self.redactor.event(event)
            with self._lock:
                self.log.stamp(event)
                self.writer.send(protocol.event(event))
        except Exception as exc:                    # noqa: BLE001
            try:
                print(f"worker: event {type_} not sent: {exc}", file=sys.stderr)
            except Exception:
                pass


# ---------------------------------------------------------------------------
# The REPL loop with the terminal swapped out
# ---------------------------------------------------------------------------

def _whole_message(_read_line, first: str, *, cont: str = "") -> str:
    """A prompt from a window is one complete message, whatever it contains."""
    return first


def _agent_class():
    """``TerminalAgent`` subclass, built on first use (the import is heavy)."""
    import signal

    from .. import repl_keys as rk
    from ..repl import HarnessText, RenderItem, TerminalAgent, TurnResult

    class _Lines(list):
        """``queued`` as the REPL loop takes from it: each line it takes
        is the next line in progress (what a cancel can refer to)."""

        def __init__(self, on_take: Callable[[], None]):
            super().__init__()
            self._on_take = on_take

        def pop(self, index=-1):
            item = super().pop(index)
            self._on_take()
            return item

    class ProtocolAgent(TerminalAgent):
        """TerminalAgent whose keyboard is a request and whose screen is
        events.

        What the loop is doing decides where a prompt goes, so the two are
        kept under one lock (``lock``/``cond``):

        * ``idle``  -- the loop waits at the prompt: a prompt goes to the
          ``inbox`` and is read as the next line;
        * ``busy``  -- inside ``turn()``: a prompt is queued behind the turn
          (or steered into it while the engine's turn is live);
        * neither   -- between: a line is being handled, or a turn just
          ended and its follow-up is being decided. A prompt is queued,
          which also drops the automatic follow-up, as a line typed
          during the turn does in the REPL.

        Every prompt is shown (``user.message``) when it arrives, so the
        window sees messages in the order they were sent; lines read from
        the inbox are therefore not shown again when the loop takes them.
        """

        def __init__(self, engine, *, broker, opts, emit, checkpoint=True):
            self._emit = emit
            super().__init__(engine, opts, out=_Sink(), err=_Sink(),
                             read_line=self._read_next, broker=broker)
            self.transcript = ProtocolTranscript(emit)
            self._stdin = _NoTerminal()
            self.lock = threading.Lock()
            self.cond = threading.Condition(self.lock)
            self.queued = _Lines(self._took_line)
            #: Lines for the idle read, oldest first.
            self.inbox: deque[str] = deque()
            self.idle = True
            self.busy = False
            #: Which line the loop is working on; a cancel names one.
            self.line_gen = 0
            #: Esc: ends the turn that runs, or the one about to start for
            #: the line ``cancel_gen``; never a line sent after it.
            self.cancel_requested = False
            self.cancel_gen = -1
            #: No automatic follow-up (length, outage, announced work)
            #: until the user sends something again.
            self.followups_off = False
            self.closing = threading.Event()
            self.leave_reason = ""
            self.leave_code = 0
            #: Set by a prompt, a cancel or shutdown: ends the wait during
            #: an endpoint outage the way a key press ends it in the REPL.
            self.wake = threading.Event()
            self._turn_seq = 0
            self._next_source = "user"
            #: The line just read was shown when it arrived.
            self._shown_already = False
            #: The next turn is a wake-up (a finished job, a message), not
            #: a follow-up: a queued line does not drop it.
            self._wake_turn = False
            #: Steers routed into the running turn (to take back unread).
            self._steered = 0
            self._checkpoint = checkpoint
            self._post_lock = threading.Lock()
            self._pumping = False
            self._state_lock = threading.Lock()
            self._state = ""
            self.ended = False

        def _leave(self) -> BaseException:
            if self.leave_reason == "shutdown":
                # Asked to end: the way ctrl+d at the prompt ends a session.
                return EOFError()
            # The window is gone, or a signal: the REPL's way out when its
            # terminal is gone -- run() stops the turn, refuses every
            # question, ends what the session started, under the deadline.
            return rk.TerminalLeft(self.leave_reason or "closed",
                                   self.leave_code or 129)

        def leave_now(self) -> None:
            """Unwind the loop wherever it waits -- a long tool call, a
            ``!command`` -- the way a hangup does in the REPL: its SIGHUP
            handler raises TerminalLeft on the main thread, and the leave
            runs under the REPL's deadline. Only while that handler is in
            place; outside the loop the worker's own handler takes it."""
            try:
                if signal.getsignal(signal.SIGHUP) == self._on_leave_signal:
                    signal.pthread_kill(threading.main_thread().ident,
                                        signal.SIGHUP)
            except (AttributeError, ValueError, OSError):
                pass

        # -- input: requests instead of the keyboard -------------------------
        def _took_line(self) -> None:
            with self.cond:
                self.line_gen += 1

        def _line_ready(self) -> bool:
            with self.cond:
                if self.queued and not self.closing.is_set():
                    # Queued after the loop last looked for queued lines and
                    # before it came here: they are older than anything the
                    # inbox holds, and were shown when they arrived.
                    stranded = list(self.queued)
                    del self.queued[:]
                    self.inbox.extendleft(reversed(stranded))
                if self.closing.is_set() or self.inbox:
                    return True
                came_back = not self.idle
                self.idle = True
                self.cancel_requested = False
                self.cancel_gen = -1
                self.followups_off = False
            if came_back:
                self._back_at_the_prompt()
                self.refresh_state()
            with self.cond:
                if not (self.closing.is_set() or self.inbox):
                    self.cond.wait(_IDLE_WAIT_S)
                return self.closing.is_set() or bool(self.inbox)

        def _back_at_the_prompt(self) -> None:
            """What a cancel armed for the line before ends here: the
            refusal of questions, and the engine's stop request (a stop
            left armed ends the next ``!command`` with 130)."""
            if self.closing.is_set() or self.broker is None:
                return
            try:
                self.broker.reset_abort()
            except Exception:
                pass
            try:
                self.engine.clear_stop()
            except Exception:
                pass

        def _read_next(self, _prompt: str = "") -> str:
            with self.cond:
                while not self.closing.is_set() and not self.inbox:
                    self.idle = True
                    self.cond.wait(_IDLE_WAIT_S)
                if self.closing.is_set():
                    raise self._leave()
                text = self.inbox.popleft()
                # The rest waits behind it, in order, like lines typed
                # during the turn it starts.
                self.queued.extend(self.inbox)
                self.inbox.clear()
                self._shown_already = True
                self.line_gen += 1
                self.idle = False
            self.refresh_state()
            return text

        def _wake_text(self, typed: str) -> str:
            woke = super()._wake_text(typed)
            if woke:
                with self.cond:
                    self.idle = False
                    self.line_gen += 1
                    self._wake_turn = True
                    # A prompt that arrived while the wake was decided goes
                    # after the turn the wake starts (it was shown already).
                    self.queued.extend(self.inbox)
                    self.inbox.clear()
                self._next_source = "harness"
                self.refresh_state()
            return woke

        def _handle_line(self, text: str):
            if self.closing.is_set():
                # Lines still queued when the session ends are dropped,
                # never run: a queued "!command" must not run on the way out.
                return ""
            # A command's output, printed line by line, is one notice.
            with self.transcript.block():
                return super()._handle_line(text)

        def _show_user_input(self, text: str, *, queued: bool = False) -> None:
            if self._shown_already:
                # A line from the inbox: shown when it arrived.
                self._shown_already = False
                return
            if queued and (self.queued or self.inbox):
                # A follow-up that a line sent meanwhile drops (turn()).
                return
            source = "harness" if (queued or self._next_source == "harness") \
                else "user"
            self._next_source = "user"
            cleaned = str(text or "").strip()
            if cleaned:
                self._emit("user.message", text=cleaned, source=source)

        def _checkpoint_session(self) -> None:
            if self._checkpoint:
                super()._checkpoint_session()

        # -- the stream: the turn's output and what other threads say --------
        def post(self, fn: Callable[[], None]) -> None:
            """Run ``fn`` in stream order: behind everything the running
            turn has already produced, on the thread that renders it; at
            once when no turn is running."""
            with self._post_lock:
                if self._pumping:
                    self._q.put(RenderItem("call", data={"fn": fn}))
                    return
            fn()

        def _own_stream(self) -> None:
            with self._post_lock:
                self._pumping = True

        def _release_stream(self) -> None:
            """Render what is still queued, then let ``post`` run directly.

            Item by item outside the lock (rendering may take the broker's
            lock, which an asking thread holds while it posts), and the
            switch happens only once the queue is empty under the lock.
            """
            while True:
                with self._post_lock:
                    try:
                        item = self._q.get_nowait()
                    except Exception:
                        self._pumping = False
                        return
                if item.kind != "done":
                    try:
                        self.transcript.render(item)
                    except Exception:
                        pass

        # -- a turn, bracketed by events ------------------------------------
        def _begin_turn_resets_the_abort(self) -> None:
            # The REPL lifts an abort at every turn start. Not for a turn
            # that starts while the session closes or is being cancelled:
            # the refusal that shutdown or Esc just put in place must hold,
            # or the turn's question waits with nobody to answer it.
            # Checked and reset under the lock a cancel sets its flag and
            # refuses under: either the reset comes first and the cancel's
            # refusal follows it, or the cancel comes first and there is no
            # reset. (The broker never takes this lock, so the order holds.)
            with self.cond:
                if self.closing.is_set() or self.cancel_requested:
                    return
                super()._begin_turn_resets_the_abort()

        def turn(self, prompt: str):
            harness = isinstance(prompt, HarnessText)
            with self.cond:
                if self.closing.is_set():
                    return TurnResult(error="the session is closing")
                wake, self._wake_turn = self._wake_turn, False
                if harness and not wake and (self.queued or self.inbox):
                    # A follow-up: a line sent while it was being set up goes
                    # first, and the follow-up is dropped -- as the REPL drops
                    # it for a line typed during the turn.
                    return TurnResult()
                self.busy = True
                self.idle = False
                skip = self.cancel_requested and self.cancel_gen == self.line_gen
                if skip:
                    self.followups_off = True
                elif not harness:
                    # Something the user sent: its own follow-ups may run.
                    self.followups_off = False
                self._steered = 0
            self.wake.clear()
            self._own_stream()
            self._turn_seq += 1
            turn_id = f"t{self._turn_seq}"
            self.transcript.begin_turn(turn_id)
            self._emit("turn.started", turnId=turn_id,
                       source="harness" if harness else "user")
            self.refresh_state()
            started = time.monotonic()
            before = self._status()
            error = ""
            try:
                if skip:
                    # Esc right after Enter: the turn never starts.
                    error = "cancelled before it started"
                    return TurnResult(error=error)
                result = super().turn(prompt)
                error = str(getattr(result, "error", "") or "")
                return result
            except BaseException as exc:            # noqa: BLE001
                error = str(exc) or type(exc).__name__
                raise
            finally:
                self._release_stream()
                self.transcript.finish()
                late = self._take_back_steers()
                after = self._status()
                self._emit(
                    "turn.ended", turnId=turn_id,
                    durationMs=int((time.monotonic() - started) * 1000),
                    inputTokens=max(0, after[0] - before[0]),
                    outputTokens=max(0, after[1] - before[1]),
                    costUsd=round(max(0.0, after[2] - before[2]), 6),
                    toolCalls=self.transcript.tool_calls,
                    error=error or None)
                if late:
                    self._emit("notice", text=(
                        "A message sent into the turn arrived after it had "
                        "ended; it goes out next."))
                with self.cond:
                    self.busy = False
                    self.cancel_requested = False
                    if late and not self.closing.is_set():
                        self.queued[:0] = late
                self.refresh_state()

        def _take_back_steers(self) -> list[str]:
            """Messages steered into this turn that it never read: they would
            surface in the next, unrelated turn. Plain text only -- a
            session's mail (HarnessText) is the REPL's to deliver."""
            if not self._steered:
                return []
            out: list[str] = []
            client = getattr(self.engine, "client", None)
            drain = getattr(client, "_drain_steer", None)
            try:
                if callable(drain):
                    leftovers = list(drain(delivered=False) or ())
                    for text in leftovers:
                        if isinstance(text, HarnessText):
                            client.push_steer(text)
                        else:
                            out.append(str(text))
                else:
                    fallback = getattr(self.engine, "drain_undelivered_steers", None)
                    if callable(fallback):
                        out = [str(t) for t in fallback() or ()]
            except Exception:
                return []
            return [t for t in out if t.strip()]

        def _status(self) -> tuple[int, int, float]:
            try:
                st = self.engine.get_status() or {}
                return (int(st.get("input_tokens", 0) or 0),
                        int(st.get("output_tokens", 0) or 0),
                        float(st.get("total_cost_usd",
                                     st.get("cost_usd", 0.0)) or 0.0))
            except Exception:
                return (0, 0, 0.0)

        def refresh_state(self) -> None:
            """The state from the facts, announced when it changed."""
            with self._state_lock:
                if self.ended:
                    state = "ended"
                elif self.broker is not None and self.broker.pending():
                    state = "needs_input"
                elif not self.idle:
                    state = "running"
                else:
                    state = "idle"
                if state != self._state:
                    self._state = state
                    self._emit("session.state", state=state)

        @property
        def state(self) -> str:
            return self._state or "idle"

        # -- follow-ups: none once the user said stop or sent something, or
        # -- the session ends ------------------------------------------------
        def _follow_ups_off(self) -> bool:
            return (self.closing.is_set() or self.followups_off
                    or bool(self.queued) or bool(self.inbox))

        def _continue_after_length(self, answered_prompt: str) -> str:
            if self._follow_ups_off():
                return ""
            return super()._continue_after_length(answered_prompt)

        def _continue_after_endpoint_failure(self) -> str:
            if self._follow_ups_off():
                self.__dict__["_endpoint_pauses"] = 0
                return ""
            return super()._continue_after_endpoint_failure()

        def _continue_after_announcement(self) -> str:
            if self._follow_ups_off():
                try:
                    self.engine.clear_turn_continuation()
                except Exception:
                    pass
                return ""
            return super()._continue_after_announcement()

        def _endpoint_read_key(self, timeout: float):
            """The outage pause: a prompt, a cancel or shutdown ends it, as
            any key does in the REPL; otherwise it runs its full length."""
            deadline = time.monotonic() + max(0.0, float(timeout))
            while True:
                with self.cond:
                    interrupted = self._follow_ups_off()
                if interrupted:
                    return "\x1b"
                left = deadline - time.monotonic()
                if left <= 0:
                    return ""
                self.wake.wait(min(left, _IDLE_WAIT_S))
                self.wake.clear()

        # -- steering from the reader thread ---------------------------------
        def steer_text(self, text: str) -> bool:
            """Into the running turn; False when the backend takes nothing
            mid-turn (the caller queues it then, as the REPL does)."""
            try:
                delivered = bool(self.engine.steer(text))
            except Exception:
                delivered = False
            if delivered:
                with self.cond:
                    self._steered += 1
                self.post(lambda: self._emit("user.message", text=text,
                                             source="steer"))
            return delivered

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


def _lines(text: str) -> list[str]:
    return [line for line in text.splitlines() if line.strip()]


def open_real_session(args: argparse.Namespace, workspace: Path):
    """Engine, banner lines and opening notes -- cmd_chat's steps, in order.

    Kept in step with ``cli.cmd_chat``: launch checks, isolation override,
    engine, resume/claim, approval posture, banner. A launch finding that
    needs a yes/no is refused here (a worker has nobody to ask yet). What
    those steps print for the person -- "session not found, starting
    fresh", a fork, a crash-recovery note, a partial restore, isolation
    that was asked for and is not available -- is captured and handed to
    the window instead of a log nobody reads.
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
    # cmd_chat's own steps, by its own helpers: the grants launch_guard
    # accepted, then the isolation flags (with the line that says so).
    _cli._apply_launch_grants(ns, report)
    said = io.StringIO()
    with contextlib.redirect_stderr(said):
        isolation_note = _cli._apply_isolation_flags(ns)
        try:
            engine = _cli._build_engine(ns)
        except SystemExit as exc:
            raise SessionRefused(said.getvalue().strip() or
                                 f"The agent engine could not start ({exc.code})."
                                 ) from None
        except Exception as exc:
            raise SessionRefused(f"The agent engine could not start: {exc}") from None
        if not _cli._open_session(engine, ns, workspace):
            raise SessionRefused(said.getvalue().strip() or
                                 "This session is being written by another "
                                 "process; open a copy with --fork-session or "
                                 "start a new one.")
        persisted, persisted_from = _cli._persisted_default_mode(workspace)
        posture, why = launch_guard.resolve_posture(
            flag_mode=ns.permission_mode or "", persisted_mode=persisted,
            unattended_opt_in=False,
            settings_path=persisted_from or "~/.delfin/settings.json")
        try:
            # The starting rule, not a switch: no turn is running to be told
            # (cmd_chat's call, with its fallback for an older engine).
            try:
                engine.set_kit_permission_mode(posture, announce=False)
            except TypeError:
                engine.set_kit_permission_mode(posture)
        except Exception:
            pass
        banner = _cli._startup_banner(engine, report, workspace, why,
                                      isolation_note,
                                      tuple(_cli._bounding_notices(ns, engine)))
    notices = _lines(said.getvalue())
    launch = report.render()
    if launch:
        notices.append(launch)
    return engine, ns, banner, notices


def end_real_session(engine, workspace: Path, title: str, *,
                     learn: bool = True) -> None:
    """cmd_chat's finally: save, report, learn, index. Never raises."""
    from .. import cli as _cli

    try:
        _cli._save_session(engine, workspace, title=title)
    except Exception:
        pass
    if not learn:
        return
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
                 name: str = "", fake: bool = False, save: bool = True):
        from ..repl import ReplOptions

        self.engine = engine
        self.writer = writer
        self.reader = reader
        self.workspace = workspace
        self.name = name
        self.fake = fake
        self.save = save
        sid = str(getattr(engine, "session_id", "") or "")
        self.notices = list(notices)
        self.log = EventLog(log_path(sid), on_error=self._log_failed)
        self.emitter = Emitter(writer, self.log, sid)
        emit = self.emitter.emit
        try:
            from ..cli import _presence_key_for
            key = _presence_key_for(name, sid)
        except Exception:
            key = name or sid
        self.broker = ProtocolBroker(
            self._emit_in_order, stop_turn=self._stop_turn,
            on_waiting=self._waiting, on_feedback=self._feedback,
            persist=lambda pat: engine.persist_kit_pattern(pat, kind="allow"),
            persist_read=lambda d: engine.add_kit_read_dir(d, persist=True),
            set_mode=engine.set_kit_permission_mode,
            session_id=sid, session_key=key)
        self._bind_broker(key)
        self.agent = _agent_class()(
            engine, broker=self.broker, emit=emit, checkpoint=save,
            opts=ReplOptions(cwd=workspace, session_name=name, banner=""))
        # Turn ids stay unique in the session's log across resumes.
        self.agent._turn_seq = self.log.turns
        self.banner = banner
        self._initialized = False
        self._seen_lock = threading.Lock()
        self._seen = self._session_facts()
        # Held while a request is being handled: the process does not end
        # under a request whose answer is still to be sent.
        self._handling = threading.Lock()

    def _log_failed(self, exc: Exception) -> None:
        self.notices.append(f"The session's event log cannot be written ({exc}); "
                            "a window that attaches later will not see "
                            "earlier turns.")

    def _bind_broker(self, key: str) -> None:
        bind = getattr(self.engine, "set_kit_confirm_callback", None)
        if callable(bind) and bind(self.broker.callback):
            perms = self.engine.kit_permissions
            perms.ask_user_callback = self.broker.ask_user
            perms.plan_approval_callback = self.broker.approve_plan
            perms.presence_key = key

    def _emit_in_order(self, type_: str, **fields: Any) -> None:
        """What the broker announces, in stream order (``ProtocolAgent.post``)."""
        def say() -> None:
            if type_ == "permission.asked":
                # Drawn below what was said before it, as a tool call is.
                self.agent.transcript.close_parts()
            self.emitter.emit(type_, **fields)

        agent = getattr(self, "agent", None)
        if agent is None:
            self.emitter.emit(type_, **fields)
        else:
            agent.post(say)

    # -- state -------------------------------------------------------------------
    def _session_facts(self) -> dict[str, str]:
        client = getattr(self.engine, "client", None)
        try:
            approval = str(getattr(self.engine.kit_permissions, "mode", "") or "")
        except Exception:
            approval = ""
        return {"approval": approval,
                "model": str(getattr(client, "model", "") or ""),
                "effort": str(getattr(client, "effort", "")
                              or getattr(self.engine, "effort", "") or "")}

    def check_session(self) -> None:
        """Announce the posture, model or effort whenever it changed,
        however it changed: a plan approval, "stop asking about writes",
        /permissions, /model, setMode."""
        with self._seen_lock:
            now = self._session_facts()
            changed = {k: v for k, v in now.items() if v != self._seen.get(k)}
            self._seen = now
        if changed:
            self.emitter.emit("session.updated", **changed)

    def _waiting(self, _waiting: bool) -> None:
        self.agent.post(self.agent.refresh_state)

    def _stop_turn(self) -> None:
        self.agent._stop_engine()

    def _feedback(self, text: str) -> None:
        if not self.agent.steer_text(text):
            self.agent.post(lambda: self.emitter.emit("notice", text=(
                "This backend takes no message mid-turn; what you said "
                "about the plan was not delivered.")))

    def info(self) -> dict[str, Any]:
        facts = self._session_facts()
        return {"sessionId": self.emitter.session_id,
                "workspace": str(self.workspace),
                "provider": str(getattr(self.engine, "provider", "") or ""),
                "model": facts["model"],
                "mode": str(getattr(self.engine, "mode", "") or ""),
                "approval": facts["approval"],
                "title": self.name or None,
                "banner": [l for l in self.banner if l is not None],
                "notices": list(self.notices)}

    def header(self) -> dict[str, Any]:
        """A ``session.info`` event as the session is now (for replays)."""
        event = {"type": "session.info", "session": self.emitter.session_id,
                 "seq": 0, "ts": round(time.time(), 3)}
        event.update({k: v for k, v in self.info().items() if v is not None})
        return self.emitter.redactor.event(event)

    # -- requests ------------------------------------------------------------------
    def handle(self, message: dict) -> None:
        msg_id = message.get("id")
        method = str(message.get("method") or "")
        if _TRACE:
            print(f"worker: {time.time():.3f} <- {method} #{msg_id}",
                  file=sys.stderr, flush=True)
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
        # A shutdown is taken at any time: the service ends a worker that is
        # still starting the same way as one that runs.
        if method not in ("initialize", "session.shutdown") and not self._initialized:
            self.writer.send(protocol.error(msg_id, protocol.NOT_ALLOWED,
                                            "initialize first"))
            return
        if method == "session.shutdown":
            # Answered first: once the session is told to close, the
            # process may be gone before a later answer is written.
            if msg_id is not None:
                self.writer.send(protocol.result(msg_id, {"closing": True}))
            self.shutdown("shutdown")
            return
        if method in _SLOW:
            threading.Thread(target=self._answer, args=(msg_id, method, params),
                             name=f"protocol-{method}", daemon=True).start()
            return
        self._answer(msg_id, method, params)

    def _answer(self, msg_id: Any, method: str, params: dict) -> None:
        # Answers carry what events carry -- a pending question's preview,
        # a persisted pattern, an error naming the command -- so they go
        # through the same redactor.
        clean = self.emitter.redactor.value
        try:
            value = getattr(self, "_m_" + method.replace(".", "_"))(params)
        except _Refusal as exc:
            self.writer.send(protocol.error(msg_id, exc.code, clean(str(exc))))
            return
        except Exception as exc:                    # noqa: BLE001
            self.writer.send(protocol.error(
                msg_id, protocol.INTERNAL_ERROR,
                clean(f"{type(exc).__name__}: {exc}")))
            return
        if msg_id is not None:
            self.writer.send(protocol.result(msg_id, clean(value)))

    def _m_initialize(self, params: dict) -> dict:
        version = int(params.get("protocolVersion", 0) or 0)
        if version != protocol.PROTOCOL_VERSION:
            raise _Refusal(protocol.VERSION_MISMATCH,
                           f"worker speaks protocol {protocol.PROTOCOL_VERSION},"
                           f" client {version}")
        self._initialized = True
        return {"protocolVersion": protocol.PROTOCOL_VERSION, "pid": os.getpid(),
                "session": self.info(), "state": self.agent.state,
                "pending": [self.broker.describe(r) for r in self.broker.pending()]}

    def _m_session_prompt(self, params: dict) -> dict:
        # Stripped as the REPL strips a typed line: "  !ls" is a shell
        # escape at the prompt and in the queue alike.
        text = str(params.get("text") or "").strip()
        if not text:
            raise _Refusal(protocol.INVALID_PARAMS, "empty prompt")
        how = str(params.get("how") or "queue")
        if how not in ("queue", "steer"):
            raise _Refusal(protocol.INVALID_PARAMS, "how must be queue or steer")
        agent = self.agent
        with agent.cond:
            if agent.closing.is_set():
                raise _Refusal(protocol.NOT_ALLOWED, "the session is closing")
            # Steered only while the engine's turn is live: once it has
            # finished, a steer would surface in the next, unrelated turn.
            if how == "steer" and agent.busy and agent._turn_active.is_set():
                route = "steer"
            elif agent.idle:
                agent.inbox.append(text)
                agent.cond.notify_all()
                route = "next"
                # Shown now, under the lock: before the loop can take it and
                # announce the turn it starts.
                self.emitter.emit("user.message", text=text, source="user")
            else:
                # During a turn, or between a turn and its follow-up: like
                # a line typed during the turn, it goes before any
                # follow-up, which is dropped.
                agent.queued.append(text)
                route = "queued"
        agent.wake.set()
        if route == "steer":
            if agent.steer_text(text):
                return {"accepted": True, "steered": True}
            with agent.cond:
                agent.queued.append(text)
            route = "queued"
        if route == "queued":
            agent.post(lambda: self.emitter.emit(
                "user.message", text=text, source="queued"))
            return {"accepted": True, "queued": True}
        return {"accepted": True}

    def _m_session_cancel(self, _params: dict) -> dict:
        """Esc: the turn that runs stops; a line just taken (or sent and not
        yet read) does not start its turn; no automatic follow-up runs.
        Lines sent after the cancel are not touched by it, and what it
        arms ends when the loop is back at the prompt."""
        agent = self.agent
        refused: list = []
        with agent.cond:
            in_turn = agent.busy
            starting = not agent.idle and not agent.busy
            unread = agent.idle and bool(agent.inbox)
            stopping = in_turn or starting or unread
            if stopping:
                agent.cancel_requested = True
                agent.followups_off = True
                # The line the cancel is about: the one in progress, or the
                # first one waiting at the prompt (read next).
                agent.cancel_gen = agent.line_gen + (1 if unread else 0)
            if in_turn or self.broker.pending():
                # Every open question is refused, and so is every one the
                # ending turn still asks (lifted at the prompt). Under the
                # lock: a turn starting now must not lift it (see
                # ProtocolAgent._begin_turn_resets_the_abort).
                refused = self.broker.abort_all()
        if in_turn:
            agent._stop_engine()
        agent.wake.set()
        return {"stopping": stopping, "refused": len(refused)}

    def _m_session_setMode(self, params: dict) -> dict:
        mode = str(params.get("mode") or "")
        if mode not in ("plan", "default", "acceptEdits"):
            raise _Refusal(protocol.INVALID_PARAMS,
                           "mode must be plan, default or acceptEdits "
                           "(an unattended session is started as one, never "
                           "switched into)")
        self.engine.set_kit_permission_mode(mode)
        self.agent.post(self.check_session)
        return {"approval": self._session_facts()["approval"]}

    def _m_session_status(self, _params: dict) -> dict:
        try:
            status = dict(self.engine.get_status() or {})
        except Exception:
            status = {}
        status.update({"state": self.agent.state,
                       "approval": self._session_facts()["approval"],
                       "queued": len(self.agent.queued)})
        return protocol_safe(status)

    def _m_session_replay(self, params: dict) -> dict:
        turns = int(params.get("turns", 0) or 0)
        return {"events": self.log.replay(turns, header=self.header(),
                                          max_bytes=REPLAY_BUDGET_BYTES)}

    def _m_permission_reply(self, params: dict) -> dict:
        from .approvals import InvalidAnswer, UnknownRequest

        try:
            landed, message = self.broker.answer(
                str(params.get("requestId") or ""),
                choice=str(params.get("choice") or ""),
                answers=list(params.get("answers") or []),
                reason=str(params.get("reason") or ""),
                persist=bool(params.get("persist", False)))
        except UnknownRequest as exc:
            raise _Refusal(protocol.UNKNOWN_REQUEST, str(exc)) from None
        except InvalidAnswer as exc:
            raise _Refusal(protocol.INVALID_PARAMS, str(exc)) from None
        self.agent.post(self.check_session)
        return {"applied": landed, "message": message}

    def _m_session_commands(self, _params: dict) -> dict:
        return {"commands": choices.command_list(self.workspace)}

    def _m_session_completePath(self, params: dict) -> dict:
        return {"paths": choices.complete_path(str(params.get("fragment") or ""),
                                               self.workspace)}

    def _m_session_toggleTasks(self, _params: dict) -> dict:
        """ctrl+t, the REPL's way: the end-of-turn task list on or off."""
        agent = self.agent
        with agent.cond:
            agent._show_tasks = not agent._show_tasks
            on = agent._show_tasks
        agent.post(lambda: agent.transcript.chrome(
            agent.transcript.theme.dim(f"task list {'on' if on else 'off'}")))
        return {"showTasks": on}

    def _m_session_choices(self, params: dict) -> dict:
        try:
            return choices.choices(self.engine, str(params.get("what") or ""))
        except ValueError as exc:
            raise _Refusal(protocol.INVALID_PARAMS, str(exc)) from None

    # -- lifecycle -------------------------------------------------------------------
    def shutdown(self, reason: str, code: int = 0) -> None:
        """End the session: no next turn, every question refused (now and
        from any thread still to ask -- for good, see ProtocolBroker.close),
        the loop out by the way ``reason`` calls for (``_leave``)."""
        agent = self.agent
        with agent.cond:
            if agent.closing.is_set():
                return
            agent.leave_reason = reason
            agent.leave_code = code
            agent.closing.set()
            agent.inbox.clear()
            # A stop that must hold, like a cancel's: the engine clears a
            # stop request when its turn begins (_look_around repeats it).
            agent.cancel_requested = True
            agent.cond.notify_all()
        agent.wake.set()
        try:
            self.broker.close()
        except Exception:
            pass
        agent._stop_engine()
        if reason != "shutdown":
            # Nobody is left to wait for: a tool call or a !command must not
            # keep the session -- and its saved state -- waiting.
            agent.leave_now()

    def _read_requests(self) -> None:
        def bad(problem: str) -> None:
            self.writer.send(protocol.error(None, protocol.PARSE_ERROR, problem))

        try:
            for message in protocol.read_messages(self.reader, on_bad_line=bad):
                with self._handling:
                    self.handle(message)
        except Exception as exc:                    # noqa: BLE001
            print(f"worker: request reader failed: {exc}", file=sys.stderr)
        finally:
            # The window side went away (EOF): a hangup, for the session.
            self.shutdown("eof", 129)

    def _look_around(self) -> None:
        """Every tick: posture or model changes; during a turn also usage,
        presence, operator mail and a cancel that must hold."""
        last = None
        agent = self.agent
        while not agent.ended:
            time.sleep(_TICK_S)
            if agent.cancel_requested and agent.busy:
                # The engine clears a stop request when its turn begins; a
                # cancel (or shutdown) that came early is repeated until
                # the turn is over.
                agent._stop_engine()
            if agent.closing.is_set():
                continue
            agent.post(self.check_session)
            if not agent.busy:
                continue
            # The REPL's pump does these only on a terminal it reads keys
            # from; the worker has none, so they are done from here. Mail
            # is posted: it renders into the stream it is steered into.
            try:
                agent._announce_presence()
            except Exception:
                pass
            agent.post(self._operator_mail)
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

    def _operator_mail(self) -> None:
        try:
            self.agent._steer_operator_mail()
        except Exception:
            pass

    def _end_reason(self, code: int) -> str:
        agent = self.agent
        if agent.leave_reason:
            return agent.leave_reason
        if getattr(agent, "_quit", False):
            return "exit"
        if code >= 128:
            return "signal"
        return "closed"

    def run(self) -> int:
        self.emitter.emit("worker.ready", protocolVersion=protocol.PROTOCOL_VERSION,
                          pid=os.getpid())
        self.emitter.emit("session.info", **self.info())
        self.agent.refresh_state()
        threading.Thread(target=self._read_requests, name="protocol-reader",
                         daemon=True).start()
        threading.Thread(target=self._look_around, name="look-around",
                         daemon=True).start()
        import signal

        def outside_the_loop(signum, _frame) -> None:
            # Before the loop installs its handlers and after it restores
            # them: end the loop when it runs, but never kill the process
            # mid-save. No lock: the main thread may hold one right now;
            # the loop sees the flag at its next look (_IDLE_WAIT_S).
            agent = self.agent
            if not agent.closing.is_set():
                agent.leave_reason = agent.leave_reason or "signal"
                agent.leave_code = agent.leave_code or 128 + int(signum)
                agent.closing.set()

        for sig in (signal.SIGHUP, signal.SIGTERM):
            try:
                signal.signal(sig, outside_the_loop)
            except (ValueError, OSError):
                pass
        code = 0
        try:
            code = int(self.agent.run() or 0)
        finally:
            with self.agent.cond:
                self.agent.closing.set()
                self.agent.cond.notify_all()
            if self.save:
                end_real_session(self.engine, self.workspace, self.name,
                                 learn=not self.fake)
            # A request being handled right now gets its answer out first.
            if self._handling.acquire(timeout=2.0):
                self._handling.release()
            self.agent.ended = True
            self.agent.refresh_state()
            self.emitter.emit("session.ended", reason=self._end_reason(code),
                              exitCode=code)
            self.log.close()
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
    p.add_argument("--fake-save", action="store_true",
                   help="save a fake session to the session store at the end "
                        "(tests of the save path)")
    p.add_argument("--fake-pause", type=float, default=0.0,
                   help="wait this long during a (fake) outage instead of the "
                        "REPL's 60 s and more (tests of the wait)")
    return p.parse_args(argv)


def guard_process(kind: str) -> list[str]:
    """What ``cli.main`` does for every agent process before its session.

    Refuse to start when the principles were tampered with; end the
    background shells this session starts when it ends; make the process
    unreadable by the commands it runs (it holds the model's key) and
    register it so ``delfin-agent stop-all`` finds it; and end it when an
    emergency stop is given. Returns the notes ``cli.main`` would print --
    a security update, a provider key exported in the environment -- for
    the window to show.
    """
    import atexit

    from .. import cli as _cli
    from .. import lifeline as _lifeline
    from .. import process_guard as _process_guard

    _cli.guard_principles_or_exit()
    atexit.register(_cli._stop_own_background_shells)
    _process_guard.protect(kind)
    said = io.StringIO()
    with contextlib.redirect_stderr(said):
        try:
            _cli._show_security_notices()
        except Exception:
            pass
    notes = _lines(said.getvalue())
    try:
        exported = _process_guard.exported_provider_keys()
        if exported:
            notes.append("Warning: " + _process_guard.exported_key_advice(exported))
    except Exception:
        pass

    def _end_with_everything() -> None:
        _cli._stop_own_background_shells()
        _lifeline.exit_now()

    _lifeline.watch(_end_with_everything)
    return notes


def main(argv: list[str] | None = None) -> int:
    from .. import repl as _repl
    from .. import state_paths

    # Before anything imports a state sink, as cli.main does.
    state_paths.scratch_state_from_environment()
    args = parse_args(argv)
    reader, raw_out = claim_stdio()
    if os.environ.get("DELFIN_TUI_STACKS", "").strip() not in ("", "0"):
        # kill -USR1 <pid>: every thread's stack into the log, for a hang.
        import faulthandler
        import signal
        faulthandler.register(signal.SIGUSR1, all_threads=True)
    writer = protocol.LineWriter(raw_out)
    workspace = Path(args.cwd or os.getcwd()).expanduser().resolve()
    # This process runs one session and has no terminal, so two of the
    # REPL's reading conventions are set for it here:
    # * a line that is just three quotes opens a block, and a trailing
    #   backslash means "more follows" -- ways to type several lines into
    #   a one-line prompt. A window sends each message whole; with the
    #   grammar on, a message ending in a backslash would leave the session
    #   waiting for input nobody sends.
    # * the idle read sleeps a quarter second between looks at the
    #   keyboard. The worker's look (``ProtocolAgent._line_ready``) waits
    #   for a message itself, so the sleep would only delay every prompt.
    _repl._read_block_tail = _whole_message
    _wakeable = _repl.read_block_wakeable

    def _read_without_sleep(read_line, wake, *, ready, prompt="> ",
                            cont="… ", tick=0.25):
        return _wakeable(read_line, wake, ready=ready, prompt=prompt,
                         cont=cont, tick=0.0)

    _repl.read_block_wakeable = _read_without_sleep
    if args.fake and args.fake_pause > 0:
        _repl.TerminalAgent._ENDPOINT_PAUSE_S = (args.fake_pause,) * 3
        _repl.TerminalAgent._ENDPOINT_PAUSE_CEIL_S = args.fake_pause

    def refuse(reason: str, code: int) -> int:
        writer.send(protocol.event({"type": "session.ended", "session": "",
                                    "ts": round(time.time(), 3),
                                    "reason": reason, "exitCode": code}))
        return code

    try:
        notes = guard_process("tui session")
    except SystemExit as exc:
        return refuse(str(exc.code or exc), 2)
    try:
        if args.fake:
            from .fake_engine import FakeEngine
            engine = FakeEngine(str(workspace), delay=args.fake_delay)
            if args.permission_mode:
                engine.set_kit_permission_mode(args.permission_mode,
                                               announce=False)
            banner = ["delfin-agent · fake/fake-model · solo",
                      f"workspace  {workspace}",
                      f"approval   {engine.kit_permissions.mode}"]
            notices: list[str] = []
        else:
            engine, _ns, banner_text, notices = open_real_session(args, workspace)
            banner = banner_text.splitlines()
    except SessionRefused as exc:
        return refuse(str(exc), 2)
    except Exception as exc:                        # noqa: BLE001
        return refuse(f"The session could not be opened: "
                      f"{type(exc).__name__}: {exc}", 3)
    worker = Worker(engine, writer=writer, reader=reader, workspace=workspace,
                    banner=banner, notices=notes + notices, name=args.name,
                    fake=args.fake, save=(not args.fake) or args.fake_save)
    return worker.run()


if __name__ == "__main__":
    sys.exit(main())
