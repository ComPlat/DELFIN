"""The terminal controller: one turn, rendered as it happens.

Division of labour, and the reason for it: ``repl_render`` decides what a
line says, this module decides when it is written and to which stream.
Keeping the two apart is what stops this file growing into the thing it
replaces — a single closure where formatting, threading and state are the
same 14 000 lines.

Stream discipline: **stdout is the answer, stderr is everything else.**
Tool lines, notices, thinking, the banner and the prompt all go to stderr,
so ``delfin-agent -p "..." > answer.txt`` produces exactly the answer and
nothing about how it was reached.
"""

from __future__ import annotations

import queue
import signal
import sys
import threading
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

from . import repl_render as rr

__all__ = [
    "RenderItem", "Transcript", "TurnResult", "run_turn", "TURN_KEYS",
    "ReplOptions", "TerminalAgent", "read_block", "HISTORY_NAME",
]

HISTORY_NAME = "agent_repl_history"
_HISTORY_LINES = 1000
_PUMP_TICK_S = 0.05
_JOIN_NOTICE_AFTER_S = 2.0


# The keys `cli._run_once` returns, and therefore what --output-format json
# promises. Named here so the two paths can be asserted equal instead of
# drifting into two shapes of the same answer.
TURN_KEYS: frozenset[str] = frozenset({
    "text", "tool_calls", "input_tokens", "output_tokens", "error",
})


@dataclass(frozen=True)
class RenderItem:
    """One thing that happened, on its way from the worker to the screen."""

    kind: str                      # text|thinking|tool_use|tool_result
                                   # |notice|denied|error|done
    text: str = ""
    name: str = ""
    data: dict | None = None


class Transcript:
    """Owns both streams and the cursor column.

    The column tracking is not decoration. The model streams
    ``"Let me check the tests"`` with no trailing newline and then calls a
    tool; without closing the line first, the tool line is appended to the
    middle of the sentence and the transcript stops being readable exactly
    when it matters most.
    """

    def __init__(self, out=None, err=None, *, theme: rr.Theme | None = None,
                 show_tools: bool = True, show_thinking: bool = False,
                 width: int | None = None) -> None:
        self.out = out if out is not None else sys.stdout
        self.err = err if err is not None else sys.stderr
        self.theme = theme if theme is not None else rr.theme_for(self.err)
        self.show_tools = show_tools
        self.show_thinking = show_thinking
        self.width = width or rr.terminal_width()
        # Two flags, because one had to answer two different questions and
        # got one of them wrong. `_out_open` is about the FILE: stdout ends
        # mid-line and needs closing before the process exits.
        # `_break_pending` is about the SCREEN: the next chrome line needs a
        # visual break first, and only the first one does. Sharing a flag
        # meant the first tool line cleared it and stdout was left without
        # its final newline.
        self._out_open = False
        self._break_pending = False

    # -- primitives ------------------------------------------------------
    def answer(self, delta: str) -> None:
        """Answer text. The only thing that ever reaches stdout."""
        if not delta:
            return
        self.out.write(delta)
        self._flush(self.out)
        self._out_open = not delta.endswith("\n")
        self._break_pending = self._out_open

    def chrome(self, line: str) -> None:
        """Everything that is not the answer, on stderr, on its own line.

        The break that closes a half-written answer line is written to
        STDERR, not stdout. Both streams share one cursor when they share a
        terminal, so the visual effect is identical — but a newline on
        stdout would land inside the answer, and `delfin-agent -p "..." >
        answer.txt` would capture a line break that the model never wrote.
        Redirecting stderr alone costs the visual break and keeps the
        answer exact, which is the right way round.
        """
        if not line:
            return
        if self._break_pending:
            self.err.write("\n")
            self._break_pending = False
        self.err.write(line + "\n")
        self._flush(self.err)

    def finish(self) -> None:
        """Close the answer stream, so a redirected stdout ends in a newline."""
        if self._out_open:
            self.out.write("\n")
            self._flush(self.out)
            self._out_open = False
            self._break_pending = False

    @staticmethod
    def _flush(stream) -> None:
        try:
            stream.flush()
        except Exception:
            pass

    # -- dispatch --------------------------------------------------------
    def render(self, item: RenderItem) -> None:
        kind = item.kind
        if kind == "text":
            self.answer(item.text)
        elif kind == "notice":
            self.chrome(rr.notice_line(item.text, theme=self.theme))
        elif kind == "thinking":
            if self.show_thinking:
                self.chrome(rr.thinking_line(
                    item.text, width=self.width, theme=self.theme))
        elif kind == "tool_use":
            if self.show_tools:
                self.chrome(rr.tool_headline(
                    item.name, item.text, width=self.width, theme=self.theme))
        elif kind == "tool_result":
            if self.show_tools:
                self.chrome(rr.tool_result_line(
                    item.name, item.text, meta=item.data,
                    width=self.width, theme=self.theme))
        elif kind == "denied":
            self.chrome(rr.denied_line(item.name, theme=self.theme))
        elif kind == "error":
            self.chrome(self.theme.red(f"[error] {item.text}"))


@dataclass
class TurnResult:
    text: str = ""
    tool_calls: list[dict] = field(default_factory=list)
    input_tokens: int = 0
    output_tokens: int = 0
    error: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "text": self.text,
            "tool_calls": self.tool_calls,
            "input_tokens": self.input_tokens,
            "output_tokens": self.output_tokens,
            "error": self.error,
        }


def _usage(engine) -> tuple[int, int]:
    usage = getattr(engine, "token_usage", {}) or {}
    return int(usage.get("input", 0) or 0), int(usage.get("output", 0) or 0)


def run_turn(engine, prompt: str, *, sink: Callable[[RenderItem], None],
             max_tokens: int = 0, memory_context: str = "") -> TurnResult:
    """One turn, with every callback turned into a RenderItem.

    The accounting mirrors ``cli._run_once`` deliberately — same token
    deltas, same result keys — so the interactive and the headless path
    report the same turn the same way.

    ``sink`` is called from whatever thread ``stream_response`` runs on.
    It must not touch the terminal; in the REPL it is a queue put, and the
    main thread does the writing.
    """
    chunks: list[str] = []
    tool_calls: list[dict] = []
    result = TurnResult()

    def _emit(item: RenderItem) -> None:
        try:
            sink(item)
        except Exception:
            pass

    def _on_token(text: str) -> None:
        if text:
            chunks.append(text)
            _emit(RenderItem("text", text=text))

    def _on_notice(text: str) -> None:
        if text:
            _emit(RenderItem("notice", text=text))

    def _on_thinking(text: str) -> None:
        if text:
            _emit(RenderItem("thinking", text=text))

    def _on_tool_use(name: str, input_json: str) -> None:
        tool_calls.append({"name": name, "input": rr._as_dict(input_json)})
        _emit(RenderItem("tool_use", name=name, text=input_json or ""))

    # One result, one line. The engine reports a tool result on two
    # callbacks — the head slice on on_tool_result, the verdict and the
    # true size on on_tool_result_meta — and rendering each as it arrives
    # drew every call twice, the second time with a body it did not have.
    # The text is held until the verdict catches up.
    pending: dict[str, Any] = {"name": "", "text": "", "held": False}

    def _flush_pending() -> None:
        if pending["held"]:
            _emit(RenderItem("tool_result", name=pending["name"],
                             text=pending["text"]))
            pending.update(name="", text="", held=False)

    def _on_tool_result(name: str, output: str) -> None:
        # A second result before the first was completed means the meta
        # callback is not firing on this path; do not lose the first line.
        _flush_pending()
        pending.update(name=name, text=output or "", held=True)

    def _on_tool_result_meta(name: str, meta: dict) -> None:
        text = pending["text"] if pending["held"] and pending["name"] == name else ""
        if pending["held"] and pending["name"] != name:
            _flush_pending()
        pending.update(name="", text="", held=False)
        _emit(RenderItem("tool_result", name=name, text=text, data=dict(meta)))

    def _on_denied(name: str) -> None:
        _emit(RenderItem("denied", name=name))

    in_before, out_before = _usage(engine)
    try:
        full_text = engine.stream_response(
            user_message=prompt,
            memory_context=memory_context,
            on_token=_on_token,
            on_notice=_on_notice,
            on_thinking=_on_thinking,
            on_tool_use=_on_tool_use,
            on_tool_result=_on_tool_result,
            on_tool_result_meta=_on_tool_result_meta,
            on_permission_denied=_on_denied,
            max_tokens=max_tokens,
        ) or ""
    except Exception as exc:
        result.error = str(exc)
        full_text = ""
        _emit(RenderItem("error", text=result.error))

    _flush_pending()
    in_after, out_after = _usage(engine)
    result.text = (full_text or "".join(chunks)).strip()
    result.tool_calls = tool_calls
    result.input_tokens = max(0, in_after - in_before)
    result.output_tokens = max(0, out_after - out_before)

    # The answer was streamed before the sanitiser ran on it, and a
    # terminal cannot un-print. Say so rather than let the screen and the
    # recorded answer quietly disagree.
    try:
        from .text_sanitize import sanitize_agent_text
        cleaned = sanitize_agent_text(result.text).text
        if cleaned != result.text:
            result.text = cleaned
            _emit(RenderItem("notice", text=(
                "The model leaked tool-channel text into the answer above; "
                "the recorded answer is the cleaned version.")))
    except Exception:
        pass

    _emit(RenderItem("done"))
    return result


# ---------------------------------------------------------------------------
# The interactive loop
# ---------------------------------------------------------------------------

_BLOCK_FENCE = '"""'


def read_block(read_line: Callable[[str], str], *, prompt: str = "> ",
               cont: str = "… ") -> str:
    """One message, which may span several lines.

    Two forms, both deterministic and both testable without a terminal: a
    line ending in a backslash continues onto the next, and a line that is
    exactly \"\"\" opens a block that ends at the next one.

    Bracketed paste is the real fix for a multi-line paste becoming N
    separate turns, and it needs raw mode — so it arrives with the key
    layer, and the fence is what the help text points at until then.
    """
    first = read_line(prompt)
    if first.strip() == _BLOCK_FENCE:
        lines: list[str] = []
        while True:
            line = read_line(cont)
            if line.strip() == _BLOCK_FENCE:
                break
            lines.append(line)
        return "\n".join(lines)

    parts = [first]
    while parts[-1].endswith("\\"):
        parts[-1] = parts[-1][:-1]
        parts.append(read_line(cont))
    return "".join(parts) if len(parts) > 1 else parts[0]


@dataclass
class ReplOptions:
    cwd: Path = field(default_factory=Path.cwd)
    max_tokens: int = 0
    show_thinking: bool = False
    show_tools: bool = True
    color: str = "auto"
    banner: str = ""


class TerminalAgent:
    """Reads, runs one turn, renders it, repeats.

    The threading model is the whole design, and it is deliberately
    lopsided:

      MAIN THREAD                        WORKER ("agent-turn", daemon)
      owns stdin, owns every write       runs engine.stream_response
      pumps the queue and renders        callbacks only put on the queue
      handles SIGINT                     never touches the terminal

    Two consequences that are the point rather than a side effect.

    A ``KeyboardInterrupt`` cannot cross ``stream_response``, because
    Python delivers signals to the main thread and the turn runs on the
    worker. That matters more than it looks: the engine's turn body
    catches ``Exception``, not ``BaseException``, so an interrupt raised
    inside it would skip the cleanup that pops the unanswered user
    message — and the alternation sanitiser then resolves two consecutive
    user messages by keeping the NEWEST, silently overwriting the question
    that was asked. An interrupt has to arrive as ``request_stop()``, and
    the turn has to unwind through its own stop path.

    And rendering happens on one thread only, so a half-written answer
    line cannot be cut in two by a tool line from somewhere else.
    """

    def __init__(self, engine, opts: ReplOptions | None = None, *,
                 out=None, err=None, read_line: Callable[[str], str] | None = None):
        self.engine = engine
        self.opts = opts or ReplOptions()
        self.out = out if out is not None else sys.stdout
        self.err = err if err is not None else sys.stderr
        self.transcript = Transcript(
            self.out, self.err,
            theme=rr.theme_for(self.err, self.opts.color),
            show_tools=self.opts.show_tools,
            show_thinking=self.opts.show_thinking,
        )
        self._read_line = read_line or self._input
        self._q: queue.Queue[RenderItem] = queue.Queue()
        self._turn_active = threading.Event()
        self._interrupts = 0
        self._idle_interrupts = 0
        self._prev_sigint = None
        self._stdin = sys.stdin
        # What the user is typing WHILE a turn runs, and the last tool
        # result, so Ctrl+O has something to expand.
        self._input_line = ""
        self._last_result_text = ""
        # Typed during a turn: queued, never injected. A queued message
        # cannot be lost and cannot land in a context nobody could see.
        self.queued: list[str] = []

    # -- input -----------------------------------------------------------
    def _input(self, prompt: str) -> str:
        """Read one line, without the prompt landing in the answer.

        ``input(prompt)`` writes the prompt to STDOUT, which is the one
        stream that must carry the answer and nothing else — so
        `delfin-agent > answer.txt` would collect a `> ` for every turn.

        When stdout is a terminal the prompt still goes through ``input``,
        because readline needs to know the prompt to redraw the line
        correctly during history search and on a wrapped line. When stdout
        is redirected there is nothing to redraw, so the prompt goes to
        stderr where the rest of the chrome lives.
        """
        try:
            interactive_out = bool(self.out.isatty())
        except Exception:
            interactive_out = False
        if interactive_out:
            return input(prompt)
        self.err.write(prompt)
        self._flush_err()
        return input()

    def _flush_err(self) -> None:
        try:
            self.err.flush()
        except Exception:
            pass

    def _notice(self, text: str) -> None:
        self.transcript.render(RenderItem("notice", text=text))

    # -- interrupts ------------------------------------------------------
    def _on_sigint(self, signum, frame) -> None:
        """Runs on the main thread, BETWEEN two bytecodes — possibly inside
        a write. So it never prints: it puts on the queue and lets the pump
        do the writing. Getting that wrong produces interleaved output once
        every few interrupts, which is the worst kind of bug to reproduce.
        """
        if not self._turn_active.is_set():
            raise KeyboardInterrupt
        self._interrupts += 1
        if self._interrupts == 1:
            self._stop_engine()
            self._q.put(RenderItem("notice", text=(
                "Interrupt — ending this turn. A tool call already running "
                "finishes first; a shell command is capped at its own "
                "timeout, so this can take a moment.")))
            return
        if self._interrupts == 2:
            self._q.put(RenderItem("notice", text=(
                "Still waiting on the running tool call. Once more abandons "
                "the process, and anything it started keeps running.")))
            return
        raise KeyboardInterrupt

    def _stop_engine(self) -> None:
        """The dashboard's stop pair, so both surfaces stop the same way."""
        try:
            self.engine.request_stop()
        except Exception:
            pass
        client = getattr(self.engine, "client", None)
        if client is not None and hasattr(client, "signal_stop"):
            try:
                client.signal_stop()
            except Exception:
                pass

    # -- one turn --------------------------------------------------------
    def turn(self, prompt: str) -> TurnResult:
        result_box: list[TurnResult] = []

        def _worker() -> None:
            try:
                result_box.append(run_turn(
                    self.engine, prompt, sink=self._q.put,
                    max_tokens=self.opts.max_tokens))
            except BaseException as exc:            # noqa: BLE001
                # The worker must never die silently: a turn that vanished
                # looks exactly like a turn that answered nothing.
                self._q.put(RenderItem("error", text=str(exc)))
                self._q.put(RenderItem("done"))

        self._interrupts = 0
        self._turn_active.set()
        worker = threading.Thread(target=_worker, name="agent-turn", daemon=True)
        worker.start()
        try:
            self._pump(worker)
        finally:
            self._turn_active.clear()
            self._settle(worker)
        self.transcript.finish()
        return result_box[0] if result_box else TurnResult(error="turn produced nothing")

    def _pump(self, worker: threading.Thread) -> None:
        from . import repl_keys as rk

        with rk.RawMode(self._stdin) as raw:
            decoder = rk.KeyDecoder()
            while True:
                if raw.active:
                    for event in decoder.feed(raw.read_ready(_PUMP_TICK_S)):
                        self._on_key(event, decoder)
                try:
                    item = self._q.get(
                        timeout=0.0 if raw.active else _PUMP_TICK_S)
                except queue.Empty:
                    if not worker.is_alive() and self._q.empty():
                        self._clear_input_line()
                        return
                    continue
                if item.kind == "done":
                    self._clear_input_line()
                    return
                self._render_around_input(item)

    # -- keys during a turn ----------------------------------------------
    def _on_key(self, event, decoder) -> None:
        from . import repl_keys as rk

        if event.kind == rk.INTERRUPT:
            self._stop_engine()
            self._clear_input_line()
            self.transcript.chrome(self.transcript.theme.yellow(
                "! Esc — ending this turn. A running tool call finishes "
                "first."))
            return
        if event.kind == rk.SUBMIT:
            text = event.text.strip()
            self._clear_input_line()
            if text:
                self.queued.append(text)
                self.transcript.chrome(self.transcript.theme.dim(
                    f"queued ({len(self.queued)}) — goes out when this turn "
                    "ends"))
            return
        if event.kind == rk.CYCLE_MODE:
            self._clear_input_line()
            self._cycle_mode()
            self._draw_input_line(decoder.buffer)
            return
        if event.kind == rk.EXPAND:
            self._clear_input_line()
            self._expand_last_result()
            self._draw_input_line(decoder.buffer)
            return
        if event.kind == rk.REDRAW:
            self._clear_input_line()
            self._draw_input_line(decoder.buffer)
            return
        if event.kind == rk.EDIT:
            self._draw_input_line(event.text)

    def _cycle_mode(self) -> None:
        """Shift+Tab, one step along the ladder — never onto bypass.

        Reaching unattended execution must stay a thing someone types on
        purpose, not a thing a key lands on while a turn is running.
        """
        order = ["plan", "default", "acceptEdits"]
        try:
            current = str(getattr(self.engine.kit_permissions, "mode", "")
                          or "default")
        except Exception:
            self.transcript.chrome(self.transcript.theme.dim(
                "this backend carries no permission gate"))
            return
        try:
            nxt = order[(order.index(current) + 1) % len(order)]
        except ValueError:
            nxt = "default"
        try:
            self.engine.set_kit_permission_mode(nxt)
        except Exception as exc:
            self.transcript.chrome(self.transcript.theme.red(
                f"could not switch mode: {exc}"))
            return
        self.transcript.chrome(self.transcript.theme.cyan(
            f"approval → {nxt}"))

    def _expand_last_result(self) -> None:
        if not self._last_result_text:
            self.transcript.chrome(self.transcript.theme.dim(
                "nothing to expand yet"))
            return
        for line in rr.strip_control(self._last_result_text).splitlines():
            self.transcript.chrome("  " + line)

    # -- the line being typed under a running turn ------------------------
    def _can_redraw(self) -> bool:
        """Only a terminal gets cursor control.

        Redrawing needs \\r and an erase sequence, and those are exactly
        the characters this codebase strips out of tool output before
        printing it. Writing them into a redirected stderr would put
        control codes in a log file for a line nobody can see anyway.
        """
        try:
            return bool(self.err.isatty())
        except Exception:
            return False

    def _draw_input_line(self, text: str) -> None:
        self._input_line = text or ""
        if not text or not self._can_redraw():
            return
        self.err.write("\r\x1b[K» " + text)
        self._flush_err()

    def _clear_input_line(self) -> None:
        if not self._input_line:
            return
        self._input_line = ""
        if not self._can_redraw():
            return
        self.err.write("\r\x1b[K")
        self._flush_err()

    def _render_around_input(self, item: RenderItem) -> None:
        """Rendering must not tear the line the user is typing.

        Erase it, write the transcript line, put it back — the classic
        bottom-line problem, and the reason typing during a turn is worth
        having rather than merely possible.
        """
        held = self._input_line
        if held:
            self._clear_input_line()
        if item.kind == "tool_result" and item.text:
            self._last_result_text = item.text
        self.transcript.render(item)
        if held:
            self._draw_input_line(held)

    def _settle(self, worker: threading.Thread) -> None:
        """Join, THEN clear the stop. Both halves are load-bearing.

        The turn gate refuses a second concurrent turn by RETURNING a
        sentence rather than raising, so returning to the prompt while a
        stopped worker is still unwinding makes the next turn render
        machinery speech as the model's answer. And ``clear_stop`` refuses
        while the owning turn is in flight, so clearing before the join is
        a silent no-op that leaves the brake armed for the next turn.
        """
        worker.join(timeout=_JOIN_NOTICE_AFTER_S)
        if worker.is_alive():
            self._notice("Waiting for the running tool call to finish…")
            worker.join()
        # Drain whatever the worker queued while we were joining.
        while True:
            try:
                item = self._q.get_nowait()
            except queue.Empty:
                break
            if item.kind != "done":
                self.transcript.render(item)
        try:
            self.engine.clear_stop()
        except Exception:
            pass

    # -- the loop --------------------------------------------------------
    def run(self, first_prompt: str = "") -> int:
        self._install_sigint()
        self._load_history()
        if self.opts.banner:
            for line in self.opts.banner.splitlines():
                self.transcript.chrome(line)
        pending = first_prompt.strip()
        try:
            while True:
                if not pending and self.queued:
                    pending = self.queued.pop(0)
                if not pending:
                    try:
                        pending = read_block(self._read_line).strip()
                    except EOFError:
                        self.transcript.chrome("")
                        return 0
                    except KeyboardInterrupt:
                        self._idle_interrupts += 1
                        if self._idle_interrupts >= 2:
                            self.transcript.chrome("")
                            return 130
                        self.transcript.chrome(
                            "(interrupt — once more to leave)")
                        continue
                self._idle_interrupts = 0
                if not pending:
                    continue
                if pending in ("/exit", "/quit"):
                    return 0
                try:
                    self.turn(pending)
                except KeyboardInterrupt:
                    self.transcript.chrome("")
                    return 130
                pending = ""
        finally:
            self._save_history()
            self._restore_sigint()

    # -- plumbing --------------------------------------------------------
    def _install_sigint(self) -> None:
        try:
            self._prev_sigint = signal.signal(signal.SIGINT, self._on_sigint)
        except (ValueError, OSError):
            # Not the main thread, or no signal support. The loop still
            # works; Ctrl+C simply behaves as the default.
            self._prev_sigint = None

    def _restore_sigint(self) -> None:
        if self._prev_sigint is not None:
            try:
                signal.signal(signal.SIGINT, self._prev_sigint)
            except (ValueError, OSError):
                pass

    @staticmethod
    def history_path() -> Path:
        from . import state_paths
        return state_paths.ensure_dir(Path.home() / ".delfin") / HISTORY_NAME

    def _load_history(self) -> None:
        try:
            import readline
        except ImportError:
            return
        try:
            path = self.history_path()
            readline.set_history_length(_HISTORY_LINES)
            if path.exists():
                readline.read_history_file(str(path))
        except Exception:
            # libedit behaves differently, and a missing file is normal on
            # a first run. Neither is worth refusing to start over.
            pass

    def _save_history(self) -> None:
        try:
            import readline
            from . import state_paths
            path = self.history_path()
            readline.write_history_file(str(path))
            # Prompts describe the user's project, so the file is theirs
            # alone -- the same 0600 every other state file here gets.
            state_paths.secure_file(path)
        except Exception:
            pass
