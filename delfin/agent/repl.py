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

import os
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
    "permission_mode",
]

HISTORY_NAME = "agent_repl_history"
_HISTORY_LINES = 1000
_PUMP_TICK_S = 0.05
_JOIN_NOTICE_AFTER_S = 2.0
# What leaving costs at most. The last step of the interrupt ladder unwinds
# through the settle, so the join there has to be bounded: an unbounded one
# waits for the tool call the interrupt exists to walk away from.
_JOIN_ABANDON_S = 1.0
# How long a session whose terminal is gone may take to save and clean up.
# That work runs on storage that may be the very thing that went wrong, and
# waiting on it is how a session outlives its terminal.
_LEAVE_DEADLINE_S = 10.0
_STATUS_TICK_S = 0.25
_SPINNER = "⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"


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
    # Kept on a result because apply_patch returns only a success object;
    # the actual diff lives in the matching call's input.  Carrying it
    # here lets the renderer show what changed without reading the file a
    # second time or guessing from the current worktree.
    tool_input: dict | None = None


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
                 width: int | None = None, color: str = "auto") -> None:
        self.out = out if out is not None else sys.stdout
        self.err = err if err is not None else sys.stderr
        self.theme = theme if theme is not None else rr.theme_for(self.err)
        self.show_tools = show_tools
        self.show_thinking = show_thinking
        # A fixed width or a live one. Read once, a terminal resized after
        # launch went on truncating to the size it had at startup — every
        # tool line and every status row cut at the wrong column for the
        # rest of the session. The fixed form stays available because the
        # tests need a stable width to assert against.
        self._fixed_width = int(width) if width else 0
        self._width_cache = self._fixed_width or rr.terminal_width()
        # A SECOND colour decision, taken about stdout. The theme above
        # is about stderr, where the chrome goes; the answer goes to
        # stdout, and the two streams are redirected independently. Asking
        # stderr about stdout is how `delfin-agent -p "..." > answer.txt`
        # in a terminal would end up writing escape codes into the file.
        self._answer_theme = rr.theme_for(self.out, color)
        self._markdown = rr.MarkdownStream(self._answer_theme)
        # stdout and stderr normally name the same terminal.  Remember
        # that fact explicitly: while an answer is streaming the turn UI
        # is painted on stderr, then the next answer delta continues on
        # stdout.  A redirected stdout is a different surface and must
        # never influence where the stderr cursor is restored.
        self._shared_terminal = self._streams_share_terminal(
            self.out, self.err)
        self._screen_answer_open = False
        self._screen_answer_column = 0
        # Two flags, because one had to answer two different questions and
        # got one of them wrong. `_out_open` is about the FILE: stdout ends
        # mid-line and needs closing before the process exits.
        # `_break_pending` is about the SCREEN: the next chrome line needs a
        # visual break first, and only the first one does. Sharing a flag
        # meant the first tool line cleared it and stdout was left without
        # its final newline.
        self._out_open = False
        self._break_pending = False

    @property
    def width(self) -> int:
        return self._fixed_width or self._width_cache

    @width.setter
    def width(self, value: int) -> None:
        """Assigning a width PINS it, which is what assigning one meant
        before this became a property. Callers that set it are saying
        "this is the width", not "start following the terminal from
        here", and a resize must not undo that."""
        self._fixed_width = int(value or 0)

    def refresh_width(self) -> int:
        """Re-ask the terminal. Called from the resize handler."""
        if not self._fixed_width:
            self._width_cache = rr.terminal_width()
        return self.width

    # -- primitives ------------------------------------------------------
    def answer(self, delta: str) -> None:
        """Answer text. The only thing that ever reaches stdout.

        The line-open flags follow the MODEL's text, never the styled
        form: an SGR reset is not a newline, and a delta that ends in one
        would otherwise be read as a closed line and let the status row
        repaint over the sentence being written.
        """
        if not delta:
            return
        rendered = self._markdown.feed(delta)
        self.out.write(rendered)
        self._flush(self.out)
        self._track_screen_answer(rendered)
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
            if self._shared_terminal:
                self._screen_answer_open = False
                self._screen_answer_column = 0
        self.err.write(line + "\n")
        self._flush(self.err)
        if self._shared_terminal:
            # Chrome always owns complete lines.  The logical stdout
            # answer may still be open (and must stay exact in a pipe),
            # but visually the next answer byte starts on a fresh row.
            self._screen_answer_open = False
            self._screen_answer_column = 0

    def finish(self) -> None:
        """Close the answer stream, so a redirected stdout ends in a newline."""
        tail = self._markdown.flush()
        if tail:
            # Held-back bytes and any style still open. Without this a
            # two-character tail that looked like the start of a marker
            # would never be printed at all.
            self.out.write(tail)
            self._flush(self.out)
            self._track_screen_answer(tail)
        if self._out_open:
            self.out.write("\n")
            self._flush(self.out)
            self._track_screen_answer("\n")
            self._out_open = False
            self._break_pending = False

    @property
    def answer_open(self) -> bool:
        """True while the logical stdout answer ends mid-line.

        This is the file/output contract.  ``screen_answer_open`` is the
        separate cursor contract used by the live terminal composer.
        """
        return self._out_open

    @property
    def screen_answer_open(self) -> bool:
        """Whether the shared terminal cursor is after answer text.

        This differs deliberately from :attr:`answer_open`.  A chrome
        line may visually break a partial answer while stdout must remain
        byte-for-byte contiguous for redirection.  The bottom composer
        cares about the former; output correctness cares about the latter.
        """
        return self._shared_terminal and self._screen_answer_open

    @property
    def screen_answer_column(self) -> int:
        return self._screen_answer_column if self.screen_answer_open else 0

    @staticmethod
    def _streams_share_terminal(out, err) -> bool:
        try:
            if not (out.isatty() and err.isatty()):
                return False
        except Exception:
            return False
        try:
            import os
            return os.fstat(out.fileno()).st_rdev == os.fstat(err.fileno()).st_rdev
        except Exception:
            # String-backed model terminals used by tests have no fd.  If
            # both claim to be TTYs they model the ordinary shared screen.
            return True

    def _track_screen_answer(self, rendered: str) -> None:
        """Track the physical column reached by styled answer bytes.

        ANSI styling costs no columns; wide and combining characters do.
        Only the column is needed: the live composer is drawn immediately
        below the answer cursor and can therefore return with one relative
        move, even if painting it made the terminal scroll.
        """
        if not self._shared_terminal or not rendered:
            return
        from .repl_box import char_width

        width = max(2, int(self.width or 80))
        clean = rr.strip_control(rendered)
        for ch in clean:
            if ch == "\n":
                self._screen_answer_column = 0
                self._screen_answer_open = False
                continue
            if ch == "\r":
                self._screen_answer_column = 0
                continue
            if ch == "\t":
                step = 8 - (self._screen_answer_column % 8)
            else:
                step = char_width(ch)
            self._screen_answer_column = (
                self._screen_answer_column + step) % width
            self._screen_answer_open = True

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
                    tool_input=item.tool_input,
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


def permission_mode(engine) -> str:
    """The posture, or "" when this backend carries no permissions object.

    create_client builds KitToolPermissions for every chat-API provider
    (kit, ollama, openai). On the CLI backends there is none, so
    an empty answer here is a fact the banner has to state rather than
    paper over with a plausible-looking default.
    """
    try:
        perms = engine.kit_permissions
    except Exception:
        return ""
    mode = getattr(perms, "mode", "") if perms is not None else ""
    return mode if isinstance(mode, str) else ""


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
    tool_inputs: list[tuple[str, dict]] = []
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
        parsed = rr._as_dict(input_json)
        tool_calls.append({"name": name, "input": parsed})
        tool_inputs.append((name, parsed))
        _emit(RenderItem("tool_use", name=name, text=input_json or ""))

    def _take_tool_input(name: str) -> dict | None:
        """Match a result to the oldest outstanding call of that tool."""
        for index, (queued_name, parsed) in enumerate(tool_inputs):
            if queued_name == name:
                tool_inputs.pop(index)
                return parsed
        return None

    # One result, one line. The engine reports a tool result on two
    # callbacks — the head slice on on_tool_result, the verdict and the
    # true size on on_tool_result_meta — and rendering each as it arrives
    # drew every call twice, the second time with a body it did not have.
    # The text is held until the verdict catches up.
    pending: dict[str, Any] = {
        "name": "", "text": "", "held": False, "tool_input": None,
    }

    def _flush_pending() -> None:
        if pending["held"]:
            _emit(RenderItem("tool_result", name=pending["name"],
                             text=pending["text"],
                             tool_input=pending["tool_input"]))
            pending.update(name="", text="", held=False, tool_input=None)

    def _on_tool_result(name: str, output: str) -> None:
        # A second result before the first was completed means the meta
        # callback is not firing on this path; do not lose the first line.
        _flush_pending()
        pending.update(name=name, text=output or "", held=True,
                       tool_input=_take_tool_input(name))

    def _on_tool_result_meta(name: str, meta: dict) -> None:
        text = pending["text"] if pending["held"] and pending["name"] == name else ""
        tool_input = (pending["tool_input"]
                      if pending["held"] and pending["name"] == name
                      else _take_tool_input(name))
        if pending["held"] and pending["name"] != name:
            _flush_pending()
        pending.update(name="", text="", held=False, tool_input=None)
        _emit(RenderItem("tool_result", name=name, text=text, data=dict(meta),
                         tool_input=tool_input))

    def _on_denied(name: str) -> None:
        _take_tool_input(name)
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
    #: ``-n/--name``. The address an operator reaches this session by, and
    #: the title it shows in the list of open sessions.
    session_name: str = ""


def _is_a_note(text: str) -> bool:
    """True when ``text`` is a "#" note rather than a pasted document.

    The check used to be ``startswith("#")`` against the whole
    submission, and a brief pasted into a terminal arrives as ONE
    submission. A task beginning "# Gemeinsamer Teil" therefore went
    into the memory store instead of to the model -- done to five
    sessions at once on 2026-09-20, each of which printed a line about
    remembering and then sat at an empty prompt with its assignment
    gone.

    A note is one line. A document has more than one, and that is the
    whole distinction: "#" plus a few words is still a note, because
    guessing otherwise would take the feature away from the people who
    use it.
    """
    body = str(text or "")
    if not body.startswith("#"):
        return False
    if not body[1:].strip():
        return False                       # a lone "#" says nothing
    return len([ln for ln in body.strip().splitlines() if ln.strip()]) == 1


def _fit_to_width(text: str, width: int) -> str:
    """*text* cut to *width* display columns, never folded.

    Measured in columns rather than characters: the spinner and the
    notices carry emoji, and one of those is two columns wide. Counting
    characters let a line that "fit" wrap anyway, which is the fault this
    exists to stop.
    """
    from .repl_box import string_width

    if width <= 0:
        return ""
    # UI strings may carry trusted SGR colour.  It occupies no terminal
    # columns; counting its bytes both shortens a line unnecessarily and
    # can cut an escape sequence in half.  Keep styling when the visible
    # text fits, and fall back to safe plain text when truncation is needed.
    visible = rr.strip_control(text)
    if string_width(visible) <= width:
        return text
    cut = visible
    while cut and string_width(cut) > max(1, width - 1):
        cut = cut[:-1]
    return cut + "…"


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
                 out=None, err=None, read_line: Callable[[str], str] | None = None,
                 broker=None):
        self.engine = engine
        self.broker = broker
        self.opts = opts or ReplOptions()
        self.out = out if out is not None else sys.stdout
        self.err = err if err is not None else sys.stderr
        self.transcript = Transcript(
            self.out, self.err,
            theme=rr.theme_for(self.err, self.opts.color),
            show_tools=self.opts.show_tools,
            show_thinking=self.opts.show_thinking,
            color=self.opts.color,
        )
        self._read_line = read_line or self._input
        self._q: queue.Queue[RenderItem] = queue.Queue()
        self._turn_active = threading.Event()
        self._interrupts = 0
        self._idle_interrupts = 0
        self._prev_sigint = None
        self._prev_sigwinch = None
        # How much history was already on disk when this session started.
        self._hist_base = 0
        # Set by the resize handler, acted on by the loop. A signal handler
        # runs between bytecodes of whatever the main thread was doing, so
        # it records and the loop repaints.
        self._width_dirty = False
        self._stdin = sys.stdin
        # What the user is typing WHILE a turn runs, and the last tool
        # result, so Ctrl+O has something to expand.
        self._input_line = ""
        self._input_cursor = 0
        self._last_result_text = ""
        # The live composer is one owned region: progress above the same
        # two-rule input box and key hint used while idle. Input and status
        # used to compete for one row, so an empty input vanished behind
        # the spinner and typing hid all progress info.
        self._bottom = ""
        self._bottom_cursor = (-1, -1)
        self._bottom_rows = 0
        self._bottom_below = 0
        self._bottom_anchor_gap = 0
        self._bottom_anchor_column = 0
        self._show_tasks = False
        self._turn_t0 = 0.0
        self._turn_base = (0, 0, 0.0)
        self._spin = 0
        self._last_paint = 0.0
        # Typed during a turn: queued, never injected. A queued message
        # cannot be lost and cannot land in a context nobody could see.
        self.queued: list[str] = []
        self._quit = False

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
            # What the third one does: raise here, unwind through turn()'s
            # finally with a bounded join, and return 130 from the loop. It
            # cannot cancel the call, so it does not claim to.
            self._q.put(RenderItem("notice", text=(
                "Still waiting on the running tool call. Once more leaves "
                "this session without it: the call is not cancelled, and "
                "anything it started keeps running.")))
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
    def _begin_turn_resets_the_abort(self) -> None:
        """An abort belongs to the turn it ended, not to the session.

        [a] abort refuses everything still in flight -- "refuse
        everything that arrives after". Without this reset the refusal
        posture survived into every later turn, and a session that
        aborted ONCE answered every question for hours with a silent
        "user denied" nobody had decided. The next turn starts from a
        clean posture: its questions are asked.
        """
        if self.broker is not None:
            try:
                self.broker.reset_abort()
            except Exception:
                pass

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

        import time as _time

        self._begin_turn_resets_the_abort()
        self._interrupts = 0
        self._turn_t0 = _time.monotonic()
        self._last_paint = 0.0
        # The engine's counters are cumulative for the SESSION, so the live
        # line differences against this baseline. Without it the first turn
        # looks right and every later one reports the whole conversation.
        try:
            _st = self.engine.get_status() or {}
            self._turn_base = (int(_st.get("input_tokens", 0) or 0),
                               int(_st.get("output_tokens", 0) or 0),
                               float(_st.get("cost_usd", 0.0) or 0.0))
        except Exception:
            self._turn_base = (0, 0, 0.0)
        compaction_before = getattr(self.engine, "last_compaction_info", None)

        from . import repl_keys as rk

        self._turn_active.set()
        worker = threading.Thread(target=_worker, name="agent-turn", daemon=True)
        worker.start()
        abandoning = False
        leaving = False
        try:
            self._pump(worker)
        except KeyboardInterrupt:
            # The last step of the ladder. It travels through the finally
            # below, so the settle there must not wait on the running call.
            abandoning = True
            raise
        except rk.TerminalLeft:
            # No next turn will run, so there is nothing to settle for: the
            # worker is told to stop and left behind, and nothing is written
            # to a terminal that may no longer take it.
            leaving = True
            self._stop_for_leaving()
            raise
        finally:
            self._turn_active.clear()
            if not leaving:
                self._settle(worker, abandon=abandoning)
        self.transcript.finish()
        self._report_compaction(compaction_before)
        self._report_tasks()
        self._offer_next_steps()
        self._report_status()
        result = result_box[0] if result_box \
            else TurnResult(error="turn produced nothing")
        # For the length-continuation: a length end WITH a tool call is
        # the tool loop's to carry on, not ours.
        self.__dict__["_turn_used_tools"] = bool(result.tool_calls)
        return result

    def _pump(self, worker: threading.Thread) -> None:
        from . import repl_keys as rk

        with rk.RawMode(self._stdin) as raw:
            decoder = rk.KeyDecoder()
            while True:
                if self.broker is not None:
                    pending = self.broker.take()
                    if pending is not None:
                        self._answer(pending, raw)
                        continue
                if self._width_dirty:
                    # The handler recorded the new size; the repaint
                    # happens here, on the thread that owns the terminal.
                    self._width_dirty = False
                    # Erase with the OLD geometry, then ask for the new
                    # width before laying the composer out again.  Reversing
                    # that order sizes the erase from a picture that is not
                    # on screen; omitting refresh_width repaints forever at
                    # the launch-time width after SIGWINCH.
                    self._clear_bottom()
                    self.transcript.refresh_width()
                    self._repaint_bottom(force=True)
                if raw.active:
                    for event in decoder.feed(raw.read_ready(_PUMP_TICK_S)):
                        self._on_key(event, decoder)
                    self._repaint_bottom()
                    # Operator mail, mid-turn: the dashboard steers a
                    # message into the running turn between rounds
                    # (push_steer); the terminal dropped it at the
                    # prompt instead, so a long turn was deaf to its
                    # supervisor for exactly as long as it ran
                    # (assignment 8, point b, 2026-09-22). Same
                    # throttle as the idle wake, same rendering (it
                    # reads as coming from another session, never as
                    # user input).
                    try:
                        self._steer_operator_mail()
                    except Exception:
                        pass
                    # The heartbeat of a session IN a turn. Presence was
                    # renewed only at the prompt and on a broker question,
                    # so a turn longer than session_presence._STALE_S --
                    # one model round after another, no question asked --
                    # dropped the session out of open_sessions(), and a
                    # peer got "no other open session 'runde2c-s4'" while
                    # s4 was working (measured 2026-09-21). The pump runs
                    # for the whole turn, so the refresh travels with it;
                    # announce itself throttles unchanged records to one
                    # write per heartbeat, never one per read.
                    try:
                        self._announce_presence()
                    except Exception:
                        pass
                try:
                    item = self._q.get(
                        timeout=0.0 if raw.active else _PUMP_TICK_S)
                except queue.Empty:
                    if not worker.is_alive() and self._q.empty():
                        self._clear_bottom()
                        return
                    continue
                if item.kind == "done":
                    self._clear_bottom()
                    return
                self._render_around_bottom(item)

    # -- keys during a turn ----------------------------------------------
    def _on_key(self, event, decoder) -> None:
        from . import repl_keys as rk

        if event.kind == rk.EOF:
            # ctrl+d on an empty line leaves, during a turn as at the
            # prompt. It is also all `script` passes down when the terminal
            # above it goes away, so dropping it here kept a session
            # running with nobody left to see it.
            raise rk.TerminalLeft("ctrl+d", 0)
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
                self._show_user_input(text, queued=True)
                self.queued.append(text)
                self._recall_later(text)
                self.transcript.chrome(self.transcript.theme.dim(
                    f"queued ({len(self.queued)}) — goes out when this turn "
                    "ends, or ctrl+g to send it into this one"))
            return
        if event.kind == rk.STEER:
            self._clear_input_line()
            if event.text.strip():
                self._show_user_input(event.text.strip())
            self._steer(event.text.strip())
            return
        if event.kind == rk.CYCLE_MODE:
            self._clear_input_line()
            self._cycle_mode()
            self._draw_input_line(decoder.buffer, cursor=decoder.cursor)
            return
        if event.kind == rk.EXPAND:
            self._clear_input_line()
            self._expand_last_result()
            self._draw_input_line(decoder.buffer, cursor=decoder.cursor)
            return
        if event.kind == rk.TASKS:
            self._show_tasks = not self._show_tasks
            self._clear_bottom()
            self.transcript.chrome(self.transcript.theme.dim(
                f"task list {'on' if self._show_tasks else 'off'}"))
            self._repaint_bottom(force=True)
            return
        if event.kind == rk.REDRAW:
            self._clear_bottom()
            # Ctrl+L is also the recovery key when a multiplexer swallowed
            # SIGWINCH.  A redraw at the stale width is not a recovery.
            self.transcript.refresh_width()
            self._repaint_bottom(force=True)
            return
        if event.kind == rk.EDIT:
            self._draw_input_line(event.text, cursor=decoder.cursor)

    def _steer(self, text: str) -> None:
        """Put what was typed into the turn that is running.

        The engine already carries the whole mechanism: the client's steer
        inbox is drained between tool rounds AND at a turn that would
        otherwise end, so a message sent late still continues the turn
        rather than being answered next time. What was missing was a way
        to reach it from a keyboard.

        A backend without the inbox is told about, never worked around:
        `steer` returns False there, and the message becomes a queued one
        with the reason on screen. Dropping it would lose something the
        user typed, and pretending it landed would be worse.
        """
        # Stripped here, not only at the call site. `push_steer` drops a
        # blank of its own accord, so an unstripped one would send nothing
        # and still print that it had — a line on screen asserting an
        # effect that did not happen.
        text = (text or "").strip()
        if not text:
            return
        try:
            delivered = bool(self.engine.steer(text))
        except Exception:
            delivered = False
        self._recall_later(text)
        if delivered:
            self.transcript.chrome(self.transcript.theme.cyan(
                "→ sent into the running turn"))
            return
        self.queued.append(text)
        self.transcript.chrome(self.transcript.theme.dim(
            f"this backend takes no message mid-turn — queued "
            f"({len(self.queued)}) instead"))

    def _posture_now(self) -> str:
        """The approval posture, or "" where this backend has no gate.

        Read fresh each time the box is drawn: shift+tab changes it, and
        a row that shows the posture you had two turns ago is worse than
        one that shows none.
        """
        try:
            return str(getattr(self.engine.kit_permissions, "mode", "") or "")
        except Exception:
            return ""

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

    # -- what a typed line means ------------------------------------------
    def _handle_line(self, line: str) -> str:
        """Returns the text to send to the model, or "" when it was handled.

        Order matters and is the dashboard's: a builtin, then a user
        command, then a skill, then a subagent command, then the model.
        The same `/deploy` has to mean the same thing on both surfaces.
        """
        from . import repl_commands as rc

        if line.startswith(rc.SHELL_PREFIX) and len(line) > 1:
            self._shell_out(line[1:].strip())
            return ""
        if _is_a_note(line):
            self._remember(line[1:].strip())
            return ""

        result = rc.dispatch(line, self._ctx())
        if not result.handled:
            return rc.expand_at_references(line, self.opts.cwd)
        if result.output:
            for out_line in result.output.splitlines():
                self.transcript.chrome(out_line)
        if result.clear:
            self._start_fresh()
        if result.quit:
            self._quit = True
        return result.prompt

    def _ctx(self):
        from . import repl_commands as rc
        return rc.ReplContext(
            engine=self.engine, workspace=self.opts.cwd,
            session_id=str(getattr(self.engine, "session_id", "") or ""))

    def _start_fresh(self) -> None:
        try:
            self.engine.messages.clear()
        except Exception:
            pass

    def _shell_out(self, command: str) -> None:
        """`!cmd` runs through the AGENT's gate, not around it.

        A shell escape that skipped the deny-list and the approval prompt
        would be a way to do from the agent's prompt exactly what the
        agent may not do, which is worse than not having the affordance.
        The output joins the session, so the model can see what happened.
        """
        if not command:
            return
        executor = getattr(self.engine, "run_gated_bash", None)
        if not callable(executor):
            self.transcript.chrome(self.transcript.theme.dim(
                "! is not available on this backend — it runs through the "
                "same gate the agent uses, and this one has none"))
            return
        out = self._off_thread(lambda: executor(command))
        if isinstance(out, Exception):
            self.transcript.chrome(
                self.transcript.theme.red(f"[error] {out}"))
            return
        for line in self._shell_lines(out):
            self.transcript.chrome("  " + line)

    @staticmethod
    def _shell_lines(out) -> list[str]:
        """What the user asked to see, not the envelope it arrived in.

        The tool returns a JSON result because that is what a MODEL needs;
        a person typing `!git status` wants the output of git status.
        """
        import json

        text = str(out or "")
        try:
            payload = json.loads(text)
        except (TypeError, ValueError):
            payload = None
        if isinstance(payload, dict):
            if payload.get("error"):
                return [f"refused: {payload['error']}"]
            parts = [str(payload.get("stdout", "") or ""),
                     str(payload.get("stderr", "") or "")]
            text = "\n".join(p for p in parts if p.strip())
            code = payload.get("exit_code")
            if code not in (0, None):
                text += f"\n(exit {code})"
        lines = rr.strip_control(text).splitlines()
        if len(lines) > 200:
            lines = lines[:200] + [f"… {len(lines) - 200} more lines"]
        return lines

    def _off_thread(self, work):
        """Run something that can hit the gate, and stay able to answer it.

        Anything routed through the permission gate may ask for approval,
        and the broker parks that question for THE MAIN THREAD to answer.
        Calling such work directly from the main thread is therefore a
        deadlock: the asker waits for the answerer, which is itself. So it
        runs on a worker while this thread keeps answering — the same
        arrangement a turn already uses, for the same reason.
        """
        box: list = []

        def _run():
            try:
                box.append(work())
            except Exception as exc:      # noqa: BLE001
                box.append(exc)

        worker = threading.Thread(target=_run, name="agent-gated", daemon=True)
        worker.start()
        while worker.is_alive() or not box:
            if self.broker is not None:
                pending = self.broker.take()
                if pending is not None:
                    # Entered for the prompt and released again as soon as
                    # it is answered: the span that needs cbreak is the
                    # keystroke, not the command that raised it.
                    raw = self._raw_for_prompt()
                    try:
                        self._answer(pending, raw)
                    finally:
                        self._release_prompt_raw(raw)
                    continue
            worker.join(timeout=_PUMP_TICK_S)
            if not worker.is_alive():
                break
        return box[0] if box else None

    def _raw_for_prompt(self):
        """A keystroke reader for a prompt raised outside a turn.

        Entered here, released by ``_release_prompt_raw``, and the pairing
        is load-bearing twice: a reader that is never left holds the
        terminal in cbreak with echo off, so the next idle ``input()``
        shows nothing of what is typed — and every ``RawMode.__enter__``
        registers an atexit hook that only ``restore`` takes back, so an
        unreleased one leaves a hook behind per approval.
        """
        from . import repl_keys as rk
        raw = rk.RawMode(self._stdin)
        return raw.__enter__()

    @staticmethod
    def _release_prompt_raw(raw) -> None:
        """Give the terminal back, whatever the reader turned out to be."""
        restore = getattr(raw, "restore", None)
        if callable(restore):
            try:
                restore()
            except Exception:
                pass

    def _remember(self, text: str) -> None:
        """`#note` writes a memory, marked as the user's own."""
        if not text:
            return
        try:
            from . import memory_store
            # The names the store actually takes. It was called with
            # kind=/workspace=/author= and raised on every single note,
            # so "#" had never once worked on this path -- the error
            # about remembering was the only thing it ever produced.
            memory_store.save_typed_memory(
                text=text, memory_type="user", repo_root=self.opts.cwd,
                scope="project", source=memory_store.SOURCE_USER)
            self.transcript.chrome(self.transcript.theme.dim("remembered"))
        except Exception as exc:
            self.transcript.chrome(
                self.transcript.theme.red(f"could not remember: {exc}"))

    # -- approvals --------------------------------------------------------
    def _answer(self, req, raw) -> None:
        """Answer one request -- or refuse it, when the terminal leaves.

        A request taken off the queue is out of reach of ``abort_all``, so
        a session leaving from inside the dialog refuses it here. Left
        alone, the asking thread would wait out the request's expiry.
        """
        from . import repl_keys as rk

        try:
            self._answer_request(req, raw)
        except rk.TerminalLeft:
            try:
                if self.broker is not None:
                    self.broker.resolve(req, self._refuse(req))
            except Exception:
                pass
            raise

    def _answer_request(self, req, raw) -> None:
        """Render one request and read the answer. Main thread only.

        The worker is blocked inside the gate waiting on this, which is
        what makes reading from stdin here safe: nothing else can be
        writing to the terminal, and nothing else is reading from it.
        """
        from . import terminal_confirm as tc

        self._clear_bottom()
        if req.kind == tc.ASK:
            self._answer_question(req, raw)
            return

        suggestion = (self.broker.suggest_pattern(req.command)
                      if req.is_shell and req.command else "")
        options = tc.options_for(req, suggestion=suggestion)
        self.transcript.chrome(tc.render_request(
            req, theme=self.transcript.theme, width=self.transcript.width,
            options=options))

        allowed = {o.key for o in options}
        while True:
            # Answered from outside (approvals approve/deny) while the
            # dialog was waiting: it is OVER. Reading on would take the
            # keys of the next prompt line -- a typed command acted as
            # [a] abort of a question long since answered, and the turn
            # it belonged to never ran. Check before every key: the
            # answer can land in the middle of the loop too.
            if getattr(req, "resolved", False):
                self._say_ended_elsewhere(req)
                return
            key = self._read_key(raw, allowed | {"\x1b"},
                                 stop=lambda: getattr(req, "resolved", False))
            if not key:
                continue                # answered while waiting: see above
            if key in ("\x1b", "n"):
                if self._apply(req, self._refuse(req)):
                    self.transcript.chrome(
                        self.transcript.theme.dim("  refused"))
                return
            if key == "?":
                self.transcript.chrome(tc.render_help(req, options))
                continue
            if key == "a":
                # This one first, the rest after: abort_all resolves
                # everything still queued, so asking it first would make
                # the answer to THIS request late by construction and the
                # line below would report an expiry that never happened.
                landed = self._apply(req, self._refuse(req))
                denied = self.broker.abort_all()
                self._stop_engine()
                extra = (f" and {len(denied)} other request(s) in flight"
                         if denied else "")
                head = ("! aborted — this was refused" if landed
                        else "! aborted")
                self.transcript.chrome(self.transcript.theme.yellow(
                    f"{head}{extra}, and the turn is ending"))
                return
            if key == "y":
                self._apply(req, self._allow(req))
                return
            if key == "d" and req.kind == tc.PLAN:
                if self._apply(req, {"approved": True, "new_mode": "default"}):
                    self.transcript.chrome(self.transcript.theme.cyan(
                        "approval → default"))
                return
            if key == "e" and req.kind == tc.PLAN:
                if self._apply(req,
                               {"approved": True, "new_mode": "acceptEdits"}):
                    self.transcript.chrome(self.transcript.theme.cyan(
                        "approval → acceptEdits"))
                return
            if key == "e":
                # The mode switch is a posture for the session, so it is
                # reported on its own terms — it happens whether or not the
                # answer to this request still had somewhere to land.
                ok, msg = self.broker.accept_edits()
                self.transcript.chrome(self.transcript.theme.cyan(f"  {msg}"))
                self._apply(req, self._allow(req))
                return
            if key in ("A", "k"):
                pattern = (self.broker.exact_pattern(req.command) if key == "A"
                           else suggestion)
                if self._confirm_persist(pattern, raw):
                    ok, msg = self.broker.persist(pattern)
                    self.transcript.chrome(
                        self.transcript.theme.cyan(f"  {msg}") if ok
                        else self.transcript.theme.red(f"  {msg}"))
                self._apply(req, self._allow(req))
                return
            # Unreachable in practice: _read_key only ever returns a key
            # from `allowed`. Kept as a hard stop rather than a fallthrough
            # so a future key added to the menu and not to this chain
            # cannot silently mean "yes".
            self.transcript.chrome(self.transcript.theme.dim(
                "  that key does nothing here"))

    def _say_ended_elsewhere(self, req) -> None:
        """Close a dialog whose question something else already ended.

        An expiry resolves the request too, and it is not an answer from
        anywhere: say which of the two ended it. Nothing is applied -- the
        request is resolved, and a second answer would only be discarded
        and reported as late.
        """
        self._clear_bottom()
        if getattr(req, "expired", False):
            self.transcript.chrome(self.transcript.theme.yellow(
                "  too late — that request expired and was refused "
                "without you"))
        else:
            self.transcript.chrome(self.transcript.theme.dim(
                "  answered elsewhere"))

    def _apply(self, req, decision) -> bool:
        """Hand the answer to the broker, and say when it arrived too late.

        ``resolve`` discards an answer to a request that had already
        expired or been aborted — the gate was told "no" on the user's
        behalf and the model has moved on. Printing "refused" or
        "approval → default" regardless puts an effect on the screen that
        nothing applied, so every line that describes an answer is printed
        only when this returned True.
        """
        if self.broker.resolve(req, decision):
            return True
        self.transcript.chrome(self.transcript.theme.yellow(
            "  too late — that request had already expired and was refused "
            "without you, so this key changed nothing"))
        return False

    def _confirm_persist(self, pattern: str, raw) -> bool:
        """A second keystroke, against the consequence spelled out.

        Both facts surprise people, so both are on screen: it is GLOBAL
        (the merge rule takes allow-patterns from the user file only, and
        ignores a repository's) and it is PERMANENT.
        """
        if not pattern:
            return False
        self.transcript.chrome(
            f"  persist  {pattern}\n"
            "  Every future command matching it runs without asking — in "
            "this project and every other, in this session and every "
            "future one.")
        self.transcript.chrome("  [y] persist   [any other key] just this once")
        return self._read_key(raw, {"y", "n", "\x1b"}) == "y"

    def _answer_question(self, req, raw) -> None:
        """ask_user_question: numbered options, one key."""
        payload = req.payload or {}
        question = str(payload.get("question", "") or "(no question given)")
        options = [str(o) for o in (payload.get("options") or [])][:9]
        self.transcript.chrome(self.transcript.theme.bold(f"? {question}"))
        for i, opt in enumerate(options, 1):
            self.transcript.chrome(f"  {i}. {rr.strip_control(opt)}")
        if not options:
            self._apply(req, {"answers": []})
            return
        allowed = {str(i) for i in range(1, len(options) + 1)} | {"\x1b"}
        while True:
            # Same as _answer_request: an answer from outside ends the
            # dialog, and later keys belong to the prompt that follows.
            if getattr(req, "resolved", False):
                self._say_ended_elsewhere(req)
                return
            key = self._read_key(raw, allowed,
                                 stop=lambda: getattr(req, "resolved", False))
            if not key:
                continue                # answered while waiting
            if key == "\x1b":
                self._apply(req, {"answers": []})
                return
            self._apply(req, {"answers": [options[int(key) - 1]]})

    def _read_key(self, raw, allowed: set[str], *, stop=None) -> str:
        """One keystroke, from the reader the key layer already owns.

        ``stop`` is asked between the 0.1 s reads; when it says the
        question is over, the wait ends with "" and no key is taken. The
        dialog checked for an outside answer only BEFORE it started to
        wait, so an answer that landed during the wait -- the usual case,
        a supervisor answers a dialog that is already on screen -- left
        it reading on: the session never got back to its prompt, mail
        waiting there was never delivered, and the next key typed into
        the pane went to a question long since answered.
        """
        from . import repl_keys as rk

        if raw is not None and getattr(raw, "active", False):
            while True:
                if stop is not None and stop():
                    return ""
                chunk = raw.read_ready(0.1)
                if not chunk:
                    continue
                for ch in chunk:
                    if ch == "\x04":
                        # ctrl+d is leaving, not an unclaimed key: a dialog
                        # that ignored it waited on a terminal that was gone.
                        raise rk.TerminalLeft("ctrl+d", 0)
                    if ch in allowed:
                        return ch
                if chunk[0] == "\x1b":
                    return "\x1b"
                # Enter, arrows, anything unclaimed: ignored, never a
                # default. A default-yes turns approval into a rhythm, and
                # rhythm is what this whole layer exists to break.
        # No terminal to read from: an unanswerable prompt is a refusal,
        # never a silent yes — and the refusal has to be a key the caller
        # offered, because that is the chain it goes back into. A hardcoded
        # "n" reached `int(key)` on the question path, where the keys are
        # digits: the prompt raised instead of being answered and the
        # asking thread waited on a request nobody would resolve.
        for refusal in ("n", "\x1b"):
            if refusal in allowed:
                return refusal
        # Every caller offers one of those two. If one ever does not, Esc
        # still reads as "refuse" in each chain rather than as a choice.
        return "\x1b"

    @staticmethod
    def _allow(req):
        from . import terminal_confirm as tc
        return True if req.kind == tc.CONFIRM else req.decision

    @staticmethod
    def _refuse(req):
        from . import terminal_confirm as tc
        if req.kind == tc.PLAN:
            return {"approved": False, "new_mode": "plan"}
        if req.kind == tc.ASK:
            return {"answers": []}
        return False

    # -- the live bottom composer -----------------------------------------
    def _can_redraw(self) -> bool:
        """Only a terminal gets cursor control.

        Redrawing needs \r and an erase sequence, and those are exactly
        the characters this codebase strips out of tool output before
        printing it. Writing them into a redirected stderr would put
        control codes in a log file for a line nobody can see anyway.
        """
        try:
            return bool(self.err.isatty())
        except Exception:
            return False

    def _set_bottom(self, text: str, status: str | None = None, *,
                    cursor: int | None = None) -> None:
        """Paint the live composer without surrendering the answer cursor.

        The active form is the idle prompt with progress above it::

            ⠴ 12s  model  mode  ↑0 ↓0  esc to interrupt
            ─────────────────────────────
            > the next message, wrapped onto as many rows as it needs
            ─────────────────────────────
              mode · esc interrupt · shift+tab approval mode · /help

        ``status=None`` retains the small two-row primitive used by a few
        callers and screen-model tests.  In the active form the cursor is
        returned to the input row, with the status and help still visible.

        If a model sentence is only half streamed, the composer starts on
        the following row.  Clearing it later walks back to the exact
        answer column before the next stdout delta is written.  Relative
        movement is intentional: unlike save/restore cursor, it remains
        correct when drawing the composer scrolls a terminal at its edge.
        """
        width = max(2, int(getattr(self.transcript, "width", 0) or 80))
        if status is None:
            # The compatibility primitive is still deliberately one line.
            # The live path below owns wrapping and its dynamic row count.
            rule = "─" * (width - 1)
            line = _fit_to_width(text, width - 1)
            rows_on_screen = [rule, line]
            cursor_row = 1
            cursor_col = min(self._display_width(line), width - 1)
        else:
            from . import repl_box as rb

            raw = text or ""
            at = len(raw) if cursor is None else int(cursor)
            at = max(0, min(at, len(raw)))
            view = rb.render_box(
                raw, at, width, _box_hint(self._posture_now()))
            box_rows = list(view.rows)
            if view.hint_row is not None:
                box_rows[view.hint_row] = self.transcript.theme.dim(
                    box_rows[view.hint_row])
            rows_on_screen = [_fit_to_width(status, width - 1), *box_rows]
            # BoxView.cursor is relative to the CONTENT, while the full
            # form has a rule before its content.  The narrow fallback has
            # only its content row and therefore no offset.
            content_offset = 0 if len(view.rows) == 1 else 1
            cursor_row = 1 + content_offset + view.cursor[0]
            cursor_col = view.cursor[1] + (2 if view.border else 0)
            cursor_col = min(max(0, cursor_col), width - 1)

        state = "\n".join(rows_on_screen)
        cursor_state = (cursor_row, cursor_col)
        if (state == self._bottom
                and cursor_state == getattr(self, "_bottom_cursor", None)):
            return
        self._clear_bottom()
        self._bottom = state
        self._bottom_cursor = cursor_state
        if not self._can_redraw():
            return

        answer_open = bool(getattr(
            self.transcript, "screen_answer_open", False))
        answer_col = int(getattr(
            self.transcript, "screen_answer_column", 0) or 0)
        out = ["\r\n" if answer_open else ""]
        for i, row in enumerate(rows_on_screen):
            if i:
                out.append("\r\n")
            out.extend(("\r\x1b[K", row))
        rows = len(rows_on_screen)
        below = rows - 1 - cursor_row
        natural_col = self._display_width(rows_on_screen[-1])
        if below or cursor_col != natural_col:
            if below:
                out.append(f"\x1b[{below}A")
            out.append("\r")
            if cursor_col:
                out.append(f"\x1b[{cursor_col}C")
        self.err.write("".join(out))
        self._bottom_rows = rows
        self._bottom_below = below
        self._bottom_anchor_gap = 1 if answer_open else 0
        self._bottom_anchor_column = min(max(0, answer_col), width - 1)
        self._flush_err()

    def _clear_bottom(self) -> None:
        """Erase the rows this zone drew, and leave the cursor on the first.

        The count is kept rather than assumed. The transcript continues
        from there, so a line rendered next lands on the rule's row rather
        than below a rule nobody erased.
        """
        if not self._bottom:
            return
        self._bottom = ""
        self._bottom_cursor = (-1, -1)
        rows = int(getattr(self, "_bottom_rows", 0) or 0)
        below = int(getattr(self, "_bottom_below", 0) or 0)
        anchor_gap = int(getattr(self, "_bottom_anchor_gap", 0) or 0)
        anchor_col = int(getattr(self, "_bottom_anchor_column", 0) or 0)
        self._bottom_rows = 0
        self._bottom_below = 0
        self._bottom_anchor_gap = 0
        self._bottom_anchor_column = 0
        if not self._can_redraw() or rows <= 0:
            return
        out = []
        if below:
            out.append(f"\x1b[{below}B")
        out.append("\r\x1b[K")
        out.append("\x1b[1A\r\x1b[K" * (rows - 1))
        if anchor_gap:
            out.append("\x1b[1A\r")
            if anchor_col:
                out.append(f"\x1b[{anchor_col}C")
        self.err.write("".join(out))
        self._flush_err()

    @staticmethod
    def _display_width(text: str) -> int:
        """Columns in a trusted UI row (ANSI styling costs none)."""
        from .repl_box import string_width
        return string_width(rr.strip_control(text or ""))

    def _show_user_input(self, text: str, *, queued: bool = False) -> None:
        """Keep a submitted message in the transcript after its box clears.

        The editable box is transient by design.  Without this permanent
        copy, pressing Enter made the user's task disappear precisely when
        a slow first turn left them waiting minutes for an answer.  Control
        sequences are stripped because pasted text is still untrusted
        terminal input; line structure remains so a pasted traceback is
        recognisable as the message that was sent.
        """
        cleaned = rr.strip_control(text or "").strip()
        if not cleaned:
            return
        lines = cleaned.splitlines() or [cleaned]
        marker = self.transcript.theme.bold("» ")
        continuation = "  "
        shown = [marker + lines[0]]
        shown.extend(continuation + line for line in lines[1:])
        if queued:
            shown[-1] += self.transcript.theme.dim("  (queued)")
        self.transcript.chrome("\n".join(shown))

    def _teardown_screen(self) -> None:
        """Give the screen back to the shell, on a row of its own.

        Leaving used to draw nothing: the prompt stayed where it was and
        the shell wrote its own into the middle of it --

            (.venv) [user@host]$ ^C────────────────────────────
            (.venv) [user@host]$  mode · /help

        Erasing from the cursor to the end of the screen takes the input
        row and the hint under it; the rule above stays in the scrollback,
        where it marks the end of the session rather than getting in the
        way of the next command.
        """
        self._clear_bottom()
        if not self._can_redraw():
            return
        try:
            self.err.write("\r\x1b[J")
            self._flush_err()
        except Exception:
            pass

    def _draw_input_line(self, text: str, *, cursor: int | None = None) -> None:
        """Record the draft and its cursor, then repaint its whole box."""
        self._input_line = text or ""
        at = len(self._input_line) if cursor is None else int(cursor)
        self._input_cursor = max(0, min(at, len(self._input_line)))
        self._repaint_bottom(force=True)

    def _clear_input_line(self) -> None:
        self._input_line = ""
        self._input_cursor = 0
        self._clear_bottom()

    def _status_line(self) -> str:
        """Elapsed, model, posture, what this turn has cost so far.

        Read from get_status() and differenced against a baseline taken at
        turn start, because the engine's counters are cumulative for the
        session — reporting them raw would show the whole conversation's
        cost as this turn's.
        """
        import time as _time

        self._spin = (self._spin + 1) % len(_SPINNER)
        elapsed = max(0.0, _time.monotonic() - self._turn_t0)
        try:
            status = self.engine.get_status() or {}
        except Exception:
            status = {}
        base_in, base_out, base_cost = self._turn_base
        tin = max(0, int(status.get("input_tokens", 0) or 0) - base_in)
        tout = max(0, int(status.get("output_tokens", 0) or 0) - base_out)
        cost = max(0.0, float(status.get("cost_usd", 0.0) or 0.0) - base_cost)
        model = str(getattr(getattr(self.engine, "client", None), "model", "")
                    or "?")
        mode = permission_mode(self.engine) or str(status.get("mode", "") or "")

        bits = [f"{_SPINNER[self._spin]} {elapsed:4.0f}s", model]
        if mode:
            bits.append(mode)
        bits.append(f"↑{tin} ↓{tout}")
        if cost > 0:
            bits.append(f"${cost:.4f}")
        bits.append("esc to interrupt")
        line = "  ".join(bits)
        return self.transcript.theme.dim(
            rr.truncate_middle(line, max(20, self.transcript.width - 1)))

    def _repaint_bottom(self, *, force: bool = False) -> None:
        """Keep input and turn information visible at the same time."""
        import time as _time

        if not self._can_redraw():
            return
        now = _time.monotonic()
        if not force and (now - self._last_paint) < _STATUS_TICK_S:
            return
        self._last_paint = now
        if self._turn_active.is_set():
            self._set_bottom(
                self._input_line, self._status_line(),
                cursor=getattr(self, "_input_cursor", len(self._input_line)))
        elif self._input_line:
            # Defensive/legacy path: _draw_input_line is a turn-time API,
            # but a direct caller still gets a visible editable row.
            self._set_bottom(self._typed_row(self._input_line))

    def _typed_row(self, text: str) -> str:
        """The line being typed, cut to one screen row like the status row.

        `_clear_bottom` erases with a single `\\r\\x1b[K`, which reaches one
        line — so a row that ran past COLUMNS and wrapped left its first
        half stranded in the transcript. The END is what survives the cut
        here rather than the middle: that is where the cursor is, and a
        person typing has to see the characters they are typing.
        """
        width = max(20, self.transcript.width - 1)
        row = "> " + (text or "")
        if len(row) <= width:
            return row
        return "…" + row[len(row) - (width - 1):]

    def _render_around_bottom(self, item: RenderItem) -> None:
        """Rendering must not tear whatever is on the bottom row.

        Erase it, write the transcript line, put it back — the classic
        bottom-line problem, and the reason typing during a turn is worth
        having rather than merely possible.
        """
        held_input = self._input_line
        held_cursor = getattr(self, "_input_cursor", len(held_input))
        self._clear_bottom()
        if item.kind == "tool_result" and item.text:
            self._last_result_text = item.text
        self.transcript.render(item)
        self._input_line = held_input
        self._input_cursor = held_cursor
        self._repaint_bottom(force=True)

    def _settle(self, worker: threading.Thread, *,
                abandon: bool = False) -> None:
        """Join, THEN clear the stop. Both halves are load-bearing.

        The turn gate refuses a second concurrent turn by RETURNING a
        sentence rather than raising, so returning to the prompt while a
        stopped worker is still unwinding makes the next turn render
        machinery speech as the model's answer. And ``clear_stop`` refuses
        while the owning turn is in flight, so clearing before the join is
        a silent no-op that leaves the brake armed for the next turn.

        ``abandon`` is the last step of the interrupt ladder, and it is the
        one case where waiting is the wrong answer: the thread being joined
        is running the tool call the user is escaping, so the join gets a
        bound and what is left is said out loud rather than waited out.
        """
        worker.join(timeout=_JOIN_ABANDON_S if abandon
                    else _JOIN_NOTICE_AFTER_S)
        if worker.is_alive() and abandon:
            self._notice(
                "Leaving the running tool call behind — it was asked to stop "
                "and has not yet, and anything it started keeps running.")
        elif worker.is_alive():
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

    # -- what a finished turn has to say ----------------------------------
    def _report_compaction(self, before) -> None:
        """A long session silently losing its early history is worth a line."""
        after = getattr(self.engine, "last_compaction_info", None)
        if not after or after == before:
            return
        try:
            line = self.engine._compaction_status_line()
        except Exception:
            return
        if line:
            self.transcript.chrome(self.transcript.theme.dim(
                line.lstrip("- ")))

    def _report_status(self) -> None:
        """The user's own status line, if they configured one.

        Reused rather than re-derived: it already never raises, and it
        already refuses to run a workspace-supplied command.

        The "if" is enforced here. `render_status_line` falls back to a
        built-in template when nothing is configured, so this printed a
        line after every turn whose two fields — mode and branch — the
        banner and the live turn line already carry. A sentence in a
        docstring is not a condition; the condition is.
        """
        try:
            from . import status_line as sl
            if not sl.has_custom_status_line(self.opts.cwd):
                return
            status = self.engine.get_status() or {}
            text = sl.render_status_line(sl.StatusContext(
                workspace=self.opts.cwd,
                mode=permission_mode(self.engine) or str(status.get("mode", "")),
                model=str(getattr(getattr(self.engine, "client", None),
                                  "model", "") or ""),
                tokens=int(status.get("input_tokens", 0) or 0)
                + int(status.get("output_tokens", 0) or 0),
                cost_usd=float(status.get("cost_usd", 0.0) or 0.0),
            ))
        except Exception:
            return
        if text:
            self.transcript.chrome(self.transcript.theme.dim(text))

    def _offer_next_steps(self) -> None:
        """What could come next, offered as numbers the user can type.

        Drawn from the tasks the agent itself left open, so a suggestion
        is never a guess about what the user wants and costs no tokens
        to produce. A turn that closed everything offers nothing, which
        is the honest answer to "what now".

        Never raises: a courtesy that can take the prompt with it is not
        a courtesy.
        """
        self._next_steps = []
        try:
            from . import task_ticker
            steps = task_ticker.next_steps(
                self.opts.cwd,
                session_id=str(getattr(self.engine, "session_id", "") or ""))
        except Exception:
            return
        if not steps:
            return
        self._next_steps = list(steps)
        theme = self.transcript.theme
        self.transcript.chrome(theme.dim("  next"))
        for i, step in enumerate(steps, 1):
            self.transcript.chrome(theme.dim(f"    {i}  {step[:100]}"))

    def _tab_suggestion(self, buffer: str):
        """The next suggestion to put on an empty line, or None.

        None means "not mine": the line has something in it, or nothing
        is on offer, and Tab goes back to completing a ``/command``. That
        separation is the whole design — the two meanings of the key
        never meet, because one only applies to an empty line.

        Pressing on past the last suggestion gives the line back empty. A
        cycle with no way out is a prompt the user has to erase by hand.
        """
        steps = getattr(self, "_next_steps", None) or []
        if not isinstance(steps, list) or not steps:
            return None
        text = buffer or ""
        if not text:
            return steps[0]
        try:
            at = steps.index(text)
        except ValueError:
            return None          # the user typed; this is not the offer
        if at + 1 >= len(steps):
            return ""            # one past the last: the line comes back
        return steps[at + 1]

    def _expand_next_step(self, text: str) -> str:
        """A bare number at the prompt means the suggestion of that number.

        Only a bare number, and only while an offer stands: "2" is a
        message in its own right in a conversation about numbers, and it
        must stay one when nothing was offered.
        """
        steps = getattr(self, "_next_steps", None) or []
        raw = (text or "").strip()
        if not steps or not raw.isdigit():
            return text
        i = int(raw)
        if 1 <= i <= len(steps):
            chosen = steps[i - 1]
            self.transcript.chrome(
                self.transcript.theme.dim(f"  → {chosen[:100]}"))
            return chosen
        return text

    def _report_tasks(self) -> None:
        """Open work, when the agent left some and the user wants to see it."""
        if not self._show_tasks:
            return
        try:
            from . import task_ticker
            text = task_ticker.render_text(
                self.opts.cwd,
                session_id=str(getattr(self.engine, "session_id", "") or ""))
        except Exception:
            return
        if text and text != "(no tasks)":
            for line in text.splitlines():
                self.transcript.chrome(self.transcript.theme.dim(line))

    # -- the loop --------------------------------------------------------
    def _checkpoint_session(self) -> None:
        """Save the conversation after a finished turn. Never raises.

        The conversation used to be saved only on the way out --
        cmd_chat's finally -- so a process killed by SIGKILL (an
        expired job, a killed node) lost the whole session and ``-r``
        came back to an empty conversation. Saving at every turn end
        makes a kill cost at most the turn it interrupts: the
        mid-turn checkpoint (session_store.save_turn_checkpoint)
        still covers that one, and every FINISHED turn is already on
        disk. The saver is cli._save_session -- the same call the
        exit path makes, so the record and the resume are one
        contract, not two.

        Locked sessions are skipped silently: the writer lock exists
        so two processes never tear the file, and a turn-end save
        that cannot take it is no reason to end the turn.
        """
        try:
            engine = getattr(self, "engine", None)
            if engine is None:
                return
            from .cli import _save_session
            _save_session(engine, Path(self.opts.cwd))
        except Exception:
            pass

    def run(self, first_prompt: str = "") -> int:
        from . import repl_keys as rk

        self._install_sigint()
        self._load_history()
        self._install_completer()
        # Before the first prompt, so a session is addressable from the
        # moment it exists rather than from its first idle poll.
        self._announce_presence()
        self._arm_presence_heartbeat()
        if self.opts.banner:
            for line in self.opts.banner.splitlines():
                self.transcript.chrome(line)
        pending = first_prompt.strip()
        # A line queued during the preceding turn is printed at the moment
        # Enter consumes it.  Do not print it a second time when its turn
        # starts; ordinary idle-prompt input has not been shown yet.
        pending_already_shown = False
        try:
            while True:
                if not pending and self.queued:
                    pending = self.queued.pop(0)
                    pending_already_shown = True
                if not pending:
                    pending_already_shown = False
                    try:
                        from .repl_keys import raw_mode_supported
                        if raw_mode_supported(self._stdin):
                            # The framed box, only on a raw terminal —
                            # a pipe or a redirect keeps the readline
                            # path exactly as it was.
                            pending = self.read_boxed().strip()
                        else:
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
                if not pending_already_shown:
                    self._show_user_input(pending)
                pending_already_shown = False
                # A bare number stands for the suggestion of that number,
                # and only while one is on offer.
                pending = self._expand_next_step(pending)
                pending = self._handle_line(pending)
                if self._quit:
                    return 0
                if not pending:
                    continue
                pending_prompt = pending          # visible past the try
                try:
                    self.turn(pending)
                except KeyboardInterrupt:
                    self.transcript.chrome("")
                    return 130
                # The turn finished: the conversation goes on disk now,
                # not only when the process leaves through its finally.
                # A SIGKILL at the next prompt must cost nothing (see
                # _checkpoint_session).
                self._checkpoint_session()
                # A turn cut at the token ceiling without a tool call is
                # a stall, not an answer: on a model that thinks
                # invisibly the ceiling goes to thinking, the visible
                # answer ends at finish_reason "length", and the session
                # stands at the prompt half-way through the edit
                # (operator point d, 2026-09-22). Continue it, at most
                # twice; the third time the operator hears it.
                continuation = self._continue_after_length(pending_prompt)
                while continuation:
                    self._show_user_input(continuation, queued=True)
                    self._checkpoint_session()
                    # A continued turn is a turn: the third Ctrl+C
                    # raises KeyboardInterrupt out of turn() here the
                    # same way it does in the first one above, and the
                    # same hatch applies -- leave with 130 instead of
                    # raising out of run() past every handler.
                    try:
                        self.turn(continuation)
                    except KeyboardInterrupt:
                        self.transcript.chrome("")
                        return 130
                    self._checkpoint_session()
                    continuation = self._continue_after_length(continuation)
                pending = ""
        except rk.TerminalLeft as left:
            # Whatever the session was doing -- at the prompt, in a turn,
            # in a dialog -- the terminal is gone and so is the session.
            # Returning (rather than dying) is what saves it for a resume.
            # The deadline comes first: everything after it can hang.
            self._arm_leave_deadline(left.code)
            self._stop_for_leaving()
            self._end_what_it_started()
            return left.code
        finally:
            try:
                self._teardown_screen()
            except OSError:
                pass            # a line that hung up takes no erase sequence
            self._save_history()
            self._withdraw_presence()
            self._restore_sigint()
            self._restore_sigwinch()

    # -- plumbing --------------------------------------------------------
    def _install_sigint(self) -> None:
        try:
            self._prev_sigint = signal.signal(signal.SIGINT, self._on_sigint)
        except (ValueError, OSError):
            # Not the main thread, or no signal support. The loop still
            # works; Ctrl+C simply behaves as the default.
            self._prev_sigint = None
        self._install_leave_signals()
        self._install_sigwinch()

    def _install_leave_signals(self) -> None:
        """SIGHUP and SIGTERM leave the way a closed terminal does.

        Left at their defaults they end the process on the spot: nothing
        refuses the question in flight, and the conversation, which is
        saved on the way out, is lost. Caught, they unwind through the
        loop's cleanup, under the same deadline as any other leave.
        """
        self._prev_leave: dict = {}
        for name in ("SIGHUP", "SIGTERM"):
            sig = getattr(signal, name, None)
            if sig is None:
                continue
            try:
                self._prev_leave[sig] = signal.signal(
                    sig, self._on_leave_signal)
            except (ValueError, OSError):
                pass

    def _on_leave_signal(self, signum, _frame) -> None:
        """Same rule as the other handlers: never print, only unwind."""
        from . import repl_keys as rk
        raise rk.TerminalLeft(signal.Signals(signum).name.lower(),
                              128 + int(signum))

    def _restore_leave_signals(self) -> None:
        for sig, prev in (getattr(self, "_prev_leave", None) or {}).items():
            try:
                signal.signal(sig, prev)
            except (ValueError, OSError, TypeError):
                pass
        self._prev_leave = {}

    def _stop_for_leaving(self) -> None:
        """Nothing may go on waiting for a terminal that is gone.

        The turn is told to stop, and every question -- in flight or still
        to come from a worker that has not noticed yet -- is refused.
        """
        self._stop_engine()
        if self.broker is not None:
            try:
                self.broker.abort_all()
            except Exception:
                pass

    def _end_what_it_started(self) -> None:
        """Background shells and MCP servers go with the session.

        An atexit hook stops the shells at a normal exit, but it runs only
        after the interpreter has joined every worker thread -- a
        subagent's among them -- and the deadline's os._exit skips it. A
        session whose terminal is gone ends them itself.
        """
        try:
            from . import bash_jobs
            bash_jobs.get_registry().stop_running()
        except Exception:
            pass
        try:
            from . import mcp_client
            mcp_client.reset_registry()
        except Exception:
            pass

    def _arm_leave_deadline(self, code: int) -> None:
        """The process ends by *code* within the deadline, cleaned up or not.

        What follows the loop -- saving the session, the report, the
        lifeline ending what the session started -- normally takes well
        under a second. When it does not, the terminal is gone either way
        and the process does not get to outlive it by waiting.
        """
        timer = threading.Timer(_LEAVE_DEADLINE_S, os._exit, args=(code,))
        timer.daemon = True
        timer.start()

    def _install_sigwinch(self) -> None:
        """Notice the window changing size.

        SIGWINCH does not exist on every platform, and the signal cannot
        be installed off the main thread — both are ordinary, so both are
        silent. What is not silent is the consequence of ignoring it: the
        width is read once at startup and everything truncates to it
        forever, so widening the window gains nothing and narrowing it
        wraps every tool line.
        """
        sig = getattr(signal, "SIGWINCH", None)
        if sig is None:
            self._prev_sigwinch = None
            return
        try:
            self._prev_sigwinch = signal.signal(sig, self._on_sigwinch)
        except (ValueError, OSError):
            self._prev_sigwinch = None

    def _on_sigwinch(self, _signum, _frame) -> None:
        """Record the new size. Never paint from a signal handler.

        A handler runs between bytecodes of whatever the main thread was
        doing — including halfway through a write to the same terminal —
        so drawing here interleaves with the write it interrupted. Same
        rule the interrupt handler follows: set state, let the loop act.
        """
        self.transcript.refresh_width()
        self._width_dirty = True

    def _restore_sigwinch(self) -> None:
        sig = getattr(signal, "SIGWINCH", None)
        if sig is None or self._prev_sigwinch is None:
            return
        try:
            signal.signal(sig, self._prev_sigwinch)
        except (ValueError, OSError):
            pass
        self._prev_sigwinch = None

    def _restore_sigint(self) -> None:
        if self._prev_sigint is not None:
            try:
                signal.signal(signal.SIGINT, self._prev_sigint)
            except (ValueError, OSError):
                pass
        self._restore_leave_signals()

    @staticmethod
    def history_path() -> Path:
        from . import state_paths
        return state_paths.ensure_dir(Path.home() / ".delfin") / HISTORY_NAME

    def _install_completer(self) -> None:
        """Tab completes commands and paths INSIDE the workspace only.

        Offering paths outside it would advertise files the agent may not
        read, which is worse than no completion: it invites a request that
        is then refused.
        """
        try:
            import readline
        except ImportError:
            return
        from . import repl_commands as rc

        def _complete(text: str, state: int):
            try:
                if text.startswith("@"):
                    hits = ["@" + p for p in
                            rc.complete_path(text, self.opts.cwd)]
                elif text.startswith("/"):
                    hits = [n for n in sorted(rc.BUILTINS)
                            if n.startswith(text.lower())]
                else:
                    return None
                return hits[state] if state < len(hits) else None
            except Exception:
                return None

        try:
            readline.set_completer(_complete)
            readline.set_completer_delims(" \t\n")
            readline.parse_and_bind("tab: complete")
            # Pinned rather than inherited. Recent GNU readline defaults
            # this on and older builds — and libedit, which macOS links —
            # do not, so whether a pasted block arrived as one message
            # depended on which library the interpreter happened to find.
            # During a turn the key layer brackets the paste itself; this
            # is the same guarantee at the idle prompt.
            readline.parse_and_bind("set enable-bracketed-paste on")
        except Exception:
            pass

    def _recall_later(self, text: str) -> None:
        """Put a line typed DURING a turn into the recallable history.

        Deliberately not called `_remember`: that name is taken by the
        `#note` handler, which writes a durable MEMORY. Two things a user
        might call remembering, and conflating them would have turned
        every queued line into a stored memory.

        Only the idle prompt goes through readline, so everything typed
        while the agent worked — most of what gets typed in a busy
        session — was absent from the up-arrow afterwards. Retyping it is
        what a history exists to prevent.
        """
        if not text:
            return
        try:
            import readline
            readline.add_history(text)
        except Exception:
            pass                     # no readline: the queue still works

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
            # Where this session's own lines begin. Everything above it
            # came off disk and is already in the file; appending it again
            # would grow the file by its own length on every exit.
            self._hist_base = readline.get_current_history_length()
        except Exception:
            # libedit behaves differently, and a missing file is normal on
            # a first run. Neither is worth refusing to start over.
            pass

    def _save_history(self) -> None:
        try:
            import readline
            from . import state_paths
            path = self.history_path()
            # Append only what THIS session added. Writing the whole
            # in-memory buffer means two terminals open in two projects
            # each load the file at start and write their own copy at
            # exit, so whichever quits last erases the other's session.
            # libedit carries no append_history_file, and there the whole
            # write is the only option — losing the merge beats losing
            # the file.
            added = 0
            try:
                added = max(
                    0, readline.get_current_history_length() - self._hist_base)
            except Exception:
                added = 0
            appended = False
            if added and hasattr(readline, "append_history_file"):
                try:
                    readline.append_history_file(added, str(path))
                    appended = True
                except Exception:
                    appended = False
            if not appended:
                readline.write_history_file(str(path))
            # Prompts describe the user's project, so the file is theirs
            # alone -- the same 0600 every other state file here gets.
            state_paths.secure_file(path)
        except Exception:
            pass



    # -- the framed idle prompt ------------------------------------------
    #: How often the idle prompt looks for work that finished. The read
    #: loop ticks ten times a second; asking the job registry that often
    #: would be a poll, and this is a look.
    _WAKE_EVERY_S = 5.0

    def _steer_operator_mail(self) -> None:
        """Take waiting mail and steer it into the running turn.

        The mirror of _operator_messages for a session that is not at
        its prompt: same inbox, same rendering (the header says the
        message is from another session, never from the user), but the
        destination is the turn -- engine.steer, the same call the
        dashboard makes -- so the model sees it between rounds instead
        of at the end. On the idle wake's throttle, never per read;
        the rendered text is shown on the transcript because a steer
        is not echoed either. A backend that takes no steer returns
        False and the message is put back, unread, so the prompt path
        delivers it when the turn ends -- a fallback, not a loss.
        """
        import time as _time
        now = _time.monotonic()
        if now - getattr(self, "_mail_steer_last", 0.0) < self._WAKE_EVERY_S:
            return
        self._mail_steer_last = now
        key = self._presence_key()
        if not key:
            return
        from . import session_messages as _msgs
        messages = _msgs.take(key)
        if not messages:
            return
        text = "\n\n".join(_msgs.render(m) for m in messages)
        delivered = False
        try:
            delivered = bool(self.engine.steer(text))
        except Exception:
            delivered = False
        if not delivered:
            # The turn could not take it: the messages go back, newest
            # last, so the prompt path sees them in the order they
            # arrived.
            for message in messages:
                try:
                    _msgs.send(key, str(message.get("text") or ""),
                               from_key=str(message.get("from") or ""),
                               from_title=str(
                                   message.get("from_title") or ""))
                except Exception:
                    pass
            return
        for message in messages:
            sender = (message.get("from_title") or message.get("from")
                      or "an operator")
            body = " ".join(str(message.get("text") or "").split())
            self.transcript.chrome(
                self.transcript.theme.dim(f"✉ {sender} → into the running "
                                          f"turn: {body[:300]}"))

    #: At most this many self-continuations of a length-cut turn; the
    #: next one is the operator's to drive, and the operator is told.
    _MAX_LENGTH_CONTINUATIONS = 2

    def _continue_after_length(self, answered_prompt: str) -> str:
        """Whether the finished turn should be continued, and with what.

        A turn that ended at ``length`` WITHOUT a tool call spent its
        ceiling -- on a model that thinks invisibly, largely on
        thinking -- and stopped mid-answer. The terminal sends
        "continue" on its own, at most twice: the same message the
        person at the keyboard would type, rendered as a continuation
        and not as their words. A tool-carrying length end is left to
        the tool loop (the next round carries the work); a normal end
        needs no help. After the second continuation a third length
        end goes to the operator over session_message -- the channel
        that reaches a person not watching six panes -- and the turn
        is left at the prompt for them.
        """
        try:
            reason = str(getattr(self.engine, "last_turn_stop_reason",
                                 "") or "")
        except Exception:
            return ""
        if reason != "length":
            self.__dict__["_length_continuations"] = 0
            return ""
        if getattr(self, "_turn_used_tools", False):
            # The model meant to go on; the tool loop owns it.
            self.__dict__["_length_continuations"] = 0
            return ""
        done = int(getattr(self, "_length_continuations", 0) or 0)
        if done >= self._MAX_LENGTH_CONTINUATIONS:
            try:
                self.engine._notify_operator_of_stall(
                    "a turn keeps ending at the token ceiling without "
                    "finishing; the session is standing at the prompt "
                    "mid-answer after two self-continuations")
            except Exception:
                pass
            return ""
        self.__dict__["_length_continuations"] = done + 1
        self.transcript.chrome(self.transcript.theme.dim(
            "↻ the answer was cut at the token ceiling — continuing "
            f"({done + 1}/{self._MAX_LENGTH_CONTINUATIONS})"))
        return "continue"

    def _presence_key(self) -> str:
        """The address an operator reaches this session by.

        The name from ``-n`` when there is one, because that is what a
        person can type; the head of the session id otherwise, so a
        session is never unreachable for want of a name.
        """
        opts = getattr(self, "opts", None)
        name = str(getattr(opts, "session_name", "") or "").strip()
        if name:
            return name
        engine = getattr(self, "engine", None)
        return str(getattr(engine, "session_id", "") or "")[:8]

    def _arm_presence_heartbeat(self) -> None:
        """Let the broker refresh presence whenever it asks something.

        The idle poll is the only other heartbeat, and a session that
        goes from one approval into the next never reaches it.
        """
        broker = getattr(self, "broker", None)
        if broker is None:
            return
        try:
            broker.on_activity = self._announce_presence
        except Exception:
            pass

    def _announce_presence(self) -> None:
        """Say this session is open. Never raises.

        A terminal session announced nothing, so it was missing from
        ``open_sessions`` and the ``session_message`` tool told other
        sessions it does not exist. Repeated from the idle poll, which is
        the heartbeat: ``announce`` writes only when the record changed
        or the interval elapsed.
        """
        try:
            from . import session_presence
            session_presence.announce(
                self._presence_key(),
                session_id=str(getattr(self.engine, "session_id", "") or ""),
                title=str(getattr(self.opts, "session_name", "") or ""),
                workspace=str(getattr(self.opts, "cwd", "") or ""))
        except Exception:
            pass

    def _withdraw_presence(self) -> None:
        """Take the session out of the list. Never raises."""
        try:
            from . import session_presence
            session_presence.withdraw(self._presence_key())
        except Exception:
            pass

    def _operator_messages(self, typed: str) -> str:
        """What was left in this session's inbox, ready for the prompt.

        ``session_messages.send`` wrote a file that nobody read: ``take``
        was called in the dashboard and nowhere else. The only way into a
        terminal session was its terminal, and that way has a trap -- an
        approval dialog reads single keys, so text arriving while one is
        up does not queue, its first character ANSWERS the dialog.

        Delivered where the job wake-up is delivered: at the idle prompt,
        between reads. Never during a turn, and never while a dialog
        holds the reader, so the trap cannot arise here. The remaining
        guard is the dashboard's -- never while something is typed,
        because a draft belongs to whoever typed it, and a message that
        waits is not a message that is lost.
        """
        if typed.strip():
            return ""
        # One guard around the whole of it. A wake-up that throws takes
        # the prompt with it, and the prompt is the thing the user is
        # standing at -- the case that proves it drives a bare object
        # with no engine and no transcript at all.
        try:
            key = self._presence_key()
            if not key:
                return ""
            from . import session_messages as _msgs
            messages = _msgs.take(key)
            if not messages:
                return ""
            for message in messages:
                sender = (message.get("from_title") or message.get("from")
                          or "an operator")
                # The text too, not only the sender: the message becomes
                # the next prompt, and a prompt is not echoed. Somebody
                # watching the pane has to be able to read what arrived.
                body = " ".join(str(message.get("text") or "").split())
                self.transcript.chrome(
                    self.transcript.theme.dim(f"✉ {sender}: {body[:300]}"))
            return "\n\n".join(_msgs.render(m) for m in messages)
        except Exception:
            return ""

    def _wake_text(self, typed: str) -> str:
        """The message a job that finished hands an idle prompt, or "".

        Nothing in the terminal used to wake an agent. A background job
        that finished DURING a turn was reported at the end of it, but
        one that finished afterwards sat there: the drain is a pull, and
        the only thing that pulls is a turn, and the only thing that
        starts a turn is the user typing. A run started before a quiet
        night was simply never looked at again (reported 2026-09-18).

        Two guards, both the dashboard's:

        * never while something is typed — a wake-up must not take a
          half-written message away;
        * never when the user has turned it off
          (``agent.wake_on_job_end``).

        Never raises: a wake-up that throws would take the prompt with
        it, and the prompt is the thing the user is standing at.
        """
        import time as _time

        if typed.strip():
            return ""
        now = _time.monotonic()
        if now - getattr(self, "_wake_last_look", 0.0) < self._WAKE_EVERY_S:
            return ""
        self._wake_last_look = now
        # The operator comes first, and is not subject to
        # ``agent.wake_on_job_end``: switching off the report of a
        # finished job must not switch off a person.
        self._announce_presence()       # the poll doubles as the heartbeat
        mail = self._operator_messages(typed)
        if mail:
            return mail
        try:
            from . import job_wake
            if not job_wake.wake_enabled():
                return ""
            seen = self.__dict__.setdefault("_wake_seen", set())
            return job_wake.wake_prompt(job_wake.finished_shells(seen))
        except Exception:
            return ""

    #: How often the line under the prompt re-reads what is running. The
    #: box is redrawn on every keystroke; asking the registries that
    #: often would be a poll, and this is a glance.
    _BACKGROUND_STATUS_EVERY_S = 3.0

    def _background_status(self) -> str:
        """What is still out, in one line, for under the input area.

        The dashboard has had a panel for this; the terminal had `/bash`,
        which you had to think of asking — so a suite started twenty
        minutes ago was remembered or it was not. Same collector as the
        panel, so the two cannot drift.

        Cached between glances, and never raises: a decoration that can
        take the prompt with it is not a decoration.

        The last view is kept, not just its line: the arrow walk under
        the prompt walks THAT list, and re-rendering the line with a
        selection needs the same rows the unselected line came from —
        two collectors could return two different worlds between two
        keys.
        """
        import time as _time

        now = _time.monotonic()
        if now - getattr(self, "_bg_status_at", 0.0) < \
                self._BACKGROUND_STATUS_EVERY_S:
            return getattr(self, "_bg_status", "")
        self._bg_status_at = now
        try:
            from . import background_view as _bgv
            view = _bgv.collect(self.opts.cwd)
            self._bg_view = view
            self._bg_now = _time.time()
            self._bg_status = _bgv.status_line(view, now=self._bg_now)
        except Exception:
            self._bg_status = ""
        return self._bg_status

    def _bg_status_selected(self, selected: int | None) -> str:
        """The status line with one row marked, from the cached view.

        ``selected`` is an index into what the line shows. None (or a
        stale index, or no view) gives the unselected line unchanged —
        a mark that outlives its row must not move onto another one.

        Rendered with the ``now`` the cached view was collected at, so
        a walk between two glances cannot jump a duration forward: the
        rows are frozen with the list they came from.
        """
        view = getattr(self, "_bg_view", None)
        if selected is None or not view:
            return self._background_status()
        try:
            from . import background_view as _bgv
            return _bgv.status_line(view, selected=selected,
                                    now=getattr(self, "_bg_now", None))
        except Exception:
            return self._background_status()

    def _bg_id_at(self, index: int) -> str:
        """The id Enter would take for the row at *index*, or ""."""
        view = getattr(self, "_bg_view", None)
        if not view:
            return ""
        try:
            from . import background_view as _bgv
            return _bgv.selected_id(view, index,
                                    now=getattr(self, "_bg_now", None))
        except Exception:
            return ""

    def read_boxed(self) -> str:
        """One message through the framed box, on a raw terminal only.

        Everything the readline prompt did is preserved: history (the
        same readline list, via _BoxHistory), completion (the SAME
        completer _install_completer put on readline), multi-line
        blocks and backslash continuations (read_block's grammar,
        applied to what the box collected), bracketed paste (the key
        layer's), and the Ctrl+C-twice-to-leave ladder.

        A pipe, a redirect or a terminal without termios never reaches
        here: run() falls back to the readline path when
        raw_mode_supported() is false. The turn path (_pump,
        _repaint_bottom) is untouched — a turn that cannot be
        interrupted is worse than a plain prompt.

        Drawing discipline: the box is always the LAST thing on the
        screen, and it is redrawn row by row with \\r\\x1b[K, so no row
        ever wraps and no absolute cursor addressing is needed — the
        cursor is placed by writing the rows and then moving UP and
        across, relative to where the writes left it. That is why this
        works without asking the terminal where it is.
        """
        from . import repl_keys as rk
        from . import repl_box as rb

        # Bound to the live RawMode inside read_boxed's `with` — a
        # detached instance has no fd and read_ready would return ""
        # forever, which is exactly the hang a first run would show.
        read_raw: list = []          # [callable] once the context opens

        #: How long the idle prompt waits for a key before looking around.
        #: It is not keystroke latency: read_ready is a select, so a key
        #: returns at once and a signal (SIGWINCH) interrupts the wait —
        #: this only sets how often an idle prompt wakes to do nothing.
        #: At 0.1 s that was ten times a second where readline blocked.
        _IDLE_TICK_S = 0.5

        def _read_chunk() -> str:
            # RawMode.read_ready: one chunk or "" after a short timeout,
            # with paste-marker stitching already in the decoder.
            return read_raw[0](_IDLE_TICK_S) if read_raw else ""

        class _BlankDecoder:
            """A decoder-shaped blank, so _clear_box can size a box."""
            buffer = ""
            cursor = 0


        def _view(decoder, selected: int | None = None) -> rb.BoxView:
            if search.active:
                # The search prompt is the box while it lasts: the query
                # where typed text went, the readline-style label as the
                # hint row, so it is visible WHICH mode the keys drive.
                return rb.viewport(
                    rb.render_box(search.query, len(search.query),
                                    self.transcript.width,
                                    search.hint(),
                                    status=self._background_status()),
                    _BOX_MAX_CONTENT_ROWS)
            return rb.viewport(
                rb.render_box(decoder.buffer, decoder.cursor,
                                self.transcript.width,
                                _box_hint(self._posture_now()),
                                status=self._bg_status_selected(selected)),
                _BOX_MAX_CONTENT_ROWS)

        # Where the last _draw left the terminal: rows below the
        # cursor to the box's end, and what was last painted. Drawing
        # is relative — the cursor is walked DOWN to the box's end and
        # UP to its top border, never addressed absolutely — because
        # the box sits wherever the transcript ended, and asking the
        # terminal where that is needs a cursor-position reply this
        # loop has no reader for.
        # What is ON SCREEN, which is the only thing the walk may be
        # measured against: how many rows the box has and how many of
        # them sit below the cursor. The first version kept only
        # "below" and climbed by the NEW view's height, so every paint
        # mixed two geometries.
        _state = {"below": 0, "rows": 0, "painted": None}

        def _draw(decoder, selected: int | None = None, *,
                  force: bool = False) -> None:
            view = _view(decoder, selected)
            if not force and _state["painted"] == (view.rows,
                                                    view.cursor):
                return                      # nothing changed on screen
            _state["painted"] = (view.rows, view.cursor)
            on_screen = _state["rows"]
            out = []
            # Walk back over WHAT IS THERE. On the first paint nothing
            # is, and climbing would have painted the box over the last
            # rows of the transcript -- the banner at start-up, the end
            # of the answer after a turn.
            if on_screen:
                if _state["below"]:
                    out.append(f"\x1b[{_state['below']}B")   # to box end
                if on_screen > 1:
                    out.append(f"\x1b[{on_screen - 1}A")     # to its top
            out.append("\r")
            for i, row in enumerate(view.rows):
                if i:
                    out.append("\r\n")
                out.append("\x1b[K" + row)
            # A box that shrank (a line unwrapped) leaves rows of the
            # old one standing below it.
            spare = max(0, on_screen - len(view.rows))
            for _ in range(spare):
                out.append("\r\n\x1b[K")
            if spare:
                out.append(f"\x1b[{spare}A")
            # Rows written: the cursor sits after the LAST row. Put it
            # on the cursor's row and column: up by the rows below it,
            # then across (the content rows start at column zero;
            # clamped so a cursor on a border column cannot wrap).
            crow, ccol = view.cursor
            below = len(view.rows) - 1 - (crow + 1)     # minus top border
            _state["below"] = below
            _state["rows"] = len(view.rows)
            # A form that draws a frame needs its width stepped over;
            # the two-rule form and the narrow row do not, and adding
            # the offset there put the cursor two columns past the text.
            offset = 2 if getattr(view, "border", True) else 0
            col = min(ccol + offset, max(1, self.transcript.width - 1))
            if below:
                out.append(f"\x1b[{below}A")
            out.append(f"\r\x1b[{col}C")
            self.err.write("".join(out))
            self._flush_err()

        def _clear_box(_decoder=None) -> None:
            """Erase the box and leave the cursor one row BELOW where
            its end was — the transcript continues from there.

            Sized from what is ON SCREEN, never from the decoder: after
            a submit the decoder has already dropped its buffer, so a
            view built from it describes an empty box while a taller one
            is still painted. The rows the erase then missed stayed in
            the scrollback as a stray border.
            """
            n = _state["rows"]
            if not n:
                return
            out = []
            if _state["below"]:
                out.append(f"\x1b[{_state['below']}B")   # to box end
            if n > 1:
                out.append(f"\x1b[{n - 1}A")             # to top border
            out.append("\r")
            out.append("\x1b[K\r\n" * n)
            # A row that WRAPPED on the physical screen occupies more
            # rows than the view counted — the status line under the
            # input can exceed the width the box was sized for — so the
            # counted erase misses its tail, and the tail shows through
            # the next answer. Erase to the end of the screen as well:
            # it reaches below the box only, and the transcript above
            # is untouched because the walk never climbs past the top.
            out.append("\x1b[J")
            _state["below"] = 0
            _state["rows"] = 0
            _state["painted"] = None
            self.err.write("".join(out))
            self._flush_err()

        def _word_before(text: str) -> str:
            i = len(text)
            while i > 0 and text[i - 1] not in " \t\n":
                i -= 1
            return text[i:]

        class _ReverseSearch:
            """Ctrl+R's mode, over the SAME readline store the arrows use.

            The decoder reports the key; the history is here, in the
            loop, where _BoxHistory lives — a decoder that knew about
            history would be a key layer reading a data store. Typing
            extends the query (chars the box would otherwise insert),
            Ctrl+R again steps to an older match, Enter accepts — and
            submits, like readline — Esc restores the draft and leaves.
            """

            def __init__(self):
                self.active = False
                self.query = ""
                self.draft = ""
                self.pos = 0          # index into the history, newest first

            def start(self, buffer: str) -> None:
                self.active = True
                self.query = ""
                self.draft = buffer
                self.pos = 0

            def _entries(self) -> list[str]:
                try:
                    import readline
                    n = readline.get_current_history_length()
                    return [readline.get_history_item(i)
                            for i in range(n, 0, -1)]
                except Exception:
                    return []

            def _match(self) -> str | None:
                if not self.query:
                    return self.draft or None
                for offset, entry in enumerate(self._entries()):
                    if offset < self.pos:
                        continue
                    if self.query in entry:
                        self.pos = offset
                        return entry
                return None

            def older(self) -> None:
                """Ctrl+R again: keep the query, take the next match up."""
                self.pos += 1

            def hint(self) -> str:
                match = self._match()
                if match is None:
                    return f"(reverse-i-search)`{self.query}': no match"
                return f"(reverse-i-search)`{self.query}': {match}"

            def accept(self) -> str:
                """What Enter submits: the match, else the query."""
                match = self._match()
                return match if match is not None else self.query

            def leave(self) -> str:
                """What Esc restores: the draft from before the search."""
                self.active = False
                return self.draft

        search = _ReverseSearch()

        def _complete_word(text: str) -> str:
            """The box's Tab: the same completer readline's Tab runs.

            Delegates to the readline module's completer — the one
            _install_completer installed — so the two prompts cannot
            drift. A unique hit completes in place; several list on the
            transcript, exactly as readline does.
            """
            try:
                import readline
                fn = readline.get_completer()
            except Exception:
                return text
            if fn is None:
                return text
            word = _word_before(text)
            hits: list[str] = []
            for state in range(200):            # a runaway completer stops
                hit = fn(word, state)
                if hit is None:
                    break
                hits.append(hit)
            if not hits:
                return text
            if len(hits) == 1:
                return text[:len(text) - len(word)] + hits[0]
            # Readline fills the longest common prefix FIRST and lists
            # only beside it. Listing alone — what this did — left the
            # box at what was typed, so "/mod" + Tab showed /mode and
            # /model below and the user retyped the four letters the
            # hits already shared.
            prefix = hits[0]
            for hit in hits[1:]:
                while not hit.startswith(prefix):
                    prefix = prefix[:-1]
            if len(prefix) > len(word):
                text = text[:len(text) - len(word)] + prefix
            _clear_box(_BlankDecoder())
            for hit in hits[:20]:
                self.transcript.chrome(self.transcript.theme.dim(hit))
            return text

        def _collect() -> str:
            """read_block's grammar over what the box collects."""
            decoder = rk.KeyDecoder()
            history = _BoxHistory()
            # The arrow walk over the status line's rows. None means no
            # walk is on; an index marks the row Enter would take. It
            # lives only while the buffer is EMPTY — one letter types a
            # message, and a message goes to the model, not to a job.
            selected: int | None = None
            while True:
                _draw(decoder, selected)    # dedup: only real change paints
                chunk = _read_chunk()
                if not chunk:
                    if self._width_dirty:
                        self._width_dirty = False
                        self.transcript.refresh_width()
                        # rows differ under the new width, so the next
                        # _draw's own comparison repaints them
                    woke = self._wake_text(decoder.buffer)
                    if woke:
                        _clear_box()
                        return woke
                    continue
                submit_text = None
                for event in decoder.feed(chunk):
                    kind = event.kind
                    if search.active and kind == rk.EDIT:
                        # While searching, printable keys extend the
                        # query; the decoder already put them in its
                        # buffer, so the query follows the buffer.
                        search.query = decoder.buffer
                    elif search.active and kind in (rk.SUBMIT, rk.STEER):
                        submit_text = search.accept()
                        search.active = False
                    elif search.active and kind == rk.INTERRUPT:
                        restored = search.leave()
                        decoder.buffer = restored
                        decoder.cursor = len(restored)
                    elif search.active and kind == rk.SEARCH:
                        search.older()
                    elif kind == rk.SEARCH:
                        selected = None     # a search is not a walk
                        search.start(decoder.buffer)
                        decoder.buffer = ""
                        decoder.cursor = 0
                    elif kind == rk.SUBMIT:
                        # Enter with a row of the status line selected
                        # walks INTO that job — the list below the input
                        # was always a report you had to retype a handle
                        # for; now Enter is the way in. Only on an empty
                        # buffer: a message is a message.
                        if selected is not None and not event.text.strip():
                            job_id = self._bg_id_at(selected)
                            if job_id:
                                _clear_box(decoder)
                                return f"/bash {job_id}"
                        submit_text = event.text
                    elif kind == rk.HISTORY_PREV:
                        # On an empty line the arrows walk the status
                        # line's jobs first — the thing under the prompt
                        # is what the empty prompt is about. Only when
                        # there is nothing to walk (or the walk is off
                        # the end) does the readline history take over,
                        # as it always did.
                        if not decoder.buffer and selected is None \
                                and self._bg_id_at(0):
                            selected = 0
                        elif not decoder.buffer and selected is not None \
                                and self._bg_id_at(selected + 1):
                            selected += 1
                        else:
                            older = history.up(decoder.buffer)
                            if older is not None:
                                decoder.buffer = older
                                decoder.cursor = len(older)
                                selected = None
                    elif kind == rk.HISTORY_NEXT:
                        if not decoder.buffer and selected is not None:
                            if selected > 0 and self._bg_id_at(selected - 1):
                                selected -= 1
                            else:
                                selected = None    # off the top: back to typing
                        else:
                            newer = history.down()
                            if newer is not None:
                                decoder.buffer = newer
                                decoder.cursor = len(newer)
                            selected = None
                    elif kind in (rk.STEER, rk.EXPAND, rk.TASKS):
                        # Turn-time keys. There is no turn here, and the
                        # decoder has already taken the line into the
                        # event — so without this branch Ctrl+G emptied
                        # the box and nothing happened, which is the one
                        # outcome a key must never have.
                        text = getattr(event, "text", "") or ""
                        if text:
                            decoder.buffer = text
                            decoder.cursor = len(text)
                    elif kind == rk.COMPLETE:
                        # On an empty line Tab walks the offer; with
                        # something typed it completes a /command as it
                        # always has. It FILLS — pressing return stays the
                        # user's, because a key that sent the turn would
                        # make an offer into a trap.
                        picked = self._tab_suggestion(decoder.buffer)
                        if picked is not None:
                            decoder.buffer = picked
                            decoder.cursor = len(picked)
                        else:
                            done = _complete_word(decoder.buffer)
                            if done != decoder.buffer:
                                decoder.buffer = done
                                decoder.cursor = len(done)
                    elif kind == rk.CYCLE_MODE:
                        _clear_box(decoder)
                        self._cycle_mode()
                    elif kind == rk.REDRAW:
                        self.transcript.refresh_width()
                    elif kind == rk.EOF:
                        raise EOFError
                    elif kind == rk.INTERRUPT:
                        # Esc at the idle prompt clears the line, like
                        # readline's Ctrl+U — not the session. It also
                        # drops a selection: the escape out of a walk
                        # leaves the walk.
                        decoder.buffer = ""
                        decoder.cursor = 0
                        selected = None

                if submit_text is None:
                    # The walk lives only on an empty buffer: one typed
                    # letter is a message, and a message goes to the
                    # model, not to a job. Plain keys never come through
                    # here as events — the decoder eats them straight
                    # into its buffer — so the check is on the buffer.
                    if decoder.buffer and selected is not None:
                        selected = None
                    if self._width_dirty:
                        self._width_dirty = False
                        self.transcript.refresh_width()
                    continue
                _clear_box(decoder)
                return submit_text

        with rk.RawMode(self._stdin) as raw:
            read_raw.append(raw.read_ready)
            text = _collect()
            if text.strip() == _BLOCK_FENCE:
                block: list[str] = []
                decoder = rk.KeyDecoder()
                while True:
                    _draw(decoder)
                    chunk = _read_chunk()
                    if not chunk:
                        continue
                    done = False
                    for event in decoder.feed(chunk):
                        if event.kind == rk.SUBMIT:
                            _clear_box(decoder)
                            if event.text.strip() == _BLOCK_FENCE:
                                done = True
                            else:
                                block.append(event.text)
                        elif event.kind == rk.EDIT:
                            pass
                        elif event.kind == rk.EOF:
                            raise EOFError
                    if done:
                        return "\n".join(block)
            # Backslash continuation, read_block's other half.
            parts = [text]
            while parts[-1].endswith("\\"):
                parts[-1] = parts[-1][:-1]
                parts.append(_collect())
            result = "".join(parts) if len(parts) > 1 else parts[0]
            _BoxHistory().add(result)
            return result

# ---------------------------------------------------------------------------
# The framed idle prompt
# ---------------------------------------------------------------------------

_BOX_HINT = "esc interrupt · shift+tab approval mode · /help"


def _box_hint(posture: str) -> str:
    """The hint row, with the approval posture in front of it.

    The posture was in the banner, and the banner scrolls away after the
    first turn. So the one fact that decides what the next command will
    do — whether it asks, or edits, or asks nothing at all — was off
    screen for the rest of the session. The hint row was always drawn, so
    saying it there costs no row at all.
    """
    name = str(posture or "").strip()
    return f"{name} · {_BOX_HINT}" if name else _BOX_HINT
#: Content rows the box shows at most. Anything taller is a viewport
#: with … markers (repl_box.viewport): a 300-line paste into a 24-row
#: window would otherwise push the frame off the top of the screen.
_BOX_MAX_CONTENT_ROWS = 6


class _BoxHistory:
    """History for the box, backed by the SAME store readline uses.

    The readline module's in-memory list is the single source of truth:
    the readline path (pipes, non-terminals) and the box path (raw
    terminal) must recall the same lines, and a second file would drift
    within one session. Only the box path can WRITE through here, so
    the readline path keeps its own save logic untouched.
    """

    def __init__(self):
        self._loaded = False
        self._index = -1          # -1: browsing; >= 0: walking the list
        self._draft = ""

    def _ensure(self) -> bool:
        """readline present; the file is ALREADY loaded.

        ``TerminalAgent.run`` calls ``_load_history`` before the first
        prompt on every path, so re-reading the file here would load
        every line twice and double the file on save. All this class
        needs is the module; the in-memory list is shared state.
        """
        if self._loaded:
            return True
        try:
            import readline  # noqa: F401
        except ImportError:
            return False
        self._loaded = True
        return True

    def up(self, draft: str) -> str | None:
        """One step back; None when there is nothing older."""
        if not self._ensure():
            return None
        try:
            import readline
            n = readline.get_current_history_length()
        except Exception:
            return None
        if n == 0:
            return None
        if self._index < 0:
            self._index = n - 1
            self._draft = draft
        elif self._index > 0:
            self._index -= 1
        else:
            return None
        try:
            import readline
            return readline.get_history_item(self._index + 1)
        except Exception:
            return None

    def down(self) -> str | None:
        """One step forward; the draft when the walk runs off the end."""
        if self._index < 0:
            return None
        self._index += 1
        try:
            import readline
            n = readline.get_current_history_length()
        except Exception:
            n = 0
        if self._index >= n:
            self._index = -1
            return self._draft
        try:
            import readline
            return readline.get_history_item(self._index + 1)
        except Exception:
            return None

    def add(self, text: str) -> None:
        """Same store, same order, as the readline path's add_history."""
        if not text:
            return
        try:
            import readline
            readline.add_history(text)
        except Exception:
            pass
        self._index = -1
        self._draft = ""
