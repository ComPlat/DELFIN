"""What a turn shows, as protocol events instead of terminal lines.

``TerminalAgent`` writes everything through its ``transcript``: the model's
answer, tool calls and results, harness notices, the line that says a turn
is waiting. ``ProtocolTranscript`` takes that place in the worker and turns
each of those into events with stable identities, so the TUI can draw a
tool call as one row that is completed later, and stream an answer into
the block it belongs to.

Grouping rules:

* answer text arrives in deltas; a run of them is one ``text`` part, closed
  by the next tool call, notice or the end of the turn (``text.ended``
  carries the whole text, which is what the durable log keeps);
* thinking is grouped the same way into ``reasoning`` parts;
* the engine reports tool results by tool NAME, not by call, so results are
  matched to calls first-in-first-out per name -- the rule ``run_turn``
  already uses -- until the engine carries call ids itself.
"""

from __future__ import annotations

import contextlib
import json
import threading
from collections import defaultdict, deque
from typing import Any, Callable

from .. import repl_render as rr

Emit = Callable[..., None]


def _jsonable(value: Any) -> Any:
    """A value that survives json.dumps; anything exotic becomes text."""
    try:
        return json.loads(json.dumps(value, default=str))
    except (TypeError, ValueError):
        return str(value)


class ProtocolTranscript:
    """The ``Transcript`` interface ``TerminalAgent`` needs, emitting events."""

    #: Notices are formatted by the REPL code to a terminal width. The TUI
    #: re-wraps them anyway; this only bounds the lines the REPL pads.
    width = 100

    def __init__(self, emit: Emit):
        self._emit = emit
        self._lock = threading.RLock()
        self.theme = rr.Theme(enabled=False)
        self.show_tools = True
        self.show_thinking = True
        self.turn_id = ""
        self._seq = 0
        self._text: tuple[str, list[str]] | None = None       # (part id, chunks)
        self._reasoning: tuple[str, list[str]] | None = None
        self._calls: dict[str, deque[str]] = defaultdict(deque)
        self.tool_calls = 0
        # Lines said inside block(), and the thread whose they are.
        self._block: list[str] | None = None
        self._block_owner = 0

    # -- turn boundaries (driven by the worker's agent) -----------------------
    def begin_turn(self, turn_id: str) -> None:
        with self._lock:
            self._close_parts()
            self.turn_id = turn_id
            self._calls.clear()
            self.tool_calls = 0

    def finish(self) -> None:
        with self._lock:
            self._close_parts()

    # -- the Transcript surface ------------------------------------------------
    def refresh_width(self) -> int:
        return self.width

    @property
    def answer_open(self) -> bool:
        return self._text is not None

    screen_answer_open = answer_open
    screen_answer_column = 0

    def chrome(self, line: str) -> None:
        """Harness speech: command output, reports, notes about the turn."""
        text = rr.strip_control(line or "").rstrip()
        with self._lock:
            if self._block is not None and \
                    threading.get_ident() == self._block_owner:
                self._block.append(text)
                return
            if not text:
                return
            self._close_parts()
            self._emit("notice", text=text)

    @contextlib.contextmanager
    def block(self):
        """What this thread says until the end is ONE notice.

        The REPL prints a command's output line by line (``/help``, ``!ls``),
        blank lines included; as a notice per line, the blank lines that
        separate its paragraphs were lost and a window could not tell where
        one output ended. Collected here, it arrives as one text with its
        lines as they were.
        """
        with self._lock:
            outer = self._block is not None
            if not outer:
                self._block = []
                self._block_owner = threading.get_ident()
        try:
            yield
        finally:
            if not outer:
                with self._lock:
                    lines, self._block = self._block or [], None
                    text = "\n".join(lines).strip("\n")
                    if text.strip():
                        self._close_parts()
                        self._emit("notice", text=text)

    def transient(self, line: str) -> None:
        text = rr.strip_control(line or "").strip()
        if text:
            self._emit("wait", text=text)

    def answer(self, delta: str) -> None:
        self._text_delta(delta)

    def render(self, item) -> None:
        kind = getattr(item, "kind", "")
        if kind == "call":
            # Something to say in stream order (the worker's post): an
            # approval asked or settled, a state change. Run outside this
            # lock: it may take the broker's, which an asking thread holds
            # while it takes this one.
            fn = (getattr(item, "data", None) or {}).get("fn")
            if callable(fn):
                fn()
            return
        with self._lock:
            if kind == "text":
                self._text_delta(item.text)
            elif kind == "thinking":
                self._reasoning_delta(item.text)
            elif kind == "tool_use":
                self._close_parts()
                self._tool_started(item)
            elif kind == "tool_result":
                self._close_parts()
                self._tool_ended(item)
            elif kind == "denied":
                self._close_parts()
                calls = self._calls.get(item.name)
                call_id = calls.popleft() if calls else ""
                self._emit("tool.denied", turnId=self.turn_id,
                           callId=call_id or None, name=item.name)
            elif kind == "notice":
                self.chrome(item.text)
            elif kind == "wait":
                self.transient(item.text)
            elif kind == "error":
                self._close_parts()
                self._emit("error", text=rr.strip_control(item.text or ""))

    # -- parts -----------------------------------------------------------------
    def _next_id(self, letter: str) -> str:
        self._seq += 1
        return f"{self.turn_id or 'idle'}.{letter}{self._seq}"

    def _text_delta(self, delta: str) -> None:
        if not delta:
            return
        with self._lock:
            if self._reasoning is not None:
                self._close_reasoning()
            if self._text is None:
                self._text = (self._next_id("t"), [])
                self._emit("text.started", turnId=self.turn_id,
                           partId=self._text[0])
            self._text[1].append(delta)
            self._emit("text.delta", partId=self._text[0], delta=delta)

    def _reasoning_delta(self, delta: str) -> None:
        if not delta:
            return
        with self._lock:
            if self._text is not None:
                self._close_text()
            if self._reasoning is None:
                self._reasoning = (self._next_id("r"), [])
                self._emit("reasoning.started", turnId=self.turn_id,
                           partId=self._reasoning[0])
            self._reasoning[1].append(delta)
            self._emit("reasoning.delta", partId=self._reasoning[0],
                       delta=delta)

    def _close_text(self) -> None:
        if self._text is not None:
            part_id, chunks = self._text
            self._text = None
            self._emit("text.ended", turnId=self.turn_id, partId=part_id,
                       text="".join(chunks))

    def _close_reasoning(self) -> None:
        if self._reasoning is not None:
            part_id, chunks = self._reasoning
            self._reasoning = None
            self._emit("reasoning.ended", turnId=self.turn_id,
                       partId=part_id, text="".join(chunks))

    def _close_parts(self) -> None:
        self._close_text()
        self._close_reasoning()

    def close_parts(self) -> None:
        """End the open text or thinking part: what comes next (an approval
        asked mid-answer) is drawn below it, as a tool call is."""
        with self._lock:
            self._close_parts()

    def _tool_started(self, item) -> None:
        call_id = self._next_id("c")
        self._calls[item.name].append(call_id)
        self.tool_calls += 1
        args = _jsonable(rr._as_dict(item.text))
        if isinstance(args, dict) and isinstance(args.get("command"), str):
            # The REPL's headline masks a credential in a command
            # (repl_render._headline_value -> cli_approve.redact); the
            # secret redactor does not know an `export API_TOKEN=...` it
            # was never told. Masked here, before any window or the log.
            from .. import cli_approve
            args = {**args, "command": cli_approve.redact(args["command"])}
        self._emit("tool.started", turnId=self.turn_id, callId=call_id,
                   name=item.name, input=args)

    def _tool_ended(self, item) -> None:
        calls = self._calls.get(item.name)
        call_id = calls.popleft() if calls else self._next_id("c")
        meta = _jsonable(dict(item.data or {}))
        self._emit("tool.ended", turnId=self.turn_id, callId=call_id,
                   name=item.name, output=rr.strip_control(item.text or ""),
                   meta=meta)
