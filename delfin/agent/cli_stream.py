"""What a CLI agent shows while it works.

The dashboard tab renders the turn live: each tool call with its
argument summary, elapsed time, a spinner, and a closing line with
tokens and duration. ``delfin-agent`` prints far less — the answer
text and one summary. This module is the missing rendering, and only
the rendering: pure functions over plain data.

Event shapes are the ones ``delfin/agent/cli.py::_run_once`` already
emits for ``--output-format stream-json`` (cli.py:288,297):

* ``{"type": "text", "text": str}``
* ``{"type": "tool_use", "name": str, "input": dict, "elapsed_s": float?}``

plus two this module adds for its own renderer:

* ``{"type": "tool_result", "name": str, "elapsed_s": float}`` —
  a call that has finished; rendered as its closing line.
* ``{"type": "tick"}`` — a spinner frame while a call is outstanding.

Non-tty behaviour: no spinner, no ANSI, one plain line per event —
piping the output into a file or another command stays readable.

No engine import, no curses, no threads. The only terminal control
is ANSI dim/bold, and only when ``is_tty``.
"""

from __future__ import annotations

import sys
import time
from typing import Any, Callable, Iterable, Iterator

# Spinner frames — braille, one frame per tick while a call is open.
SPINNER_FRAMES = ("⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏")

# How much of each argument value to show before ellipsis.
_ARG_MAX = 40

_ANSI_DIM = "\x1b[2m"
_ANSI_BOLD = "\x1b[1m"
_ANSI_RESET = "\x1b[0m"


def _dim(text: str, is_tty: bool) -> str:
    """Dim for a terminal, verbatim for a pipe or file."""
    if not is_tty:
        return text
    return f"{_ANSI_DIM}{text}{_ANSI_RESET}"


def _summarize_args(args: dict[str, Any]) -> str:
    """``path=delfin/agent/cli.py limit=40`` — keys in order, values cut."""
    parts: list[str] = []
    for key, value in args.items():
        raw = str(value).replace("\n", " ")
        if len(raw) > _ARG_MAX:
            raw = raw[:_ARG_MAX - 1] + "…"
        parts.append(f"{key}={raw}")
    return " ".join(parts)


def render_tool_call(name: str, args: dict[str, Any] | None,
                     state: dict[str, Any]) -> str:
    """One line for one tool call.

    ``state`` is plain data: ``{"is_tty": bool, "elapsed_s": float?,
    "spinner": int?}``. In a tty the line leads with a spinner frame
    (only when one is given, i.e. the call is still open) and the
    elapsed time is bold; in a pipe it is ``· name k=v …`` with the
    duration in plain parentheses, so grepping still works.
    """
    is_tty = bool(state.get("is_tty"))
    arg_summary = _summarize_args(args or {})
    elapsed = state.get("elapsed_s")
    spinner_index = state.get("spinner")

    head = ""
    if is_tty and spinner_index is not None:
        frame = SPINNER_FRAMES[int(spinner_index) % len(SPINNER_FRAMES)]
        head = f"{frame} "
    elif not is_tty:
        head = "· "

    body = name
    if arg_summary:
        body = f"{name} {arg_summary}"
    if is_tty:
        body = f"{_ANSI_BOLD}{body}{_ANSI_RESET}"

    if elapsed is None:
        return f"{head}{body}"

    return f"{head}{body} {_dim(f'({elapsed:.1f}s)', True)}" if is_tty else f"{head}{body} ({elapsed:.1f}s)"


def render_turn_footer(stats: dict[str, Any]) -> str:
    """The closing line of a turn: tokens, tool calls, duration.

    ``stats``: ``{"input_tokens": int?, "output_tokens": int?,
    "tool_calls": int?, "duration_s": float?, "is_tty": bool?}``.
    Missing fields are omitted rather than guessed at zero.
    """
    is_tty = bool(stats.get("is_tty"))
    parts: list[str] = []
    if stats.get("input_tokens") is not None:
        parts.append(f"in={stats['input_tokens']}")
    if stats.get("output_tokens") is not None:
        parts.append(f"out={stats['output_tokens']}")
    if stats.get("tool_calls") is not None:
        parts.append(f"tools={stats['tool_calls']}")
    if stats.get("duration_s") is not None:
        parts.append(f"{stats['duration_s']:.1f}s")
    if not parts:
        parts.append("turn ended")
    line = "→ " + " · ".join(parts)
    return _dim(line, is_tty) if not is_tty else f"{_ANSI_BOLD}{line}{_ANSI_RESET}"


class StreamRenderer:
    """Turn an event stream into rendered lines.

    Feed it any iterable of event dicts (a generator from the engine
    adapter, a list in a test) and iterate the renderer to get one
    string per line to print. ``is_tty`` defaults to ``sys.stdout``
    and can be forced for tests and pipes; ``now`` is a callable
    returning seconds, so tests drive elapsed time without sleeping.
    """

    def __init__(self, events: Iterable[dict[str, Any]], *,
                 is_tty: bool | None = None,
                 now: Callable[[], float] = time.monotonic) -> None:
        self._events = events
        self._is_tty = (sys.stdout.isatty() if is_tty is None
                        else bool(is_tty))
        self._now = now
        self._tool_calls = 0

    @property
    def tool_calls(self) -> int:
        """How many ``tool_use`` events were rendered — for the footer."""
        return self._tool_calls

    def __iter__(self) -> Iterator[str]:
        started: dict[str, float] = {}
        tick = 0
        for event in self._events:
            kind = event.get("type")
            if kind == "text":
                text = str(event.get("text", ""))
                if text.strip():
                    yield text
            elif kind == "tool_use":
                self._tool_calls += 1
                name = str(event.get("name", ""))
                started[name] = self._now()
                yield render_tool_call(
                    name, event.get("input"),
                    {"is_tty": self._is_tty,
                     "spinner": tick if self._is_tty else None,
                     "elapsed_s": event.get("elapsed_s")})
            elif kind == "tool_result":
                name = str(event.get("name", ""))
                start = started.pop(name, None)
                elapsed = event.get("elapsed_s")
                if elapsed is None and start is not None:
                    elapsed = max(0.0, self._now() - start)
                yield render_tool_call(
                    name, event.get("input"),
                    {"is_tty": self._is_tty,
                     "elapsed_s": elapsed})
            elif kind == "tick":
                tick += 1
                if self._is_tty and started:
                    # A spinner frame only while a call is outstanding —
                    # a tick with nothing open renders nothing.
                    name = next(iter(reversed(list(started))))
                    yield render_tool_call(
                        name, None,
                        {"is_tty": self._is_tty, "spinner": tick,
                         "elapsed_s": self._now() - started[name]})
            elif kind == "turn_end":
                yield render_turn_footer(
                    {"input_tokens": event.get("input_tokens"),
                     "output_tokens": event.get("output_tokens"),
                     "tool_calls": event.get("tool_calls", self._tool_calls),
                     "duration_s": event.get("duration_s"),
                     "is_tty": self._is_tty})
