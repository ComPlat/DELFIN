"""What a CLI agent shows while it works.

The dashboard tab renders the turn live: each tool call with its
argument summary, elapsed time, a spinner, and a closing line with
tokens and duration. ``delfin-agent`` prints far less — the answer
text and one summary. This module is the missing rendering, and only
the rendering: pure functions over plain data.

Event shapes are the ones ``delfin/agent/cli.py::_run_once`` already
emits for ``--output-format stream-json`` (cli.py:288,297):

* ``{"type": "text", "text": str}``
* ``{"type": "notice", "text": str}`` — harness speech (retry banners,
  a stop, a cost ceiling); shown, never counted as the answer.
* ``{"type": "wait", "text": str}`` — the "waiting for the model … Ns"
  tick while the endpoint has not answered. Progress, not history: each
  one replaces the one before it, so on a tty it is rendered as a
  rewrite of the current row (no trailing newline, so the CALLER must
  print it with ``end=""``), and in a pipe only the first is written.
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

import datetime
import sys
import time
from typing import Any, Callable, Iterable, Iterator

from delfin.agent.turn_summary import format_turn_summary

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
    """The closing line of a turn: the turn summary, plus tokens.

    The summary line comes from :mod:`delfin.agent.turn_summary`
    ("Surfaced after 31s · done 11:11 · 1 shell still running"); the
    token figures stay as a dimmed suffix when the event carries them.
    ``stats``: ``{"input_tokens", "output_tokens", "tool_calls",
    "duration_s", "ended_at", "outcome", "running_jobs", "is_tty"}``;
    missing fields are omitted rather than guessed at zero.
    """
    is_tty = bool(stats.get("is_tty"))
    ended_at = stats.get("ended_at")
    if isinstance(ended_at, str):
        # The turn_end event carries the end time as an ISO string;
        # the summary wants a datetime. Naive parse: the event is in
        # local time already, and a malformed string renders no time
        # rather than breaking the footer.
        try:
            ended_at = datetime.datetime.fromisoformat(ended_at)
        except ValueError:
            ended_at = None
    summary = format_turn_summary(
        duration_s=stats.get("duration_s"),
        ended_at=ended_at,
        outcome=stats.get("outcome", "ok"),
        running_jobs=stats.get("running_jobs"),
    )
    parts: list[str] = []
    if stats.get("input_tokens") is not None:
        parts.append(f"in={stats['input_tokens']}")
    if stats.get("output_tokens") is not None:
        parts.append(f"out={stats['output_tokens']}")
    if stats.get("tool_calls") is not None and stats["tool_calls"] > 0:
        parts.append(f"tools={stats['tool_calls']}")
    if parts:
        suffix = " · ".join(parts)
        return summary + f" · {_dim(suffix, True)}"
    return summary


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
        # One failure anywhere in the turn marks the outcome the footer
        # reports: a turn that errored never renders a success verb.
        self._saw_error = False

    @property
    def tool_calls(self) -> int:
        """How many ``tool_use`` events were rendered — for the footer."""
        return self._tool_calls

    def __iter__(self) -> Iterator[str]:
        started: dict[str, float] = {}
        tick = 0
        waited = False   # a pipe gets the first wait tick only
        for event in self._events:
            kind = event.get("type")
            if kind == "error":
                # An error event is a failed turn: the footer must not
                # call it "Surfaced" later. Rendered as a notice line.
                self._saw_error = True
                text = str(event.get("text", ""))
                if text.strip():
                    yield f"⚠ {text}"
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
            elif kind == "notice":
                # Harness speech (retry banners, a stop, a cost ceiling)
                # — shown, never mistaken for answer text. Each one is a
                # separate event in the transcript; see "wait" below for
                # the one that is not.
                text = str(event.get("text", ""))
                if text.strip():
                    yield text
            elif kind == "wait":
                # The first-byte tick, "waiting for the model … Ns". It
                # is progress and not history: every tick replaces the
                # one before it, so it must occupy ONE row however long
                # the endpoint takes. The producer emits one per tenth of
                # the request deadline capped at ten seconds, so the
                # default 600 s deadline is sixty ticks, and appending
                # them filled the screen with one sentence.
                #
                # On a tty: return to column one and erase to the end of
                # the line, and end without a newline, so the next tick
                # overwrites this one. The caller prints what it is
                # given; a trailing newline here would leave every tick
                # standing.
                #
                # In a pipe there is no cursor to move, so the first tick
                # is written (a log has to say the endpoint went quiet)
                # and the rest are dropped. Dropping the later ones loses
                # the elapsed figure; keeping them is the flood this
                # change exists to remove, and the round's real duration
                # is in the turn footer either way.
                text = str(event.get("text", ""))
                if not text.strip():
                    continue
                if self._is_tty:
                    yield f"\r\x1b[K{text}"
                elif not waited:
                    waited = True
                    yield text
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
                     "ended_at": event.get("ended_at"),
                     "running_jobs": event.get("running_jobs"),
                     "outcome": ("stopped" if self._saw_error
                                 else event.get("outcome", "ok")),
                     "is_tty": self._is_tty})
