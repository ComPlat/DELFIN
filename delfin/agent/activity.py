"""What a running turn is doing right now, as one state both surfaces read.

Input: calls from the engine's stream loop (``begin_turn``, ``output``,
``tool_start``, ``tool_end``, ``end_turn``) and from whoever asks the
person (``awaiting``). Output: ``snapshot()`` -- a frozen ``State`` --
and two pure renderers over it, ``line`` (the indicator text) and
``slow_hint`` (one sentence when a phase has run past its threshold, or
"").

Phases::

    model_wait  a request is out and no token of it has arrived
    streaming   the model is producing (text, reasoning or a tool call)
    tool        a tool call is executing
    person      the turn is blocked on an approval or a question
    done        the turn ended (``outcome``: done | stopped | failed)
    idle        no turn has run on this engine

Every phase carries its own start time, so the clock a surface shows is
the time spent in THIS phase -- "waiting for glm … 45 s" counts from the
request, not from the turn. ``awaiting`` nests over a tool call (an
approval happens inside the tool that asked) and puts the interrupted
phase back, start time included, when the person has answered.

The surfaces repaint from a timer of their own and read the snapshot; no
update depends on the model producing anything. That is the point: the
indicators that existed before were rewritten by model events (a token,
a tool call, a first-byte tick every tenth of the deadline), so a
connected stream with no token, a running test suite, or the request
after a tool result showed a frozen label or none.

Rejected: inferring the phase from the rendered output on each surface
(the dashboard read "Thinking"/"Writing" out of its own label text to
pick an icon). Two surfaces guessing from two vocabularies is how the
terminal and the dashboard came to say different things about the same
turn. Rejected: a wrapper around the gate's confirm callback to detect
approvals -- the gate reads ``__self__`` off that callback, so the
asking surface marks the wait itself.

Cost: a lock and a handful of attribute writes per stream event; the
renderers are string formatting. No thread, no I/O.
"""

from __future__ import annotations

import os
import threading
import time
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Any, Callable, Iterator

IDLE = "idle"
MODEL_WAIT = "model_wait"
STREAMING = "streaming"
TOOL = "tool"
PERSON = "person"
DONE = "done"

APPROVAL = "approval"
ANSWER = "answer"

# Seconds before a phase earns its slow hint. The first-token figure is
# below the measured KIT GLM cold start (77-136 s) on purpose: the hint
# says what is slow and how to stop it, it does not say anything failed.
FIRST_TOKEN_SLOW_S = 60.0
SILENT_STREAM_SLOW_S = 60.0
TOOL_SLOW_S = 30.0
# Tools whose normal run is long: a shell command, a test suite, a
# sub-agent. Their hint waits for what such a call usually takes.
_LONG_TOOLS = {
    "bash": 120.0, "bash_background": 120.0, "run_tests": 120.0,
    "apply_patch": 60.0, "bash_status": 120.0, "bash_output": 60.0,
    "subagent": 600.0, "orchestrate": 600.0, "Agent": 600.0, "Task": 600.0,
    "web_fetch": 60.0, "web_search": 60.0, "Bash": 120.0,
}

_ARG_KEYS = ("command", "file_path", "path", "pattern", "query", "url",
             "description", "notebook_path", "name")
_ARG_MAX = 48


def _env_s(name: str, default: float) -> float:
    try:
        value = float(os.environ.get(name, "") or 0)
    except (TypeError, ValueError):
        value = 0.0
    return value if value > 0 else default


def bare_tool_name(name: str) -> str:
    """``mcp__kit-coding__bash`` -> ``bash``; other names unchanged."""
    name = str(name or "")
    if name.startswith("mcp__"):
        return name.rsplit("__", 1)[-1]
    return name


def short_arg(tool_input: Any) -> str:
    """The one argument a person recognises a call by, cut to one row.

    Accepts the dict or its JSON text (the engine's callbacks carry the
    text). The first of command/path/pattern/query/url wins; a path is
    reduced to its last two segments.
    """
    args = tool_input
    if isinstance(args, str):
        try:
            import json
            args = json.loads(args)
        except (TypeError, ValueError):
            args = {}
    if not isinstance(args, dict):
        return ""
    for key in _ARG_KEYS:
        value = args.get(key)
        if value is None or value == "":
            continue
        text = " ".join(str(value).split())
        if key in ("file_path", "path", "notebook_path") and "/" in text:
            text = "/".join(text.rstrip("/").split("/")[-2:])
        if len(text) > _ARG_MAX:
            text = text[:_ARG_MAX - 1] + "…"
        return text
    return ""


def fmt_s(seconds: float) -> str:
    """``7 s``, ``45 s``, ``2 min 05 s`` -- a clock a person reads."""
    total = max(0, int(seconds))
    if total < 60:
        return f"{total} s"
    return f"{total // 60} min {total % 60:02d} s"


@dataclass(frozen=True)
class State:
    phase: str = IDLE
    model: str = ""
    tool: str = ""
    arg: str = ""
    person: str = ""            # APPROVAL | ANSWER while phase == PERSON
    phase_s: float = 0.0        # seconds in this phase
    turn_s: float = 0.0         # seconds since the turn began (or its length)
    silent_s: float = 0.0       # seconds since the model last produced
    first: bool = True          # no token yet in the current request
    outcome: str = ""           # done | stopped | failed, when phase == DONE
    seq: int = 0                # changes whenever the phase changes

    @property
    def active(self) -> bool:
        return self.phase not in (IDLE, DONE)


class Activity:
    """The engine's live phase. Thread-safe; written by the turn thread
    and the asking thread, read by any surface's timer."""

    def __init__(self, now: Callable[[], float] = time.monotonic) -> None:
        self._now = now
        self._lock = threading.RLock()
        self._phase = IDLE
        self._model = ""
        self._tool = ""
        self._arg = ""
        self._person = ""
        self._phase_t0 = 0.0
        self._turn_t0 = 0.0
        self._last_output = 0.0
        self._first = True
        self._outcome = ""
        self._turn_len = 0.0
        self._depth = 0
        self._tools_open = 0
        self._held: list[tuple] = []
        self._seq = 0

    # -- writers ---------------------------------------------------------
    def _enter(self, phase: str, now: float) -> None:
        self._phase = phase
        self._phase_t0 = now
        self._seq += 1

    def begin_turn(self, model: str = "") -> None:
        """A turn starts with a request to the model. Nested calls (a
        continuation run inside the same turn) count as a new request,
        not a new turn: the turn clock keeps running."""
        with self._lock:
            now = self._now()
            self._depth += 1
            if model:
                self._model = str(model)
            if self._depth == 1:
                self._turn_t0 = now
                self._tools_open = 0
                self._held = []
                self._outcome = ""
            self._first = True
            self._last_output = now
            self._enter(MODEL_WAIT, now)

    def request(self) -> None:
        """A new request is out (the round after a tool result)."""
        with self._lock:
            now = self._now()
            self._first = True
            self._last_output = now
            self._enter(MODEL_WAIT, now)

    def output(self) -> None:
        """The model produced something: a token, reasoning, a call."""
        with self._lock:
            now = self._now()
            self._last_output = now
            self._first = False
            if self._phase == MODEL_WAIT:
                self._enter(STREAMING, now)

    def tool_start(self, name: str, tool_input: Any = None) -> None:
        with self._lock:
            now = self._now()
            self._last_output = now
            self._first = False
            self._tools_open += 1
            self._tool = bare_tool_name(name)
            self._arg = short_arg(tool_input)
            if self._phase == PERSON:
                return          # the answer puts the tool back
            self._enter(TOOL, now)

    def tool_end(self, name: str = "") -> None:
        """A call returned; the next thing is the request carrying its
        result, so the indicator goes back to waiting for the model."""
        with self._lock:
            self._tools_open = max(0, self._tools_open - 1)
            if self._phase == PERSON or self._tools_open:
                return
            if self._phase in (TOOL, STREAMING, MODEL_WAIT):
                self.request()

    @contextmanager
    def awaiting(self, kind: str = APPROVAL) -> Iterator[None]:
        """The turn is blocked on the person for the duration."""
        with self._lock:
            now = self._now()
            self._held.append((self._phase, self._phase_t0, self._person))
            self._person = ANSWER if kind == ANSWER else APPROVAL
            self._enter(PERSON, now)
        try:
            yield
        finally:
            with self._lock:
                phase, t0, person = (self._held.pop() if self._held
                                     else (IDLE, self._now(), ""))
                if self._phase == PERSON:
                    # The interrupted phase with its own start time: a
                    # tool that asked keeps counting from its dispatch.
                    self._phase, self._phase_t0 = phase, t0
                    self._person = person
                    self._seq += 1

    def end_turn(self, *, error: str = "", stopped: bool = False) -> None:
        with self._lock:
            self._depth = max(0, self._depth - 1)
            if self._depth:
                return
            now = self._now()
            self._turn_len = max(0.0, now - self._turn_t0)
            self._outcome = ("failed" if error else
                             "stopped" if stopped else "done")
            self._tools_open = 0
            self._held = []
            self._enter(DONE, now)

    # -- reader ----------------------------------------------------------
    def snapshot(self) -> State:
        with self._lock:
            now = self._now()
            turn_s = (self._turn_len if self._phase == DONE
                      else max(0.0, now - self._turn_t0)
                      if self._phase != IDLE else 0.0)
            return State(
                phase=self._phase, model=self._model, tool=self._tool,
                arg=self._arg, person=self._person,
                phase_s=max(0.0, now - self._phase_t0) if self._phase_t0 else 0.0,
                turn_s=turn_s,
                silent_s=max(0.0, now - self._last_output)
                if self._last_output else 0.0,
                first=self._first, outcome=self._outcome, seq=self._seq)


def resolve(source: Any) -> "Activity | None":
    """An ``Activity`` from what a surface holds: the object itself, a
    zero-argument callable returning one (the dashboard rebuilds its
    engine), or anything carrying ``.activity``. None when there is none."""
    try:
        if isinstance(source, Activity):
            return source
        if callable(source):
            got = source()
            return got if isinstance(got, Activity) else resolve_attr(got)
        return resolve_attr(source)
    except Exception:
        return None


def resolve_attr(holder: Any) -> "Activity | None":
    got = getattr(holder, "activity", None) if holder is not None else None
    return got if isinstance(got, Activity) else None


@contextmanager
def awaiting_person(source: Any, kind: str = APPROVAL) -> Iterator[None]:
    """``Activity.awaiting`` for a surface that may hold no activity."""
    target = resolve(source)
    if target is None:
        yield
        return
    with target.awaiting(kind):
        yield


def line(state: State) -> str:
    """The indicator text for one snapshot; "" when idle.

    ``waiting for glm-4.6 … 45 s`` (no token of this request yet),
    ``glm-4.6 is answering … 12 s``, ``running bash: pytest -q … 12 s``,
    ``waiting for your approval … 8 s``, ``done in 31 s``.
    """
    model = state.model or "the model"
    clock = fmt_s(state.phase_s)
    if state.phase == MODEL_WAIT:
        return f"waiting for {model} … {clock}"
    if state.phase == STREAMING:
        return f"{model} is answering … {clock}"
    if state.phase == TOOL:
        what = f"{state.tool}: {state.arg}" if state.arg else state.tool
        return f"running {what or 'a tool'} … {clock}"
    if state.phase == PERSON:
        whom = "answer" if state.person == ANSWER else "approval"
        return f"waiting for your {whom} … {clock}"
    if state.phase == DONE:
        verb = {"stopped": "stopped after", "failed": "failed after"}.get(
            state.outcome, "done in")
        return f"{verb} {fmt_s(state.turn_s)}"
    return ""


def tool_expected_s(tool: str) -> float:
    """What a call of this tool usually takes before it counts as slow."""
    return _env_s("DELFIN_SLOW_TOOL_S",
                  _LONG_TOOLS.get(bare_tool_name(tool), TOOL_SLOW_S))


def slow_hint(state: State, surface: str = "terminal") -> str:
    """One sentence once a phase is past its threshold: what is slow and
    how to stop it on THIS surface. "" otherwise.

    ``surface`` is ``terminal`` (Esc / Ctrl-C) or ``dashboard`` (Stop).
    A tool call cannot be cut by either: the sentence says the stop takes
    effect when the call returns, because that is what it does.
    """
    stop = ("press Esc to end the turn" if surface == "terminal"
            else "press Stop to end the turn")
    model = state.model or "the model"
    if state.phase == MODEL_WAIT:
        limit = _env_s("DELFIN_SLOW_FIRST_TOKEN_S", FIRST_TOKEN_SLOW_S)
        if state.phase_s >= limit:
            return (f"slow: {model} has sent no token for "
                    f"{fmt_s(state.phase_s)} (the endpoint is queueing or "
                    f"prefilling) — wait, or {stop}.")
    elif state.phase == STREAMING:
        limit = _env_s("DELFIN_SLOW_FIRST_TOKEN_S", SILENT_STREAM_SLOW_S)
        if state.silent_s >= limit:
            return (f"slow: {model} went silent {fmt_s(state.silent_s)} ago "
                    f"mid-answer — wait, or {stop}.")
    elif state.phase == TOOL:
        limit = tool_expected_s(state.tool)
        if state.phase_s >= limit:
            after = ("; Ctrl-C twice more leaves without waiting for it"
                     if surface == "terminal" else "")
            return (f"slow: {state.tool or 'the tool'} has run for "
                    f"{fmt_s(state.phase_s)}, longer than the usual "
                    f"{fmt_s(limit)} — {stop} once it returns{after}.")
    return ""
