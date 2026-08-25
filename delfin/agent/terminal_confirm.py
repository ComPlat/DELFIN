"""Asking the person at the terminal, from the thread doing the work.

Without a confirm callback the gate does NOT simply refuse. Writes inside
the workspace are allowed silently, and `exit_plan_mode` approves itself —
so a terminal agent with no broker is not "no prompts", it is unattended
editing. This module is what makes the terminal attended.

**The callback must stay a bound method.** ``_confirm_timed_out`` reads
``perms.confirm_callback.__self__.last_timed_out`` to tell an expired
dialog from a refusal. Wrapping it in a lambda destroys that distinction
silently, and the gate then records every timeout as a real denial in
``denied_actions`` — permanently closing paths the user never saw.

**Threading.** Requests arrive on whatever thread is executing tools: the
turn worker, a subagent, a background job. They queue here and block on
their own Event. Exactly one consumer — the main thread, which owns the
terminal — renders and answers. The lock is never held across a read or
across ``event.set()``, and an answer arriving after its request expired
is discarded rather than applied to whatever the agent has moved on to.

**No timeout by default.** ``KitConfirmBroker`` waits 300 s because a
browser tab can be abandoned; a blinking terminal prompt is attended by
definition. Worse, an expiry sets ``last_timed_out``, which makes the gate
NOT record the refusal — so the model retries and blocks again.
"""

from __future__ import annotations

import itertools
import threading
from dataclasses import dataclass, field
from typing import Any, Callable

from . import repl_render as rr

__all__ = [
    "ConfirmRequest", "TerminalConfirmBroker", "Option",
    "options_for", "render_request", "CONFIRM", "ASK", "PLAN",
]

CONFIRM = "confirm"
ASK = "ask"
PLAN = "plan"

_WRITE_TOOLS = frozenset({
    "write_file", "edit_file", "multi_edit", "apply_patch", "notebook_edit",
    "edit_sheet", "create_docx", "create_pdf", "fill_pdf_form",
    "fill_docx_template", "fill_series", "merge_pdfs", "split_pdf",
})
_SHELL_TOOLS = frozenset({"bash", "bash_background"})

# Markers the gate puts at the head of its own preview. They decide which
# options may be offered, so they are read rather than re-derived.
_GUARD_MARKER = "[SELF-MODIFICATION GUARD]"
_CALC_MARKER = "[CALC EDIT]"
_OUTSIDE_MARKER = "[OUTSIDE-WORKSPACE READ]"


@dataclass(frozen=True)
class Option:
    key: str
    label: str
    detail: str = ""


@dataclass
class ConfirmRequest:
    kind: str
    tool: str = ""
    args: dict = field(default_factory=dict)
    preview: str = ""
    payload: dict = field(default_factory=dict)
    seq: int = 0
    event: threading.Event = field(default_factory=threading.Event)
    decision: Any = None
    resolved: bool = False
    expired: bool = False

    @property
    def command(self) -> str:
        return str(self.args.get("command", "") or "")

    @property
    def is_shell(self) -> bool:
        return self.tool in _SHELL_TOOLS

    @property
    def is_write(self) -> bool:
        return self.tool in _WRITE_TOOLS

    @property
    def is_protected(self) -> bool:
        return (_GUARD_MARKER in self.preview) or (_CALC_MARKER in self.preview)

    @property
    def is_outside_read(self) -> bool:
        return _OUTSIDE_MARKER in self.preview


def options_for(req: ConfirmRequest, *, suggestion: str = "") -> list[Option]:
    """What may be offered for this request, and nothing more.

    Two absences are the point rather than an oversight.

    An outside-workspace READ gets no "always". The only persistable form
    of that grant is an extra workspace directory, and that is WRITABLE —
    in this session and every future one. One keystroke on a read prompt
    must not be able to hand over write access.

    An unknown MCP tool gets no "always" either: there is no persist kind
    for one. ``allow_patterns`` are shell regexes, so an "always" here
    would have to be a session allowlist this module invented — new state,
    new bypass, and nothing asking for it.
    """
    if req.kind == PLAN:
        return [
            Option("d", "approve → default",
                   "ask before each write and each command"),
            Option("e", "approve → acceptEdits",
                   "writes go through; shell still asked"),
            Option("n", "reject", "stay in plan and say why"),
        ]

    opts = [Option("y", "yes", "this one action"),
            Option("n", "no", "refused, and the model is told not to retry")]

    if req.is_shell and req.command:
        opts.append(Option("A", "always this exact command",
                           "persisted, globally, for every future session"))
        if suggestion and suggestion != _exact_pattern(req.command):
            opts.append(Option("k", "always commands like this",
                               f"persists {suggestion}"))
    if req.is_write and not req.is_protected:
        opts.append(Option("e", "stop asking about writes this session",
                           "acceptEdits until this session ends"))

    opts.append(Option("a", "abort", "deny this and end the turn"))
    opts.append(Option("?", "explain", "what each option does here"))
    return opts


def _exact_pattern(command: str) -> str:
    import re
    return r"^\s*" + re.escape(command.strip()) + r"\s*$"


def render_request(req: ConfirmRequest, *, theme: rr.Theme | None = None,
                   width: int = 100, options: list[Option] | None = None,
                   body_lines: int = 24) -> str:
    """The frame the user reads before deciding.

    Everything from the tool is stripped of control characters first. A
    file whose contents draw a convincing approval prompt, or clear the
    screen and print "Approved.", must be text and nothing else.
    """
    theme = theme or rr.Theme()
    options = options if options is not None else options_for(req)
    name = rr.short_tool_name(req.tool) or req.tool or "?"

    head = f"┌─ {theme.bold(name)}"
    if req.is_protected:
        head = theme.red(f"┌─ {name}  PROTECTED PATH")
    elif req.is_outside_read:
        head = theme.yellow(f"┌─ {name}  OUTSIDE the workspace")

    lines = [head]
    body = rr.strip_control(req.preview or "").splitlines()
    shown, hidden = body[:body_lines], max(0, len(body) - body_lines)
    for line in shown:
        lines.append("│  " + rr.truncate_middle(line, max(20, width - 4)))
    if hidden:
        lines.append(theme.dim(f"│  … {hidden} more lines"))

    if req.is_outside_read:
        lines.append("│")
        lines.append(theme.dim(
            "│  Approving reads this file and opens its directory for READS "
            "for the rest of this session. Nothing is saved."))

    keys = "  ".join(f"[{o.key}] {o.label}" for o in options)
    lines.append("└─ " + rr.truncate_middle(keys, max(20, width - 3)))
    return "\n".join(lines)


def render_help(req: ConfirmRequest, options: list[Option]) -> str:
    out = [f"  [{o.key}] {o.label}" + (f" — {o.detail}" if o.detail else "")
           for o in options]
    if req.is_outside_read:
        out.append("  no 'always' here: the only persistable form of this "
                   "grant is a writable directory, now and in every future "
                   "session.")
    return "\n".join(out)


class TerminalConfirmBroker:
    """The three UI callbacks, answered by whoever is at the terminal."""

    def __init__(self, *, timeout_s: float = 0.0,
                 persist: Callable[[str], tuple[bool, str]] | None = None,
                 set_mode: Callable[[str], Any] | None = None,
                 on_abort: Callable[[], None] | None = None) -> None:
        self.timeout_s = float(timeout_s or 0.0)
        self._persist = persist
        self._set_mode = set_mode
        self._on_abort = on_abort
        self._lock = threading.Lock()
        self._queue: list[ConfirmRequest] = []
        self._seq = itertools.count(1)
        self.last_timed_out = False
        self.aborted = False

    # -- the three bindings ----------------------------------------------
    def callback(self, tool_name: str, args: dict, preview: str) -> bool:
        """Bound on purpose: the gate reads last_timed_out off __self__."""
        req = self._enqueue(ConfirmRequest(
            kind=CONFIRM, tool=str(tool_name or ""),
            args=dict(args or {}), preview=str(preview or "")))
        return bool(self._wait(req))

    def ask_user(self, payload: dict) -> dict:
        req = self._enqueue(ConfirmRequest(kind=ASK, payload=dict(payload or {})))
        answer = self._wait(req)
        return answer if isinstance(answer, dict) else {"answers": []}

    def approve_plan(self, plan: str) -> dict:
        req = self._enqueue(ConfirmRequest(kind=PLAN, preview=str(plan or "")))
        answer = self._wait(req)
        if isinstance(answer, dict):
            return answer
        return {"approved": False, "new_mode": "plan"}

    # -- producer side ----------------------------------------------------
    def _enqueue(self, req: ConfirmRequest) -> ConfirmRequest:
        req.seq = next(self._seq)
        with self._lock:
            if self.aborted:
                req.resolved = True
                req.decision = self._refusal_for(req)
                return req
            self._queue.append(req)
        return req

    def _wait(self, req: ConfirmRequest) -> Any:
        if req.resolved:
            return req.decision
        got = req.event.wait(self.timeout_s or None)
        with self._lock:
            if not got and not req.resolved:
                req.expired = True
                req.resolved = True
                req.decision = self._refusal_for(req)
                self.last_timed_out = True
                if req in self._queue:
                    self._queue.remove(req)
            else:
                self.last_timed_out = False
        return req.decision

    @staticmethod
    def _refusal_for(req: ConfirmRequest) -> Any:
        if req.kind == PLAN:
            return {"approved": False, "new_mode": "plan"}
        if req.kind == ASK:
            return {"answers": []}
        return False

    # -- consumer side (the main thread, and only it) ---------------------
    def take(self) -> ConfirmRequest | None:
        with self._lock:
            for req in self._queue:
                if not req.resolved:
                    return req
        return None

    def resolve(self, req: ConfirmRequest, decision: Any) -> bool:
        """Apply an answer. False when the request had already expired."""
        with self._lock:
            if req.resolved:
                return False
            req.resolved = True
            req.decision = decision
            if req in self._queue:
                self._queue.remove(req)
        req.event.set()
        return True

    def abort_all(self) -> list[ConfirmRequest]:
        """Deny everything in flight, and everything that arrives after."""
        with self._lock:
            self.aborted = True
            pending = [r for r in self._queue if not r.resolved]
            for req in pending:
                req.resolved = True
                req.decision = self._refusal_for(req)
            self._queue.clear()
        for req in pending:
            req.event.set()
        if self._on_abort:
            try:
                self._on_abort()
            except Exception:
                pass
        return pending

    def reset_abort(self) -> None:
        with self._lock:
            self.aborted = False

    # -- the two side effects an answer can have --------------------------
    def suggest_pattern(self, command: str) -> str:
        """Never invents one: kit_settings already encodes the policy.

        Its rules are what stop `git push -u origin feature` persisting as
        `^\\s*git\\s+push\\b`, which would also cover
        `git push --force-with-lease origin main`.
        """
        try:
            from . import kit_settings
            return str(kit_settings.suggest_pattern_for_command(command) or "")
        except Exception:
            return ""

    def exact_pattern(self, command: str) -> str:
        return _exact_pattern(command)

    def persist(self, pattern: str) -> tuple[bool, str]:
        if not self._persist:
            return False, "nothing to persist through"
        try:
            return self._persist(pattern)
        except Exception as exc:
            return False, f"persist failed: {exc}"

    def accept_edits(self) -> tuple[bool, str]:
        if not self._set_mode:
            return False, "this backend has no permission gate"
        try:
            self._set_mode("acceptEdits")
        except Exception as exc:
            return False, f"could not switch mode: {exc}"
        return True, "writes go through for the rest of this session"
