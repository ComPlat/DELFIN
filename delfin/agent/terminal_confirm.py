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
import json
import os
import socket
import threading
import time
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

from . import repl_render as rr

__all__ = [
    "ConfirmRequest", "TerminalConfirmBroker", "Option",
    "options_for", "render_request", "CONFIRM", "ASK", "PLAN",
    "pending_at_terminals",
]

#: Where a question waiting at a terminal is published, whole, for
#: whoever is supervising -- and where an answer for it may be left. The
#: room is file_confirm's, so its checks apply unchanged and the terminal
#: still wins a race: resolve() is the single point that decides.
_PENDING_DIR = Path.home() / ".delfin" / "terminal_confirmations"

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
    #: Where this question was published for a supervisor to read, or
    #: None. Carried on the request so the answer can take it away again.
    published: Any = None

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


def _publish_pending(req: ConfirmRequest, session_id: str,
                     session_key: str) -> "Path | None":
    """Write the whole question where a supervisor can read it.

    The pane renders it capped at 24 lines and cut to the pane's width,
    which is how an operator came to face a compound command with its
    middle missing and a guard diff saying "… 30 more lines". The comment
    above ``render_request`` gives the rule this restores: an approval you
    cannot read is one you cannot give. Here the preview goes out whole.

    Written under a temporary name and moved into place, so a reader
    never sees half a question.
    """
    room = _PENDING_DIR
    room.mkdir(parents=True, exist_ok=True)
    try:
        os.chmod(room, 0o700)
    except OSError:
        pass
    record = {
        "id": f"{int(time.time())}-{uuid.uuid4().hex[:8]}",
        "session_id": str(session_id or ""),
        "session_key": str(session_key or ""),
        "kind": str(req.kind or ""),
        "tool": str(req.tool or ""),
        "command": req.command,
        "preview": str(req.preview or ""),
        "protected": bool(req.is_protected),
        "outside_read": bool(req.is_outside_read),
        "asked_at": time.time(),
        "pid": os.getpid(),
        # A pid is a number the system hands out again. The start settles
        # which process it was -- one reader for that fact, the one
        # extracted for exactly this class of mistake.
        "proc_start": _proc_start_of_this_process(),
        "host": socket.gethostname(),
        # Answerable from both ends now. The terminal still wins a race:
        # resolve() is the single point that decides, and the first
        # answer takes it.
        "answer_at": "terminal or supervisor",
    }
    # The same filename file_confirm uses, so its reader, its lister and
    # its answer writer apply here unchanged -- with every careful part
    # they carry: not a symlink, this user's, not group-writable, naming
    # this very question, not older than it.
    path = room / f"{record['id']}.request.json"
    tmp = room / f".{record['id']}.partial"
    tmp.write_text(json.dumps(record, indent=2), encoding="utf-8")
    try:
        os.chmod(tmp, 0o600)
    except OSError:
        pass
    tmp.replace(path)
    return path


def _proc_start_of_this_process() -> str:
    """When this process began, or "". Never raises."""
    try:
        from . import proc_identity
        return proc_identity.process_start(os.getpid())
    except Exception:
        return ""


def _asker_is_gone(record: dict) -> bool:
    """True only when the asking process is provably no longer there.

    Where it cannot be asked -- another machine, a kernel that will not
    say, a record with no pid -- the question STAYS. Guessing a session
    dead costs somebody their running work, and that is the expensive
    direction to be wrong in.
    """
    try:
        pid = int(record.get("pid") or 0)
    except (TypeError, ValueError):
        return False
    if pid <= 0:
        # No pid is not a dead pid. ``alive`` answers False for a number
        # it cannot use, which is the right answer to ITS question and
        # the wrong one to this one.
        return False
    try:
        from . import proc_identity
        return proc_identity.alive(pid,
                                   str(record.get("proc_start") or ""),
                                   str(record.get("host") or "")) is False
    except Exception:
        return False


def _withdraw_pending(path) -> None:
    """Take a published question away once it has been answered."""
    if not path:
        return
    try:
        Path(path).unlink()
    except OSError:
        pass


def pending_at_terminals() -> list[dict]:
    """Every question a terminal session is waiting on, read whole.

    Never raises: a supervisor's view of the world must not be the thing
    that ends the supervisor.
    """
    try:
        from . import file_confirm as _fc
        rows = list(_fc.pending(_PENDING_DIR))
    except Exception:
        return []
    out: list[dict] = []
    for record in rows:
        if _asker_is_gone(record):
            # Nobody will ever answer it, and read as open it is worse
            # than not published at all -- eight of these stood in the
            # list within minutes of the watch going live, the oldest
            # for 3190 seconds.
            _withdraw_pending(record.get("_file"))
            continue
        out.append(record)
    return out


def _note_outside_answer(req: ConfirmRequest, session_key: str,
                         decision) -> None:
    """Record that this one was answered from outside the pane.

    An approval given where nobody at the terminal saw it is exactly the
    kind of thing the containment panel exists to show.
    """
    try:
        from .api_client import _record_security_event
        _record_security_event(
            "approval_from_outside", req.tool or "?",
            f"{session_key or '?'}: {'approved' if decision else 'refused'} "
            f"{'(PROTECTED) ' if req.is_protected else ''}from outside the "
            f"terminal", blocked=False)
    except Exception:
        pass


def answer_waiting(request_id: str, decision: bool, *, by: str = "",
                   reason: str = "") -> bool:
    """Answer a question a terminal session is stuck on. Never raises.

    The decision this reverses is the one the first published record
    stated in itself -- "answer_at: terminal" -- which was right while
    there was no checked way in. There is one now, and it is not new:
    file_confirm answers headless sessions with exactly these checks.

    Reading beats guessing. The supervisor decides on the WHOLE preview,
    which is more than the pane was ever able to show. A refusal may say
    what to do instead; the reason reaches the model in the same turn as
    the refusal, so it does not have to guess the same file differently.
    """
    try:
        from . import file_confirm as _fc
        return bool(_fc.answer(str(request_id), bool(decision),
                               room=_PENDING_DIR, by=str(by or ""),
                               reason=str(reason or "")))
    except Exception:
        return False


def _answer_from_outside(req: ConfirmRequest) -> "tuple[bool, str] | None":
    """The answer somebody left for *req*, with a refusal reason, or None.

    Never raises. The reason comes from file_confirm's single checked
    reader -- there is deliberately no laxer second path for text.
    """
    published = getattr(req, "published", None)
    if not published:
        return None
    try:
        from . import file_confirm as _fc
        path = Path(published)
        request_id = path.name[:-len(".request.json")]
        answer_path = path.with_name(f"{request_id}.answer.json")
        return _fc._read_answer(answer_path, request_id,
                                path.stat().st_mtime)
    except Exception:
        return None


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
    # Masked before anything reaches the screen: this broker printed the
    # command whole, so a line carrying a token showed it to the terminal
    # and to the scrollback behind it. The value goes, the command stays
    # readable -- an approval you cannot read is one you cannot give.
    try:
        from .cli_approve import redact_preview as _redact_preview
        _preview = _redact_preview(req.preview or "")
    except Exception:
        _preview = rr.strip_control(req.preview or "")
    body = _preview.splitlines()
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
                 on_abort: Callable[[], None] | None = None,
                 session_id: str = "", session_key: str = "",
                 poll_s: float = 0.25,
                 on_activity: "Callable[[], Any] | None" = None) -> None:
        self.timeout_s = float(timeout_s or 0.0)
        # How often the waiting thread looks for an answer left outside
        # the pane. It is the tool thread that waits, so the look costs
        # nothing the terminal would otherwise be doing.
        self.poll_s = max(0.01, float(poll_s or 0.25))
        # Asking is proof of life. Presence was refreshed only from the
        # idle poll, and a session going from one approval into the next
        # never reaches it -- so it aged out of the overview exactly
        # while it was waiting for the supervisor. The broker does not
        # learn about presence for this: it calls back.
        self.on_activity = on_activity
        # Only so a published question can name the session it belongs
        # to. Empty is fine: the record is then anonymous, not absent.
        self.session_id = str(session_id or "")
        self.session_key = str(session_key or "")
        self._persist = persist
        self._set_mode = set_mode
        self._on_abort = on_abort
        self._lock = threading.Lock()
        self._queue: list[ConfirmRequest] = []
        self._seq = itertools.count(1)
        # Per THREAD, not per broker. Requests arrive on whichever
        # thread is running tools -- the turn worker, a subagent, a
        # background job -- and they overlap. With one flag for all of
        # them, an answer to one request cleared the expiry of another,
        # and the expired one was then recorded as a REFUSAL: a path
        # permanently closed that the user never saw. The flag is read
        # off ``__self__`` by the gate, so it stays an attribute and
        # keeps its name; only its storage moved.
        self._timed_out = threading.local()
        # Per thread for the same reason as _timed_out: requests overlap
        # across the turn worker, subagents and background jobs, and a
        # reason from one thread's refusal must not land on another's.
        self._refusal_reason = threading.local()
        self.aborted = False

    # -- the three bindings ----------------------------------------------
    @property
    def last_timed_out(self) -> bool:
        """Whether THIS thread's last dialog expired rather than being
        refused. See the note in __init__ for why it is per thread."""
        return bool(getattr(self._timed_out, "value", False))

    @last_timed_out.setter
    def last_timed_out(self, value: bool) -> None:
        self._timed_out.value = bool(value)

    @property
    def last_refusal_reason(self) -> str:
        """Why THIS thread's last dialog was refused, if it said why.
        Set by the thread that asked, read by the gate in the same turn,
        and cleared by every new dialog -- a reason never outlives its
        refusal. Per thread like last_timed_out, and for the same reason."""
        return str(getattr(self._refusal_reason, "value", "") or "")

    @last_refusal_reason.setter
    def last_refusal_reason(self, value: str) -> None:
        self._refusal_reason.value = str(value or "")

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
        try:
            req.published = _publish_pending(req, self.session_id,
                                             self.session_key)
        except Exception:
            # A supervisor's convenience does not get to fail an
            # approval: the prompt is asked either way.
            req.published = None
        if self.on_activity is not None:
            try:
                self.on_activity()
            except Exception:
                pass            # a heartbeat may never cost a question
        return req

    def _wait(self, req: ConfirmRequest) -> Any:
        # A reason never outlives its refusal: every new dialog starts
        # from none, whatever the last one carried -- also when it was
        # resolved before waiting, by an abort.
        self.last_refusal_reason = ""
        if req.resolved:
            return req.decision
        import time as _time
        deadline = (_time.monotonic() + self.timeout_s) if self.timeout_s else None
        got = False
        while True:
            slice_s = self.poll_s
            if deadline is not None:
                slice_s = min(slice_s, max(0.0, deadline - _time.monotonic()))
            got = req.event.wait(slice_s if slice_s > 0 else 0.001)
            if got or req.resolved:
                got = True
                break
            outside = _answer_from_outside(req)
            if outside is not None and self.resolve(req, outside[0]):
                if outside[0] is False and outside[1]:
                    self.last_refusal_reason = outside[1]
                _note_outside_answer(req, self.session_key, outside[0])
                self._audit_dialog(req, by="outside")
                return req.decision
            if deadline is not None and _time.monotonic() >= deadline:
                break
        with self._lock:
            if not got and not req.resolved:
                req.expired = True
                req.resolved = True
                req.decision = self._refusal_for(req)
                self.last_timed_out = True
                if req in self._queue:
                    self._queue.remove(req)
                try:
                    _withdraw_pending(req.published)
                except Exception:
                    pass
                req.published = None
            else:
                self.last_timed_out = False
        self._audit_dialog(req, by="expired" if req.expired else "terminal")
        return req.decision

    def _audit_dialog(self, req: ConfirmRequest, *, by: str) -> None:
        """One audit record per dialog a person had to answer.

        The audit log recorded what ran and what was refused, never that
        someone was ASKED: a round's operator load could only be counted
        by hand (night run 2026-09-25). Carries the session name (-n),
        so a report can filter by it. Never raises.
        """
        try:
            from . import audit_log as _audit
            decision = req.decision
            _audit.append({
                "event": "dialog",
                "session_id": self.session_id,
                "session_key": self.session_key,
                "tool": str(getattr(req, "tool", "") or ""),
                "answer": ("expired" if req.expired else
                           "approved" if decision is True else
                           "denied" if decision is False else "answered"),
                "by": by,
            })
        except Exception:
            pass

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
        try:
            _withdraw_pending(req.published)
        except Exception:
            pass        # a stale record is a nuisance; a raise here is a bug
        req.published = None
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
            try:
                _withdraw_pending(req.published)
            except Exception:
                pass
            req.published = None
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
