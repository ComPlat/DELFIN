"""The session's approval questions, asked and answered over the protocol.

The terminal REPL asks through ``TerminalConfirmBroker``: the gate's three
callbacks queue a request and block the asking thread until the main thread
answers it from the keyboard. Everything that matters about asking lives in
that broker -- per-thread expiry and refusal reasons, the audit record of
every dialog, questions published for ``delfin-agent approvals``, the
persisted "always" patterns -- so the worker uses the same broker and only
changes who answers:

* every request is announced as a ``permission.asked`` event, with exactly
  the options ``terminal_confirm.options_for`` gives the REPL;
* a ``permission.reply`` from a window resolves it, with the same meaning
  each key has in ``TerminalAgent._answer_request``;
* every way a request can end -- answered here, answered from outside,
  aborted -- is announced as ``permission.resolved``.

Nothing is ever answered on the user's behalf. A request with no window
attached waits; ``timeout_s`` stays 0, as it is for a terminal.

The REPL answers on its main thread, which takes each request off the
queue (``take``) and holds the turn's output until it is answered. Here a
window answers from the protocol's reader thread instead, so ``take``
never hands anything to the main thread: it keeps streaming the turn --
and any other thread's output -- while a question is open. What the
broker announces goes through ``emit`` in stream order (the worker queues
it behind the turn's output), so a question never arrives ahead of the
tool call it is about.
"""

from __future__ import annotations

import itertools
import secrets
import threading
from typing import Any, Callable

from .. import terminal_confirm as tc

Emit = Callable[..., None]


class UnknownRequest(LookupError):
    """No open request has that id (never asked, or already resolved)."""


class InvalidAnswer(ValueError):
    """The answer is not one the request offered; nothing was applied."""


def _preview(req: tc.ConfirmRequest) -> str:
    """The request's preview with secrets masked, as the REPL shows it."""
    try:
        from ..cli_approve import redact_preview
        return redact_preview(req.preview or "")
    except Exception:
        from .. import repl_render as rr
        return rr.strip_control(req.preview or "")


def _purpose(req: tc.ConfirmRequest) -> str:
    """The call's own one-line description, masked twice: the preview's
    line masking and the service redactor, which also finds a bare key."""
    text = str((req.args or {}).get("description") or "")
    if not text.strip():
        return ""
    from .redact import Redactor
    masked = Redactor().text(_preview(tc.ConfirmRequest(kind=req.kind, preview=text)))
    return " ".join(masked.split())[:200]


def _choices(payload: dict) -> list[tuple[str, str]]:
    """``(label, description)`` per option of an ask_user_question.

    The engine sends ``{"label", "description"}`` objects (plain strings are
    accepted too); a window shows the label and answers with it, which is
    what ``delfin-agent approvals`` and the model expect.
    """
    out: list[tuple[str, str]] = []
    for option in (payload.get("options") or [])[:9]:
        if isinstance(option, dict):
            label = str(option.get("label", "") or "")
            detail = str(option.get("description", "") or "")
        else:
            label, detail = str(option), ""
        if label:
            out.append((label, detail))
    return out


class ProtocolBroker(tc.TerminalConfirmBroker):
    """``TerminalConfirmBroker`` whose questions are answered by a window."""

    def __init__(self, emit: Emit, *, stop_turn: Callable[[], None] | None = None,
                 on_waiting: Callable[[bool], None] | None = None,
                 on_feedback: Callable[[str], None] | None = None, **kwargs):
        super().__init__(timeout_s=0.0, **kwargs)
        self._emit = emit
        self._stop_turn = stop_turn
        self._on_waiting = on_waiting
        self._on_feedback = on_feedback
        self._ids = itertools.count(1)
        # Ids are unique per session, across restarts of its worker: a window
        # holding an old id must never answer a newer question with it.
        self._id_prefix = secrets.token_hex(3)
        self._by_id: dict[str, tc.ConfirmRequest] = {}
        # Held while a request is registered AND announced, and while one is
        # withdrawn AND announced: "asked" and "resolved" for one request
        # can then only be emitted in that order.
        self._ids_lock = threading.RLock()
        self._answered_by = threading.local()

    # -- announcing -----------------------------------------------------------
    def request_id(self, req: tc.ConfirmRequest) -> str:
        return str(getattr(req, "protocol_id", "") or "")

    def pending(self) -> list[tc.ConfirmRequest]:
        with self._ids_lock:
            return [r for r in self._by_id.values() if not r.resolved]

    def take(self) -> tc.ConfirmRequest | None:
        """Nothing for the main thread: a window answers (module doc)."""
        return None

    # -- closing: refused for good -----------------------------------------------
    closed = False

    def close(self) -> list[tc.ConfirmRequest]:
        """Refuse everything open and everything asked from now on, for the
        rest of the process: unlike an abort, no new turn lifts it."""
        self.closed = True
        return self.abort_all()

    def reset_abort(self) -> None:
        # A turn that starts while the session closes must not reopen the
        # gate the shutdown just closed (its question would wait forever).
        if not self.closed:
            super().reset_abort()

    def _enqueue(self, req: tc.ConfirmRequest) -> tc.ConfirmRequest:
        req = super()._enqueue(req)
        if req.resolved:              # refused on arrival by an abort
            return req
        rid = f"p{self._id_prefix}-{next(self._ids)}"
        with self._ids_lock:
            req.protocol_id = rid
            self._by_id[rid] = req
            self._emit("permission.asked", **self.describe(req))
        if req.resolved:
            # An abort landed between queueing and announcing: the window
            # must hear that this question is already over.
            self._announce_resolved(req, "abort")
        self._waiting_changed()
        return req

    def describe(self, req: tc.ConfirmRequest) -> dict[str, Any]:
        """Everything a window needs to draw and answer ``req``."""
        suggestion = ""
        if req.kind == tc.CONFIRM and req.is_shell and req.command:
            try:
                suggestion = self.suggest_pattern(req.command)
            except Exception:
                suggestion = ""
        # Frozen at ask time: "k" persists exactly what the window showed.
        req.protocol_suggestion = suggestion
        options = tc.options_for(req, suggestion=suggestion)
        fields: dict[str, Any] = {
            "requestId": self.request_id(req),
            "kind": req.kind,
            "tool": req.tool,
            "preview": _preview(req),
            "options": [{"key": o.key, "label": o.label, "detail": o.detail}
                        for o in options],
            "flags": {"shell": req.is_shell, "write": req.is_write,
                      "protected": req.is_protected,
                      "outsideRead": req.is_outside_read},
        }
        if suggestion:
            fields["suggestion"] = suggestion
        editable = tc.editable_command(req)
        if editable:
            # What key c puts in the prompt to change; absent when the
            # command holds a secret (tc.editable_command).
            fields["command"] = editable
        purpose = _purpose(req)
        if purpose:
            fields["purpose"] = purpose
        read_dir = tc.read_dir_for(req) if req.is_outside_read else ""
        if read_dir:
            fields["readDir"] = read_dir
        if req.kind == tc.ASK:
            payload = req.payload or {}
            fields["question"] = str(payload.get("question", "") or "")
            header = str(payload.get("header", "") or "")
            if header:
                fields["header"] = header
            choices = _choices(payload)
            fields["choices"] = [label for label, _ in choices]
            fields["choiceDetails"] = [detail for _, detail in choices]
            fields["multiSelect"] = bool(payload.get("multiSelect", False))
        return fields

    # -- every way a request ends ---------------------------------------------
    def resolve(self, req: tc.ConfirmRequest, decision: Any) -> bool:
        landed = super().resolve(req, decision)
        if landed:
            by = getattr(self._answered_by, "value", "") or "outside"
            self._announce_resolved(req, by)
        return landed

    def abort_all(self) -> list[tc.ConfirmRequest]:
        pending = super().abort_all()
        for req in pending:
            self._announce_resolved(req, "abort")
        return pending

    def _wait(self, req: tc.ConfirmRequest) -> Any:
        decision = super()._wait(req)
        if getattr(req, "expired", False):
            self._announce_resolved(req, "expired")
        reason = str(getattr(req, "protocol_reason", "") or "")
        if decision is False and reason:
            # Set on the ASKING thread, where the gate reads it.
            self.last_refusal_reason = reason
        return decision

    def _audit_dialog(self, req: tc.ConfirmRequest, *, by: str) -> None:
        # Runs on the asking thread inside the base _wait, before _wait above
        # copies the reason over: set it here so the audit record has it.
        reason = str(getattr(req, "protocol_reason", "") or "")
        if req.decision is False and reason:
            self.last_refusal_reason = reason
        super()._audit_dialog(req, by=by)

    def _announce_resolved(self, req: tc.ConfirmRequest, by: str) -> None:
        rid = self.request_id(req)
        if not rid:
            return
        decision = req.decision
        if req.kind == tc.CONFIRM:
            word = "approved" if decision is True else "denied"
        elif req.kind == tc.PLAN:
            word = ("approved" if isinstance(decision, dict)
                    and decision.get("approved") else "denied")
        else:
            word = "answered" if isinstance(decision, dict) and \
                decision.get("answers") else "dismissed"
        if by == "expired":
            word = "expired"
        elif getattr(req, "changed_to", ""):
            word = "changed"
        with self._ids_lock:
            if self._by_id.pop(rid, None) is None:
                return            # already announced
            self._emit("permission.resolved", requestId=rid, decision=word,
                       by=by)
        self._waiting_changed()

    def _waiting_changed(self) -> None:
        if self._on_waiting is not None:
            try:
                self._on_waiting(bool(self.pending()))
            except Exception:
                pass

    # -- answering --------------------------------------------------------------
    def answer(self, request_id: str, *, choice: str = "",
               answers: list[str] | None = None,
               reason: str = "", persist: bool = False,
               text: str = "", command: str = "") -> tuple[bool, str]:
        """Apply a window's answer the way the REPL applies a key.

        ``text`` is a line typed for a question (ask): resolved as the
        REPL's dialog resolves a typed line (``typed_answer``). Any other
        request refuses it -- typed text is never an approval.

        ``command`` is the user's own version of a shell command, and only
        with ``choice`` "c" (``TerminalConfirmBroker.change``): the command
        as asked is refused, this one approved once, exactly as written.

        Returns ``(landed, message)``: ``landed`` is False only when the
        request ended some other way in the meantime (answered from
        outside, aborted). Raises ``UnknownRequest`` for an id that is not
        open and ``InvalidAnswer`` for a key or answer the request did not
        offer -- a key that was never on screen must not mean anything,
        least of all "yes".
        """
        with self._ids_lock:
            req = self._by_id.get(str(request_id))
        if req is None or req.resolved:
            raise UnknownRequest(f"no open request {request_id!r}")
        self._answered_by.value = "client"
        try:
            return self._apply(req, choice=str(choice or ""),
                               answers=list(answers or []),
                               reason=str(reason or ""),
                               persist=bool(persist), text=str(text or ""),
                               command=str(command or ""))
        finally:
            self._answered_by.value = ""

    def _apply(self, req: tc.ConfirmRequest, *, choice: str,
               answers: list[str], reason: str,
               persist: bool = False, text: str = "",
               command: str = "") -> tuple[bool, str]:
        if command and choice != "c":
            raise InvalidAnswer("a changed command goes with the key 'c' only")
        if req.kind == tc.ASK:
            # Answered with labels; no answers at all is Esc: dismissed.
            if choice:
                raise InvalidAnswer("a question is answered with 'answers' "
                                    "(the chosen labels), not a 'choice'")
            payload = req.payload or {}
            if text.strip():
                # A typed line: an option when it names one, else the
                # answer as typed -- what the REPL's dialog and the
                # dashboard's chat input send for the same line.
                if answers:
                    raise InvalidAnswer("a question is answered with "
                                        "'answers' or 'text', not both")
                from delfin.agent.approval_answers import typed_answer
                return self.resolve(req, typed_answer(text, payload)), "answered"
            offered = [label for label, _ in _choices(payload)]
            if any(a not in offered for a in answers):
                raise InvalidAnswer("not one of the offered answers")
            if len(set(answers)) != len(answers):
                raise InvalidAnswer("an answer was given twice")
            if len(answers) > 1 and not payload.get("multiSelect"):
                raise InvalidAnswer("this question takes one answer")
            return self.resolve(req, {"answers": list(answers)}), \
                ("answered" if answers else "dismissed")

        if text:
            # Checked before the choice: a request carrying typed text is
            # refused whole, whatever key came with it.
            raise InvalidAnswer("an approval is answered with one of its "
                                "keys, never with typed text")
        options = tc.options_for(req, suggestion=self._suggestion(req))
        if choice not in {o.key for o in options}:
            raise InvalidAnswer(f"{choice!r} is not an option for this request")

        if req.kind == tc.PLAN:
            if choice == "d":
                return self.resolve(req, {"approved": True, "new_mode": "default"}), \
                    "approval → default"
            if choice == "e":
                return self.resolve(req, {"approved": True,
                                          "new_mode": "acceptEdits"}), \
                    "approval → acceptEdits"
            if not self.resolve(req, {"approved": False, "new_mode": "plan"}):
                return False, "the request ended before this answer arrived"
            if reason.strip() and self._on_feedback is not None:
                # The plan callback has no room for a reason (the engine
                # answers the model with "rejected" and the mode); what the
                # user said goes into the running turn instead, the way a
                # message steered into it does.
                try:
                    self._on_feedback(f"Plan rejected: {reason.strip()}")
                except Exception:
                    pass
            return True, "plan rejected"

        # CONFIRM: the keys of TerminalAgent._answer_request.
        if choice == "c":
            try:
                landed = self.change(req, command)
            except ValueError as exc:
                raise InvalidAnswer(str(exc)) from None
            if not landed:
                return False, "the request ended before this answer arrived"
            return True, "changed: your command runs next, once, exactly as typed"
        if choice == "y":
            return self.resolve(req, True), "approved"
        if choice == "n":
            req.protocol_reason = reason
            return self.resolve(req, False), "refused"
        if choice == "a":
            req.protocol_reason = reason
            landed = self.resolve(req, False)
            others = self.abort_all()
            if self._stop_turn is not None:
                try:
                    self._stop_turn()
                except Exception:
                    pass
            extra = f" and {len(others)} other request(s)" if others else ""
            return landed, f"aborted{extra}; the turn is ending"
        # The answer lands first, the side effect after -- and only when it
        # landed: an answer from outside that arrives in between must not be
        # followed by a pattern persisted on top of it.
        if choice == "e":
            if not self.resolve(req, True):
                return False, "the request ended before this answer arrived"
            _ok, message = self.accept_edits()
            return True, message
        if choice == "R":
            if not self.resolve(req, True):
                return False, "the request ended before this answer arrived"
            _ok, message = self.persist_read(tc.read_dir_for(req))
            return True, message
        if choice in ("A", "k"):
            # The REPL asks a second time before persisting: the pattern is
            # global and permanent. Over the protocol that second "yes" is the
            # explicit ``persist`` flag; without it the answer is "just this
            # once" -- never a pattern the user did not confirm.
            pattern = (self.exact_pattern(req.command) if choice == "A"
                       else str(getattr(req, "protocol_suggestion", "") or ""))
            if not self.resolve(req, True):
                return False, "the request ended before this answer arrived"
            if not persist:
                return True, "approved just this once (persist was not confirmed)"
            _ok, message = self.persist(pattern) if pattern else \
                (False, "nothing to persist")
            return True, message
        if choice == "?":
            raise InvalidAnswer(tc.render_help(req, options))
        raise InvalidAnswer(f"{choice!r} does nothing here")

    def _suggestion(self, req: tc.ConfirmRequest) -> str:
        if req.kind == tc.CONFIRM and req.is_shell and req.command:
            try:
                return self.suggest_pattern(req.command)
            except Exception:
                return ""
        return ""
