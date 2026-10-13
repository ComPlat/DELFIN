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
"""

from __future__ import annotations

import itertools
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


class ProtocolBroker(tc.TerminalConfirmBroker):
    """``TerminalConfirmBroker`` whose questions are answered by a window."""

    def __init__(self, emit: Emit, *, stop_turn: Callable[[], None] | None = None,
                 on_waiting: Callable[[bool], None] | None = None, **kwargs):
        super().__init__(timeout_s=0.0, **kwargs)
        self._emit = emit
        self._stop_turn = stop_turn
        self._on_waiting = on_waiting
        self._ids = itertools.count(1)
        self._by_id: dict[str, tc.ConfirmRequest] = {}
        self._ids_lock = threading.Lock()
        self._answered_by = threading.local()

    # -- announcing -----------------------------------------------------------
    def request_id(self, req: tc.ConfirmRequest) -> str:
        return str(getattr(req, "protocol_id", "") or "")

    def pending(self) -> list[tc.ConfirmRequest]:
        with self._ids_lock:
            return [r for r in self._by_id.values() if not r.resolved]

    def _enqueue(self, req: tc.ConfirmRequest) -> tc.ConfirmRequest:
        req = super()._enqueue(req)
        if req.resolved:              # refused on arrival by an abort
            return req
        rid = f"p{next(self._ids)}"
        req.protocol_id = rid
        with self._ids_lock:
            self._by_id[rid] = req
        self._emit("permission.asked", **self.describe(req))
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
        read_dir = tc.read_dir_for(req) if req.is_outside_read else ""
        if read_dir:
            fields["readDir"] = read_dir
        if req.kind == tc.ASK:
            payload = req.payload or {}
            fields["question"] = str(payload.get("question", "") or "")
            fields["choices"] = [str(o) for o in (payload.get("options") or [])][:9]
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

    def _announce_resolved(self, req: tc.ConfirmRequest, by: str) -> None:
        rid = self.request_id(req)
        if not rid:
            return
        with self._ids_lock:
            if self._by_id.pop(rid, None) is None:
                return            # already announced
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
        self._emit("permission.resolved", requestId=rid, decision=word, by=by)
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
               reason: str = "") -> tuple[bool, str]:
        """Apply a window's answer the way the REPL applies a key.

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
                               reason=str(reason or ""))
        finally:
            self._answered_by.value = ""

    def _apply(self, req: tc.ConfirmRequest, *, choice: str,
               answers: list[str], reason: str) -> tuple[bool, str]:
        if req.kind == tc.ASK:
            offered = [str(o) for o in ((req.payload or {}).get("options") or [])]
            picked = [a for a in answers if a in offered]
            if len(picked) != len(answers):
                raise InvalidAnswer("not one of the offered answers")
            return self.resolve(req, {"answers": picked}), "answered"

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
            req.protocol_reason = reason
            return self.resolve(req, {"approved": False, "new_mode": "plan"}), \
                "plan rejected"

        # CONFIRM: the keys of TerminalAgent._answer_request.
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
        # The side effect first, then the answer -- the REPL's order. The
        # request is approved either way, as in the REPL; the message says
        # whether the posture or pattern was actually saved.
        if choice == "e":
            _ok, message = self.accept_edits()
            return self.resolve(req, True), message
        if choice == "R":
            _ok, message = self.persist_read(tc.read_dir_for(req))
            return self.resolve(req, True), message
        if choice in ("A", "k"):
            pattern = (self.exact_pattern(req.command) if choice == "A"
                       else self._suggestion(req))
            _ok, message = self.persist(pattern) if pattern else \
                (False, "nothing to persist")
            return self.resolve(req, True), message
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
