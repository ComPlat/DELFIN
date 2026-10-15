"""What one turn ended as, in fields a program can read.

Input: the gate's decisions as they happen (a refusal, a request for a
human answer, a file a write tool changed) and, after the turn, the
engine's own stop reason and error.

Output: ``stop_reason`` (one of :data:`STOP_REASONS`), ``denials``
(``[{tool, reason_code}]``), ``files_changed`` (workspace-relative paths)
and, for ``--output-format stream-json``, one ``permission_request`` /
``permission_denied`` event per decision.

Semantics:

* The decisions are reported by the code that makes them -- the
  containment recorder (:mod:`security_events`), the executor's central
  refusals, the approvals-directory broker and the write tools' change
  capture -- through :func:`denied`, :func:`requested` and
  :func:`file_changed`. Nothing here reads the refusal text the model
  sees. Rejected: parsing tool results for "denied"/"refused", which is
  the vocabulary a rubric scores and changes with every reworded message.
* ``reason_code`` is a containment kind from
  ``security_events.known_kinds()`` or one of :data:`LOCAL_REASON_CODES`,
  the refusals the panel deliberately does not record (a plan-mode or a
  role refusal is routine, not a containment event).
* ``stop_reason`` precedence: an error, then the engine's own ceiling
  (timeout, budget, max_turns), then a refusal in the turn's LAST tool
  round (``needs_approval`` when a human could lift it, ``refused``
  otherwise), then ``done``. A refusal the model worked around in a later
  round leaves the turn ``done``; it is still listed in ``denials``.
  Rejected: any refusal makes the turn ``refused`` -- the bash allow-list
  alone refuses in most long runs, and every one of them would exit
  non-zero after doing the job.
* Observers are process-wide, so a delegate's refusal inside the turn is
  reported too. Listener exceptions are swallowed: reporting must never
  change a gate decision.

Measured cost: one lock and a list walk per gate decision; no tool schema
or prompt text.
"""

from __future__ import annotations

import threading
from typing import Any, Callable, Iterable

__all__ = [
    "SCHEMA_VERSION", "STOP_REASONS", "EXIT_CODES", "REASON_CODES",
    "LOCAL_REASON_CODES", "APPROVABLE_REASON_CODES",
    "PERMISSION_REQUEST_KEYS", "PERMISSION_DENIED_KEYS", "TurnObserver",
    "denied", "requested", "file_changed", "classify",
]

#: Bumped when a key changes meaning or disappears; adding a key is not a
#: bump. Readers compare it before trusting the rest.
SCHEMA_VERSION = 1

#: stop_reason -> process exit code. 2 (usage) and 3 (the engine could not
#: be built) keep their meaning and are not here: no turn ran, so there is
#: no result object.
EXIT_CODES: dict[str, int] = {
    "done": 0,
    "error": 1,
    "needs_approval": 4,
    "refused": 5,
    "max_turns": 6,
    "budget": 7,
    "timeout": 8,
}
STOP_REASONS: tuple[str, ...] = tuple(EXIT_CODES)

#: Refusals with no containment kind, named after the MCP kinds that cover
#: the same decision on the other dispatch path.
LOCAL_REASON_CODES: frozenset[str] = frozenset({
    "role_denied",        # the session's role may not use the tool
    "session_denied",     # --allowed-tools / --disallowed-tools
    "plan_mode",          # plan mode refuses everything not read-only
    "blocked_by_hook",    # a PreToolUse hook blocked the call
    "not_allowlisted",    # bash outside the auto-allow list, nobody to ask
    "needs_approval",     # a dialog was required and none is bound
    "backend_denied",     # the CLI backend's own gate (no category given)
    "running_job",        # a file a running background job of the session uses
})

try:
    from .security_events import known_kinds as _known_kinds
    _KINDS = frozenset(_known_kinds())
except Exception:          # pragma: no cover - the module is in-tree
    _KINDS = frozenset()

REASON_CODES: frozenset[str] = _KINDS | LOCAL_REASON_CODES

#: A refusal a human could have lifted by answering: the gate wanted to
#: ask, or asked and nobody answered in time. ``denied_by_user`` is not
#: here -- somebody answered, and the answer was no.
APPROVABLE_REASON_CODES: frozenset[str] = frozenset({
    "needs_approval", "not_allowlisted", "approval_timeout", "self_mod",
    "calc_edit", "egress", "git_config_program", "mcp_side_effect",
    "classifier_ask", "shell_write_unasked", "outside_ws",
    # An ask rule (the administrator's or the user's) met with nobody to
    # ask: a person answering would have decided it.
    "managed_policy_ask", "permission_rule_ask",
})

#: Containment kinds recorded with ``blocked=True`` that are not a refusal
#: of the call in hand (a session-wide state, a lifted refusal, a grant).
_NOT_A_REFUSAL: frozenset[str] = frozenset({
    "isolation", "isolation_missing", "key_exported", "settings_updated",
    "approval_from_outside", "read_grant", "refusal_lifted", "auto_verify",
    "auto_verify_exhausted", "test_tamper", "read_only_shape",
    "outside_read_bypass", "classifier_error",
})

PERMISSION_REQUEST_KEYS: frozenset[str] = frozenset({
    "type", "request_id", "tool", "shown", "preview", "answered_by",
    "mode", "approvals_dir"})
PERMISSION_DENIED_KEYS: frozenset[str] = frozenset({
    "type", "tool", "reason_code", "needs_approval", "text"})

_PREVIEW_MAX = 4000
_TEXT_MAX = 600

# Engine stop reasons (api_client message_delta, engine run-budget gate)
# that are not "done".
_ENGINE_STOP: dict[str, str] = {
    "max_tool_rounds": "max_turns",
    "cost_cap": "budget",
    "max_turn_output": "budget",
    "run_budget_usd": "budget",
    "run_budget_turns": "budget",
    "run_budget_s": "timeout",
    "error": "error",
    "stopped": "error",
    "max_context": "error",
    "no_progress_loop": "error",
    "consecutive_identical_errors": "error",
}

_LOCK = threading.Lock()
_OBSERVERS: list["TurnObserver"] = []


def _redact(text: str) -> str:
    try:
        from .cli_approve import redact_preview
        return redact_preview(text)
    except Exception:
        return "<redacted: preview unavailable>"


def _observers() -> list["TurnObserver"]:
    with _LOCK:
        return list(_OBSERVERS)


def _bare(name: str) -> str:
    name = str(name or "")
    return name.rsplit("__", 1)[-1] if name.startswith("mcp__") else name


def denied(tool: str, reason_code: str, text: str = "") -> None:
    """The gate refused ``tool``. Never raises."""
    if reason_code in _NOT_A_REFUSAL:
        return
    for obs in _observers():
        try:
            obs._on_denied(_bare(tool), str(reason_code or ""), str(text or ""))
        except Exception:
            pass


def requested(*, request_id: str, tool: str, args: dict | None,
              preview: str, approvals_dir: bool) -> None:
    """A call is waiting for a human answer. Never raises."""
    for obs in _observers():
        try:
            obs._on_requested(str(request_id or ""), _bare(tool),
                              dict(args or {}), str(preview or ""),
                              bool(approvals_dir))
        except Exception:
            pass


def file_changed(path: str) -> None:
    """A write tool changed ``path`` (absolute or workspace-relative)."""
    for obs in _observers():
        try:
            obs._on_file_changed(str(path or ""))
        except Exception:
            pass


def classify(*, error: str, engine_stop_reason: str,
             denials: Iterable[dict], ended_on_denial: bool) -> str:
    """The turn's ``stop_reason``; see the module docstring for the order."""
    if error:
        return "error"
    mapped = _ENGINE_STOP.get(str(engine_stop_reason or ""), "")
    if mapped:
        return mapped
    if ended_on_denial:
        last = list(denials)
        if last and last[-1].get("reason_code") in APPROVABLE_REASON_CODES:
            return "needs_approval"
        return "refused"
    return "done"


class TurnObserver:
    """Collect one turn's gate decisions; optionally stream them.

    ``emit`` receives the ``permission_request`` / ``permission_denied``
    event dicts; ``mode`` returns the permission mode at the time of a
    request; ``workspace`` makes ``files_changed`` relative.
    """

    def __init__(self, *, emit: Callable[[dict], None] | None = None,
                 mode: Callable[[], str] | None = None,
                 workspace: Any = None) -> None:
        self._emit = emit
        self._mode = mode
        self._workspace = workspace
        self._lock = threading.Lock()
        self.denials: list[dict] = []
        self.requests: list[dict] = []
        self._files: list[str] = []
        self._denied_since_tool_use = False
        self._current_tool = ""

    # -- lifetime --------------------------------------------------------
    def __enter__(self) -> "TurnObserver":
        with _LOCK:
            _OBSERVERS.append(self)
        return self

    def __exit__(self, *exc) -> None:
        with _LOCK:
            if self in _OBSERVERS:
                _OBSERVERS.remove(self)

    # -- what the turn reports ------------------------------------------
    def saw_tool_use(self, name: str = "") -> None:
        """A new tool call started; earlier refusals were not the last.

        ``name`` stands in for the tool of a refusal reported without one
        (a shared read gate does not know which tool called it).
        """
        with self._lock:
            self._denied_since_tool_use = False
            self._current_tool = _bare(name)

    @property
    def ended_on_denial(self) -> bool:
        """A refusal arrived after the turn's last tool call started."""
        return self._denied_since_tool_use

    @property
    def files_changed(self) -> list[str]:
        with self._lock:
            return list(self._files)

    def summary(self, *, error: str, engine_stop_reason: str) -> dict:
        """The four result fields this module owns."""
        with self._lock:
            denials = [{"tool": d["tool"], "reason_code": d["reason_code"]}
                       for d in self.denials]
        return {
            "schema_version": SCHEMA_VERSION,
            "stop_reason": classify(
                error=error, engine_stop_reason=engine_stop_reason,
                denials=denials, ended_on_denial=self.ended_on_denial),
            "denials": denials,
            "files_changed": self.files_changed,
        }

    # -- the gate's side -------------------------------------------------
    def _send(self, event: dict) -> None:
        if self._emit is None:
            return
        try:
            self._emit(event)
        except Exception:
            pass

    def _on_denied(self, tool: str, code: str, text: str) -> None:
        with self._lock:
            tool = tool or self._current_tool
            self.denials.append({"tool": tool, "reason_code": code})
            self._denied_since_tool_use = True
        self._send({
            "type": "permission_denied",
            "tool": tool,
            "reason_code": code,
            "needs_approval": code in APPROVABLE_REASON_CODES,
            "text": _redact(text)[:_TEXT_MAX],
        })

    def _on_requested(self, request_id: str, tool: str, args: dict,
                      preview: str, approvals_dir: bool) -> None:
        shown = str(args.get("command") or args.get("path") or "")
        try:
            mode = str(self._mode() if self._mode else "")
        except Exception:
            mode = ""
        event = {
            "type": "permission_request",
            "request_id": request_id,
            "tool": tool,
            "shown": _redact(shown)[:_TEXT_MAX],
            "preview": _redact(preview)[:_PREVIEW_MAX],
            # The answer is a person's. A driver reading this is told who
            # decides, never asked to decide.
            "answered_by": "human",
            "mode": mode,
            "approvals_dir": approvals_dir,
        }
        with self._lock:
            self.requests.append(event)
        self._send(event)

    def _on_file_changed(self, path: str) -> None:
        rel = path
        if self._workspace:
            try:
                from pathlib import Path
                rel = str(Path(path).resolve().relative_to(
                    Path(self._workspace).resolve()))
            except Exception:
                rel = path
        with self._lock:
            if rel and rel not in self._files:
                self._files.append(rel)
