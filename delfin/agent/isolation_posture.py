"""Whether a shell command is held by a filesystem layer, and what follows.

Input: the host's three layers (bubblewrap, Landlock, Seatbelt; the probes
in ``api_client``), the shell-isolation mode (``--isolate``/``--no-isolate``,
``agent.bash_isolation``; ``api_client._bash_isolation_setting``) and the
session's permissions. Output: an ``IsolationStatus`` for the banner,
``/status`` and the turn status line, and ``apply_at_start``, the one
start check the interactive (``cmd_chat``) and one-shot (``cmd_run``)
paths both call.

Semantics. Without a layer, what stops a command is a reading of its
text; a command that assembles its words at run time is not read. So:

* ``bypassPermissions`` asks nothing per command, so it needs a layer:
  without one the session starts in ``default`` with a notice and the
  permissions marked ``unisolated_fallback``, attended or not. The run
  continues, and a command outside the auto-allow list is asked, or with
  nobody to ask refused (``api_client._unisolated_fallback_verdict``; no
  classifier verdict lets it through, since without a layer nothing but a
  reading of its text would hold it).
* ``acceptEdits`` starts as asked, with the notice: file writes stay inside
  the workspace through the path gates, and commands are still asked.
* ``default`` and ``plan`` start as asked: every command is asked anyway.
* A session that promises confinement (locked scope, ``--isolate``) on a
  host with no layer refuses its commands rather than running them; that
  is ``refuses`` here, needs no start check, and is asked per command in
  the bash gate (``[NO ISOLATION]``).

Rejected: exiting with an error for an unattended bypass start without a
layer -- long unattended runs on such hosts would stop working; running
them in ``default`` keeps reads, workspace writes and the test runner.
Rejected (2026-11-13): keeping bypass on a person's yes. Bypass is the
mode whose damage nobody reviews, so it runs only where a layer confines
what a command can change; one confirmation at the start does not confine
anything. A host without a layer runs the modes that ask instead.
Rejected: one confirmation per mode switch -- shift+tab reaches bypass
only where a layer holds each command (posture_cycle.bypass_refusal), so
the unisolated case is decided here, at start, alone.

Cost: the three probes are cached per process; one settings read.
"""

from __future__ import annotations

from dataclasses import dataclass

_RISKY_MODE = "bypassPermissions"
_FALLBACK_MODE = "default"


@dataclass(frozen=True)
class IsolationStatus:
    """isolated: a layer holds each command. layer: its name, or "".
    refuses: commands are refused (confinement promised, no layer).
    reason: why not isolated ("" when isolated)."""

    isolated: bool
    layer: str = ""
    refuses: bool = False
    reason: str = ""

    def line(self) -> str:
        """One line, as the banner and /status print it after the label."""
        if self.isolated:
            return (f"on ({self.layer}) — a command the agent runs writes "
                    "only inside the workspace (--no-isolate)")
        if self.refuses:
            return ("NOT isolated, and commands are refused — "
                    f"{self.reason}")
        return (f"NOT isolated — {self.reason}, so a command the agent "
                "runs can write outside the workspace")

    def short(self) -> str:
        """A word for the turn status line: "" when isolated."""
        if self.isolated:
            return ""
        return "commands refused" if self.refuses else "unisolated"


def _layer() -> str:
    try:
        from .doctor import _isolation_mechanism
        return _isolation_mechanism()
    except Exception:
        return ""


def isolation_status(perms=None) -> IsolationStatus:
    """The isolation a shell command of this session gets. Never raises."""
    try:
        from . import api_client as A
        mode, source = A._bash_isolation_setting()
    except Exception:
        return IsolationStatus(False, reason="the isolation setting could "
                                             "not be read")
    try:
        if A._isolation_refusal_reason(perms, mode):
            return IsolationStatus(
                False, refuses=True,
                reason=("confinement is promised (locked folder or "
                        "--isolate) and neither bubblewrap, Landlock nor "
                        "Seatbelt works on this host"))
    except Exception:
        pass
    if mode == "off":
        why = ("switched off by --no-isolate" if source == "--no-isolate"
               else 'switched off by agent.bash_isolation = "off"')
        return IsolationStatus(False, reason=why)
    if mode not in ("auto", "bwrap"):
        return IsolationStatus(
            False, reason=f'agent.bash_isolation = "{mode}" names no layer')
    layer = _layer()
    if layer:
        return IsolationStatus(True, layer=layer)
    return IsolationStatus(
        False, reason=("neither bubblewrap, Landlock nor Seatbelt works on "
                       "this host"))


def _record(detail: str, *, blocked: bool) -> None:
    try:
        from .api_client import _record_security_event
        _record_security_event("isolation_missing", "session", detail,
                               blocked=blocked)
    except Exception:
        pass


def apply_at_start(engine, *, can_ask: bool) -> list[str]:
    """Settle the starting permission mode against the isolation status.

    Input: the constructed engine (its ``kit_permissions.mode`` is the mode
    asked for); ``can_ask`` whether a person can answer per command. Output:
    notice lines for the banner or stderr, empty where nothing needs saying.
    Side effect: a bypass start without a layer is lowered to ``default``
    and the permissions are marked ``unisolated_fallback``. Called after an
    engine is built and after a live switch of its mode. A backend without
    a permissions object has no mode and is left alone (the CLI backend's
    commands are caged by its own gate, ``cli_gate``).
    """
    perms = getattr(engine, "kit_permissions", None)
    if perms is None:
        return []
    mode = str(getattr(perms, "mode", "") or "")
    status = isolation_status(perms)
    if status.isolated or status.refuses:
        return []
    head = f"isolation  {status.line()}"
    if mode == "acceptEdits":
        return [head + "; file writes stay inside the workspace, commands "
                       "are still asked"]
    if mode != _RISKY_MODE:
        return []
    _lower(perms)
    _record("bypassPermissions needs a layer that holds each command; "
            f"started in {_FALLBACK_MODE} ({status.reason})", blocked=True)
    tail = ("any other command is REFUSED, because nobody can be asked"
            if not can_ask else "any other command is asked")
    return [head + "; bypassPermissions runs only where a layer holds each "
                   f"command, so this session runs in {_FALLBACK_MODE}: "
                   "reads, workspace writes, the test runner and the "
                   f"auto-allowed commands go through; {tail} (use a host "
                   "with bubblewrap, Landlock or Seatbelt, and isolation "
                   "on, for a run that asks nothing)"]


def _lower(perms) -> None:
    """Default mode and the fallback mark, on the permissions object.

    Not through ``engine.set_kit_permission_mode``: that one tells a
    running turn the user switched mid-turn, and before the first turn
    nobody switched."""
    try:
        perms.mode = _FALLBACK_MODE
        perms.unisolated_fallback = True
    except Exception:
        pass


__all__ = ["IsolationStatus", "isolation_status", "apply_at_start"]
