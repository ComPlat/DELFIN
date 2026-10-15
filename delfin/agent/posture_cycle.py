"""The cycle key's ladder (shift+tab in the REPL and in the window).

Input: the engine (its ``kit_permissions``), the mode it is in, and
whether bypass was confirmed in this session. Output: the next mode, or a
confirmation to ask for.

Order: default -> acceptEdits -> plan -> bypassPermissions -> default.

``bypassPermissions`` is on the ladder only when

* the administrator's managed policy allows it (``mode_refusal``), and
* a filesystem layer holds every command (``isolation_posture
  .isolation_status``: isolated, not "refused" or "unisolated"), and the
  session is not running in the unisolated fallback: bypass asks for no
  command, and without the cage a command that assembles its words at run
  time is not read by anything;

otherwise the ladder steps over it. The first step onto it in a session is
not taken by the key: the caller asks once (y in the REPL, y or a click
on the window's confirm row), and only that answer switches. An unanswered
or refused offer leaves the mode where it was. After a yes the key passes
through bypass without asking again in that session.

The switch itself is the engine's ``set_kit_permission_mode`` (which also
refuses a mode the policy disallows) and is recorded as a security event.
Only a key or a click reaches this module: the service's
``session.setMode`` is an ACTING method, refused to every process an agent
session started.

Rejected: a repeated shift+tab as the confirmation -- a key held down or
pressed twice would land on bypass by rhythm, not by decision. Rejected:
keeping bypass off the ladder entirely (the earlier rule) -- the user asked
for parity with the start flag's reach, with the cage as the condition.

Cost: one policy read and the cached isolation probes per key press.
"""

from __future__ import annotations

from typing import Optional

LADDER = ("default", "acceptEdits", "plan", "bypassPermissions")
BYPASS = "bypassPermissions"

#: A y sooner than this after the offer is typing, not an answer
#: (the window's ANSWER_PAUSE_MS, footer/approval.tsx).
_DELIBERATE_S = 0.35

#: What the confirmation says, in both front ends.
OFFER = ("bypassPermissions asks for no command and no write from here on "
         "(isolation holds each command)")


def bypass_refusal(engine) -> Optional[str]:
    """Why the ladder steps over bypass for this session, or None when it
    may stop there. Never raises: an error is a refusal."""
    perms = getattr(engine, "kit_permissions", None)
    if perms is None:
        return "this backend carries no permission gate"
    try:
        from . import managed_policy as _mp
        refused = _mp.mode_refusal(BYPASS)
    except Exception:
        return "the managed policy could not be read"
    if refused:
        return refused
    if getattr(perms, "unisolated_fallback", False):
        return "this session runs in the unisolated fallback"
    try:
        from .isolation_posture import isolation_status
        status = isolation_status(perms)
    except Exception:
        return "the isolation status could not be read"
    if not status.isolated:
        return f"isolation is not active ({status.reason or 'no layer'})"
    return None


def ladder(engine) -> list[str]:
    """The modes the key steps through, in order."""
    try:
        from . import managed_policy as _mp
        order = [m for m in LADDER if m == BYPASS or _mp.mode_refusal(m) is None]
    except Exception:
        order = list(LADDER)
    if bypass_refusal(engine) is not None:
        order = [m for m in order if m != BYPASS]
    return order or ["default"]


def next_mode(current: str, order: list[str]) -> str:
    """The mode after ``current``; from one off the ladder, default."""
    if current not in order:
        return "default" if "default" in order else order[0]
    return order[(order.index(current) + 1) % len(order)]


def step(engine, current: str, *, confirmed: bool) -> tuple[str, bool]:
    """``(mode, ask)``: the mode the key moves to, and whether that move
    is bypass and must be confirmed first (``ask`` True: nothing has been
    switched)."""
    nxt = next_mode(current, ladder(engine))
    return nxt, (nxt == BYPASS and not confirmed)


def enter_bypass(engine, *, via: str) -> Optional[str]:
    """Switch to bypass after the person confirmed it: None on success,
    else why not (re-checked here: the cage may have gone since the
    offer). Recorded as a security event."""
    why = bypass_refusal(engine)
    if why:
        return why
    setter = getattr(engine, "set_kit_permission_mode", None)
    try:
        ok = bool(callable(setter) and setter(BYPASS))
    except Exception as exc:                    # noqa: BLE001
        return f"could not switch: {exc}"
    if not ok:
        return "the engine refused the switch"
    try:
        from .api_client import _record_security_event
        _record_security_event("permission_mode", "session",
                               f"bypassPermissions entered by the user ({via})",
                               blocked=False)
    except Exception:
        pass
    return None


__all__ = ["LADDER", "BYPASS", "OFFER", "bypass_refusal", "ladder", "next_mode",
           "step", "enter_bypass"]
