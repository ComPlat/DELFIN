"""SA1 pilot: a measured, narrow read-only explore helper.

This module is a MEASUREMENT HARNESS, not a new capability. DELFIN already
has sub-agents (``delfin/agent/subagents.py``, preset ``explore`` at
:1260-1276: mode "plan", read-only tool set). The wave task forbids them by
default; this pilot measures whether ONE narrow use pays off, and the
helper becomes a recommendation (a skill) only if the measurement says so.

Phase 1 — the guard. The pilot setting allows ONLY:
- the ``explore`` sub-agent type (read-only filesystem + grep + web search),
- at most one sub-agent running per session,
- no nesting (a delegate may not start another delegate),
- no background runs (the harness measures one blocking run per arm).

Everything else is refused with the reason stated — a refusal names a
ROUTE, so every refusal here says what is not part of the pilot setting.

The real enforcement point for sub-agent calls is the engine's
``subagent_runner`` binding (engine.py:961-992, a protected file). This
guard lives at the wrapper the pilot harness calls; the wiring question
(engine hook vs. wrapper-only) is tracked in the package's PR notes.

Phase 2 — the measurement harness (``run_pilot``) runs the same small set
of read-heavy tasks twice: once plain, once allowed to call the explore
sub-agent, and compares tokens, wall time, tool calls and correctness.
"""
from __future__ import annotations

from types import SimpleNamespace
from typing import Any

from . import subagents as _sa

# The one type the pilot setting allows.
PILOT_TYPE = "explore"

#: Background runs (reserve_running) are long-lived: a reserved slot counts
#: as running until the delegate is reaped, and is exactly where two
#: delegates could overlap.
BOTH_CHANNELS = "pilot setting allows at most one explore sub-agent running per session"

#: The real runner, bound at import so tests can monkeypatch it and so the
#: refusal paths below demonstrably never reach it.
_real_run_subagent = _sa.run_subagent


def check_delegate(subagent_type: str, *, slot: Any = None,
                   depth: int = 0, isolation: str = "") -> tuple[bool, str]:
    """Should this delegation be allowed under the pilot setting?

    Returns ``(ok, reason)``. ``reason`` is "" when allowed; otherwise a
    sentence stating what the pilot setting does not allow. Pure: no I/O,
    no state of its own — the one-running slot is passed in.
    """
    if subagent_type != PILOT_TYPE:
        return False, (
            f"pilot setting allows only the {PILOT_TYPE!r} sub-agent type; "
            f"{subagent_type!r} is refused")
    if depth > 0:
        return False, (
            "pilot setting allows no nesting: a delegate may not start "
            "another delegate")
    if isolation:
        return False, (
            f"pilot setting allows no isolation={isolation!r}: the pilot "
            "runs read-only in the harness's own workspace")
    if slot is not None and not _slot_free(slot):
        return False, (
            f"{BOTH_CHANNELS}; one is already running")
    return True, ""


# ---------------------------------------------------------------------------
# The one-running-per-session slot
# ---------------------------------------------------------------------------

def pilot_slot() -> Any:
    """A fresh slot object for one pilot run (context manager)."""
    return _PilotSlot()


def _slot_free(slot: Any) -> bool:
    return not getattr(slot, "claimed", False)


def claim_slot(slot: Any) -> None:
    """Mark the slot as taken (one sub-agent running)."""
    slot.claimed = True


def reap_slot(slot: Any) -> None:
    """Release the slot (the delegate finished or was reaped)."""
    slot.claimed = False


class _PilotSlot:
    """Claimed for the duration of one sub-agent run; released on exit."""

    def __init__(self) -> None:
        self.claimed = False

    def __enter__(self) -> "_PilotSlot":
        return self

    def __exit__(self, *exc) -> None:
        self.claimed = False


# ---------------------------------------------------------------------------
# The guarded entry the harness calls
# ---------------------------------------------------------------------------

def run_pilot_subagent(*, subagent_type: str, description: str, prompt: str,
                       parent_client: Any, parent_perms: Any,
                       model: str = "", inherit_context: bool = True,
                       depth: int = 0, _pilot_slot: Any = None,
                       **kwargs: Any) -> Any:
    """``run_subagent`` under the pilot setting.

    Every refusal returns a ``SubagentResult`` with the reason in
    ``error`` — the same shape the real runner returns for its own
    refusals — and the real runner is never called.
    """
    slot = _pilot_slot if _pilot_slot is not None else _PilotSlot()
    ok, why = check_delegate(subagent_type, slot=slot, depth=depth)
    if not ok:
        return SimpleNamespace(subagent_type=subagent_type,
                               description=description,
                               final_text="", error=why)
    claim_slot(slot)
    try:
        return _real_run_subagent(
            subagent_type=subagent_type, description=description,
            prompt=prompt, parent_client=parent_client,
            parent_perms=parent_perms, model=model,
            inherit_context=inherit_context, **kwargs)
    finally:
        slot.claimed = False
