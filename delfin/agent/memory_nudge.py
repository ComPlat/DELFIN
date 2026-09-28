"""Mid-session memory nudges — triggered by work, never by message count.

Package 7, phase 3. Hermes keeps its long-term memory small by curating
it DURING the session: at intervals the agent is prompted to ask itself
whether what it just learned is worth keeping, instead of everything
piling up until the session-end distill. DELFIN adopts the nudge with
its own hard rule carried over from a measured failure: the trigger is
WORK — tool calls executed and text produced since the last nudge —
never the number of messages. A message-counter trigger once fired at
15% context fill, when nothing of substance had happened yet.

The module is a pure decision function plus an in-memory counter state.
It writes nothing to disk (so it needs no plan-mode gate: a nudge is
text played into the context, nothing more) and the caller owns the
state object. ``agent.memory_nudge: false`` in the user settings
switches it off entirely.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

# Work thresholds since the last nudge. Both must be met: a hundred
# trivial tool calls with no output is churn, and a wall of text with
# no tool calls is a monologue — neither is the sustained
# tool-work-plus-findings pattern the nudge is for.
NUDGE_TOOL_CALLS = 12
NUDGE_CHARS = 4000

# Rough chars-per-token estimate (English/code mix, the usual 4:1).
_CHARS_PER_TOKEN = 4

_NUDGE_TEXT = (
    "[DELFIN note] An automatic system note, not the user speaking: a "
    "stretch of work has passed since anything was kept. Is anything "
    "here worth keeping — a preference you confirmed, a decision and "
    "its reason, a path or command that worked? If so, save it now "
    "with the remember tool (or /memorize at the end). The hard memory "
    "budget means the store stays curated, not accumulated: what is "
    "not kept is distilled away at compaction.")


@dataclass
class NudgeState:
    """In-memory counters the caller holds between nudge checks.

    ``tool_calls_since_nudge`` / ``chars_since_nudge`` accumulate work;
    ``last_nudge_turn`` / ``nudged_this_turn`` enforce the one-nudge-
    per-turn cap. Nothing is persisted: a crashed session loses the
    counters, at worst delaying one nudge to the next session.
    """

    tool_calls_since_nudge: int = 0
    chars_since_nudge: int = 0
    last_nudge_turn: Optional[int] = None
    nudged_this_turn: bool = False


def compose_nudge(store: Optional[Path] = None) -> str:
    """The nudge text, extended with the tidy hint when one is pending.

    Phase-3 fix on work/j2-memory-layers (2026-09-29): memory_tidy closes
    the duplication gap in the fact store, but it was reachable ONLY
    through a manual CLI command — measured 0 uses across 25 real
    sessions. The nudge fires exactly when near-duplicate facts may have
    just been written, so it carries the pointer when memory_tidy.hint
    has something to say. Read-only by contract: /tidy shows proposals
    and changes nothing until they are explicitly accepted.

    ``store`` is the memory store directory; None lets the caller decide
    later (then no tidy hint is composed). Never raises — a broken tidy
    pass must not take the nudge with it.
    """
    text = _NUDGE_TEXT
    if store is not None:
        try:
            from .memory_tidy import hint as _tidy_hint
            extra = _tidy_hint(Path(store))
            if extra:
                text = text + " " + extra.strip()
        except Exception:
            pass
    return text


def maybe_nudge(
    state: NudgeState,
    *,
    tool_calls: int = 0,
    chars: int = 0,
    turn: Optional[int] = None,
    settings: Optional[dict[str, Any]] = None,
) -> str:
    """Decide whether to nudge the agent to keep what it learned.

    Called from the turn loop with the work done since the last call
    (the caller feeds deltas or totals; the state accumulates whatever
    is passed). Returns the nudge text when the work thresholds are
    both met and no nudge was delivered this turn, else "". On a nudge
    the counters reset; on silence they keep accumulating, so partial
    work is never lost.

    ``settings`` is a user-settings dict (``{"agent": {...}}``); when
    absent the settings are read from disk like every other agent
    setting. Failures reading settings degrade to the nudge being ON —
    a broken settings file must not silently disable behaviour the user
    never asked to switch off.
    """
    state.tool_calls_since_nudge += max(0, int(tool_calls or 0))
    state.chars_since_nudge += max(0, int(chars or 0))
    if turn is not None and turn != state.last_nudge_turn:
        # A new turn opens: the per-turn cap resets with it.
        state.nudged_this_turn = False
        state.last_nudge_turn = turn

    if state.nudged_this_turn:
        return ""

    if not _enabled(settings):
        return ""

    if (state.tool_calls_since_nudge < NUDGE_TOOL_CALLS
            or state.chars_since_nudge < NUDGE_CHARS):
        return ""

    state.tool_calls_since_nudge = 0
    state.chars_since_nudge = 0
    state.nudged_this_turn = True
    return compose_nudge(store=None)


def estimated_tokens(chars: int) -> int:
    """Rough token estimate for a produced-text length."""
    return max(0, int(chars)) // _CHARS_PER_TOKEN


def _enabled(settings: Optional[dict[str, Any]]) -> bool:
    if settings is not None:
        try:
            agent = (settings or {}).get("agent") or {}
            return bool(agent.get("memory_nudge", True))
        except Exception:
            return True
    try:
        from delfin.user_settings import load_settings
        agent = (load_settings() or {}).get("agent") or {}
        return bool(agent.get("memory_nudge", True))
    except Exception:
        return True
