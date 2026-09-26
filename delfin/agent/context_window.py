"""The effective context window: ONE source for engine and client.

The cap logic (``DELFIN_CONTEXT_WINDOW_CAP`` / ``agent.context_window_cap``)
lived in engine.py only. The client's in-turn elision budgeted against
the RAW model window and never saw the cap — measured 2026-09-26: a
40 000-token cap on a 131 000-token model gave a 235 800-char elision
budget, three times the whole capped window, so a capped turn overflowed
before the elision fired. Both sides now read :func:`capped_window`;
engine.py keeps a thin forwarder for its existing callers.
"""

from __future__ import annotations

import os


def capped_window(raw_window: int) -> int:
    """The model's context window, capped by ``agent.context_window_cap``
    or ``DELFIN_CONTEXT_WINDOW_CAP`` when either is set (tokens, > 0).

    Compaction fires at a share of the window. On a 131k model a whole
    supervised session never reached it (measured 2026-09-26: four long
    sessions, zero compactions), so whether a session survives several
    compactions could not be observed at all -- and a user on a small
    local model gets no way to make DELFIN compact earlier. A cap only
    ever lowers the window; it never claims more than the model has.
    """
    cap = 0
    try:
        cap = int(os.environ.get("DELFIN_CONTEXT_WINDOW_CAP", "") or 0)
    except ValueError:
        cap = 0
    if cap <= 0:
        try:
            from delfin import user_settings
            ag = (user_settings.load_settings().get("agent", {}) or {})
            cap = int(ag.get("context_window_cap", 0) or 0)
        except Exception:
            cap = 0
    return min(raw_window, cap) if cap > 0 else raw_window
