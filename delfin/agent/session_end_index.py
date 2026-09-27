"""The session-end re-index step (Paket 5, Phase 4).

Paket 2 (nacht-s12) owns the shared session-end function; this module is
the step DELFIN runs there: index the just-finished session so it is
searchable immediately. Kept as its own function so the session-end
call site is one line, whatever shape Paket 2's function takes.

Contract:
- reindex_finished_session(session_id) never raises -- a broken index
  must not break session end (proven by test, including a corrupt index
  file).
- It reads only the session's own archived sources and writes only the
  index.
"""

from __future__ import annotations


def reindex_finished_session(session_id: str) -> bool:
    """Index the session that just ended. Never raises.

    Returns True when something was indexed, False otherwise (unknown
    session, nothing changed, or an index failure -- all swallowed: the
    session is over either way, and the next index pass catches up).
    """
    try:
        if not session_id:
            return False
        from . import session_index
        return bool(session_index.index_session(session_id))
    except Exception:
        return False
