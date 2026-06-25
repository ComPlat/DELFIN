"""The subagent panel must be session-scoped.

The subagent telemetry file is global, so the dashboard panel previously kept
showing the previous session's (and old test) runs after 'New Session' (bug
2026-06-25: "neue session und ich seh immer noch explore … · 0 calls"). The
panel now filters telemetry by a per-session start timestamp.
"""

from pathlib import Path

_SRC = (Path(__file__).resolve().parent.parent / "delfin" / "dashboard"
        / "tab_agent.py").read_text(encoding="utf-8")


def test_panel_filters_by_session_start():
    # the refresh reads more than 3 and filters by the session-start stamp
    assert '_session_start_ts' in _SRC
    assert 'rec.get("ts", 0)' in _SRC
    assert 'read_telemetry(last_n=8)' in _SRC


def test_new_session_stamps_session_start():
    # New Session updates the stamp so old subagents stop showing
    assert 'state["_session_start_ts"] = _t_ns.time()' in _SRC
    # and it is initialised at build so the first session is scoped too
    assert 'state.setdefault("_session_start_ts"' in _SRC
