"""The subagent panel must be session-scoped.

The subagent telemetry file is global, so the dashboard panel previously kept
showing the previous session's (and old test) runs after 'New Session' (bug
2026-06-25: "neue session und ich seh immer noch explore … · 0 calls"). The
panel then filtered telemetry by a per-session start timestamp.

Since report 20260915-132613 it shows no finished run at all -- finished
reports live in the chat -- and a running sub-agent only when this session
started it (background_view.collect checks the owner).
"""

import inspect
from pathlib import Path
from types import SimpleNamespace

_SRC = (Path(__file__).resolve().parent.parent / "delfin" / "dashboard"
        / "tab_agent.py").read_text(encoding="utf-8")


def test_panel_shows_no_finished_run_and_no_other_sessions_agent(monkeypatch, tmp_path):
    from delfin.agent import background_view as BV
    from delfin.agent import bash_jobs as BJ
    from delfin.agent import job_monitor as JM
    from delfin.agent import scheduler as SCH
    from delfin.agent import subagents as SA

    monkeypatch.setattr(JM, "_AGENT_WATCH_INDEX_PATH", tmp_path / "index.json")
    monkeypatch.setattr(BJ, "live_jobs", lambda: [])
    monkeypatch.setattr(SCH, "get_scheduler",
                        lambda: SimpleNamespace(list_entries=lambda: []))
    monkeypatch.setattr(SA, "read_running", lambda **_: {
        "old": {"type": "explore", "description": "previous session",
                "started_at": 1.0, "owner_pid": 1}})
    monkeypatch.setattr(SA, "_entry_owned_by_us", lambda entry: False)
    assert BV.render_html(BV.collect(tmp_path)) == ""
    assert "read_telemetry" not in inspect.getsource(BV)


def test_new_session_stamps_session_start():
    # New Session updates the stamp so old subagents stop showing
    assert 'state["_session_start_ts"] = _t_ns.time()' in _SRC
    # and it is initialised at build so the first session is scoped too
    assert 'state.setdefault("_session_start_ts"' in _SRC
