"""The Background panel says when it could not read a list, instead of
showing nothing running.

Review 2026-09-16: every group in background_view.collect swallowed its
read error, and the panel's own exception handler blanked it -- a failed
read looked exactly like an idle agent.
"""
import inspect
from pathlib import Path

from delfin.agent import background_view as BV


def test_a_failed_group_becomes_an_error_row(tmp_path, monkeypatch):
    import delfin.agent.bash_jobs as bj

    def boom(*a, **k):
        raise OSError("registry unreadable")
    monkeypatch.setattr(bj, "live_jobs", boom)
    view = BV.collect(str(tmp_path))
    assert view["errors"] and view["errors"][0]["group"] == "shells"
    assert "registry unreadable" in view["errors"][0]["error"]
    rows = [r for r in BV.rows(view) if r["group"] == "errors"]
    assert rows and rows[0]["kind"].startswith("⚠")
    assert "shells" in rows[0]["label"]
    assert rows[0]["tip"] == BV._STOP_TIPS["errors"]


def test_a_clean_read_has_no_error_row(tmp_path):
    view = BV.collect(str(tmp_path))
    assert view["errors"] == []
    assert not [r for r in BV.rows(view) if r["group"] == "errors"]


def test_the_panel_says_it_could_not_read():
    src = Path(inspect.getfile(__import__("delfin.dashboard.tab_agent", fromlist=["x"]))).read_text()
    i = src.index("def _refresh_subagent_panel(")
    body = src[i:i + 5000]
    assert "the list could not be read" in body
    assert 'if _row["group"] == "errors":' in body
