"""SWV phase 2 (red control): the dashboard watch pane's pure logic.

The dashboard tab keeps SW2a's pattern (tab_swarm.py docstring, lines 9-11):
testable logic lives in module-level pure helpers that run without
ipywidgets; ``create_watch_tab`` is the thin ipywidgets half. This file
pins the watch formatters (table rows, member feed) built on the observe
API, plus the read-only guarantee: the watch can render a member but never
answer an approval or type into a session.

RED before the fix: these helpers do not exist yet in tab_swarm.py.
"""

from __future__ import annotations

import importlib

from delfin.agent.swarm import observe as _observe

# --- helpers reused from the phase-1 contract (real shape) ---------------
# A member row exactly as swarm_overview produces it (observe.py) and a feed
# entry exactly as member_feed produces it --- real writers, not hand-typed.
_REAL_ROW = {
    "package": "swv", "role": "builder", "name": "swv-builder",
    "model": "<model>", "state": "NEEDS APPROVAL",
    "commits_on_branch": 2, "last_tool_call": {"tool": "write_file",
                                               "arg": "delfin/x.py"},
    "last_commit_subject": "fix the thing", "last_activity_age_s": 12.0,
    "pending_approvals_count": 1,
}

_REAL_FEED = [
    {"at": "2026-06-01T10:00:05", "kind": "tool",
     "text": "run_tests tests/t.py"},
    {"at": "2026-06-01T10:00:10", "kind": "refusal",
     "text": "write_file delfin/x.py └─ outside scope"},
]


def _tab():
    return importlib.import_module("delfin.dashboard.tab_swarm")


def test_watch_table_one_row_per_member_escaped():
    tab = _tab()
    rows = tab.watched_rows_html([_REAL_ROW])
    assert isinstance(rows, list) and len(rows) == 1
    html = rows[0]
    # the field values are HTML-escaped (hostile input stays inert)
    assert "&lt;model&gt;" in html
    assert "swv-builder" in html
    assert "NEEDS APPROVAL" in html


def test_watch_table_highlights_needs_approval():
    tab = _tab()
    html = tab.watched_rows_html([_REAL_ROW])[0]
    # a member awaiting approval is visually flagged, not hidden
    assert "approval" in html.lower()


def test_watch_table_is_read_only_no_controls():
    tab = _tab()
    html = "".join(tab.watched_rows_html([_REAL_ROW]))
    # no approve button, no input box, no answer link -- observing never acts
    assert "<button" not in html and "<input" not in html
    assert "approve(" not in html and "answer(" not in html


def test_watch_feed_newest_last_with_refusal_highlighted():
    tab = _tab()
    entries = tab.watched_feed_html(_REAL_FEED)
    assert len(entries) == 2
    # newest last (the refusal is the newest)
    assert "outside scope" in entries[-1]
    # the refusal is visually flagged
    assert "refusal" in entries[-1].lower() or "denied" in entries[-1].lower()


def test_watch_feed_is_read_only_no_controls():
    tab = _tab()
    html = "".join(tab.watched_feed_html(_REAL_FEED))
    assert "<button" not in html and "<input" not in html


def test_watch_wrappers_route_to_observe_api():
    """The tab's watch calls the SAME observe functions the TUI and an SW3
    lead call -- one data path, not a second implementation."""
    tab = _tab()
    assert tab.watch_overview is _observe.swarm_overview
    assert tab.watch_feed is _observe.member_feed
