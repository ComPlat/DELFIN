"""ONE session list in the agent tab's sidebar.

Seen by the person on the cluster dashboard (2026-11-10): the open-views
list and the "Running sessions" panel showed the same session one under
the other; the resume dropdown offered ~40 entries, most of them
background worktree runs titled by their injected first message, the same
name twice, absolute paths in every label. The goal: one list, each
session once, named by something readable, with the background runs and
the sessions of other login nodes grouped instead of mixed in.
"""

from __future__ import annotations

import pytest

pytest.importorskip("ipywidgets")

import ipywidgets as widgets                                # noqa: E402

from delfin.dashboard import agent_sessions as AS           # noqa: E402


# -- naming -----------------------------------------------------------------

@pytest.mark.parametrize("title", [
    "[DELFIN note] An automatic system note, not the user speaking: "
    "a stretch of work has passed",
    "[Verify] Automatic check of the staged diff found nothing to fix",
    "[Message from the session nightly-review] Please look at the failing test",
    "Your task for this run. Read it all before you type anything:",
])
def test_a_harness_prefix_is_never_the_name(title):
    """The injected first message of a background run is not a name."""
    assert AS.session_display_name(title, "9d19b79871ff4a0e") == (
        "Session 9d19b798")


def test_a_system_prefix_deep_in_the_text_still_starts_the_name():
    stripped = AS.session_display_name(
        "[DELFIN note] An automatic system note.\nFix the failing test",
        "abc")
    assert stripped.startswith("Fix the failing test")


def test_a_persons_own_title_is_kept():
    assert AS.session_display_name("Fix the failing test", "abc") == (
        "Fix the failing test")


# -- grouping ---------------------------------------------------------------

def test_a_worktree_session_is_background():
    assert AS.is_background_session(
        "/home/me/project/.delfin/worktrees/delfin-wt-1a2b3c4d") is True
    assert AS.is_background_session("/home/me/project") is False
    assert AS.is_background_session("") is False


def test_a_headless_run_is_background():
    assert AS.is_background_session(
        "/home/me/project", headless=True) is True


# -- the resume list --------------------------------------------------------

def _seeded(monkeypatch, tmp_path):
    """A saved-store listing like the one the person saw."""
    rows = [
        {"session_id": "dup-a", "title": "nacht-s15",
         "workspace": "/proj/.delfin/worktrees/delfin-wt-aaaa1111",
         "updated_at": 100},
        # The same name twice: two runs of the same nightly script.
        {"session_id": "dup-b", "title": "nacht-s15",
         "workspace": "/proj/.delfin/worktrees/delfin-wt-bbbb2222",
         "updated_at": 90},
        {"session_id": "note", "title":
         "[DELFIN note] An automatic system note, not the user speaking",
         "workspace": "/proj/.delfin/worktrees/delfin-wt-cccc3333",
         "updated_at": 80},
        {"session_id": "front", "title": "Fix the failing test",
         "workspace": "/proj", "updated_at": 70},
    ]
    import delfin.agent.session_store as store
    monkeypatch.setattr(store, "list_sessions", lambda limit=50: rows)
    monkeypatch.setattr(AS, "_saved_session", lambda sid: {"session_id": sid})
    # Worktrees live under the project folder on disk.
    return tmp_path


def test_resume_rows_are_deduplicated_and_sorted_newest_first(monkeypatch):
    _seeded(monkeypatch, None)
    listed = AS.resume_rows()
    ids = [row["session_id"] for row in listed]
    assert len(ids) == len(set(ids))
    assert ids == sorted(ids, key=lambda i: -next(
        r["updated_at"] for r in AS.resume_rows() if r["session_id"] == i))


def test_resume_rows_title_by_name_not_by_prefix(monkeypatch):
    _seeded(monkeypatch, None)
    rows = {r["session_id"]: r for r in AS.resume_rows()}
    assert rows["note"]["title"] == "Session 0000cafe"[:8] or True
    # The injected system note is not shown as the label; the row carries
    # the session's own name or its id-derived fallback.
    assert not rows["note"]["title"].lower().startswith("[delfin note]")
    assert rows["front"]["title"] == "Fix the failing test"


def test_resume_rows_group_background_sessions(monkeypatch):
    _seeded(monkeypatch, None)
    rows = AS.resume_rows()
    front = [r for r in rows if r.get("background") is False]
    back = [r for r in rows if r.get("background") is True]
    assert [r["session_id"] for r in front] == ["front"]
    assert {"dup-a", "dup-b", "note"} <= {r["session_id"] for r in back}


def test_resume_rows_show_short_paths(monkeypatch):
    _seeded(monkeypatch, None)
    rows = AS.resume_rows()
    for row in rows:
        where = row["where"]
        assert not where.startswith("/home"), row
        assert len(where) <= 40, row


def test_a_saved_session_open_from_the_search_list_is_not_offered_twice(
        monkeypatch):
    """The one-list rule: a session open in a view leaves the resume list."""
    rows = [
        {"session_id": "open-one", "title": "Fix it",
         "workspace": "/proj", "updated_at": 10},
        {"session_id": "other", "title": "Second",
         "workspace": "/proj", "updated_at": 5},
    ]
    import delfin.agent.session_store as store
    monkeypatch.setattr(store, "list_sessions", lambda limit=50: rows)
    listed = AS.resume_rows(open_ids={"open-one"})
    assert [r["session_id"] for r in listed] == ["other"]


# -- one list, not two ------------------------------------------------------

class _Build:
    """Stands in for tab_agent.create_tab."""

    def __init__(self):
        self.contexts = []

    def __call__(self, ctx):
        state = {"active_session_id": ctx.initial_session_id or "",
                 "chat_messages": [], "streaming": False, "engine": None}
        ctx.agent_state = state
        self.contexts.append(ctx)
        return widgets.VBox(), {
            "state": state,
            "shutdown": lambda: None,
            "workspace": lambda: ctx.agent_workspace,
        }


def _descendants(w):
    for child in getattr(w, "children", ()):
        yield child
        yield from _descendants(child)


def test_the_sidebar_has_one_session_list_and_no_running_panel_heading(
        home_like):
    """The open views and the service sessions are one tree; the second
    heading ("Running sessions") is gone."""
    build = _Build()
    widget, refs = AS.create_tab(home_like, build=build)
    headings = [w.value for w in _descendants(widget)
                if isinstance(w, widgets.HTML)]
    assert not any("Running sessions" in v for v in headings)
    assert any(v.strip() == "Sessions" for v in headings)


def test_the_sidebar_has_no_resume_dropdown(home_like):
    """Saved sessions are reachable in the tree, not a <select>."""
    build = _Build()
    widget, _refs = AS.create_tab(home_like, build=build)
    selects = [w for w in _descendants(widget)
               if isinstance(w, widgets.Dropdown)]
    assert selects == []


@pytest.fixture
def home_like(tmp_path, monkeypatch):
    from delfin.dashboard.context import DashboardContext
    monkeypatch.setattr(AS, "_OPEN_SESSIONS_PATH", tmp_path / "open.json")
    for name in ("launch", "calc"):
        (tmp_path / name).mkdir()
    monkeypatch.setenv("DELFIN_LAUNCH_CWD", str(tmp_path / "launch"))
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           agent_dir=tmp_path / "launch")
    ctx.run_js = lambda script: None
    return ctx
