"""A session's RESTED delegates are in its archive the moment the tab opens.

Measured on the base tree with ``scratch/probe_two_sessions.py`` (two
sessions, each delegating through the real writer ``run_subagent``, then the
real tab built per session):

* sess-B's delegate never appeared in sess-A's archive -- the owner scoping
  already on main holds (``read_running(session=...)`` /
  ``list_finished(session=...)``, used at tab_agent.py:9862/9875);
* but neither did sess-A's own: right after ``create_tab`` the dropdown held
  only ``🗃 Subagent archive…``, and its own rested run appeared only after
  a LIVE delegate made the 1.5 s watcher tick once. ``/agents`` in the same
  tab listed it, so the record was there and only the archive was stale.

That is the defect these tests name: the archive was refreshed from four
places, none of which runs when the tab opens, so a session whose delegates
had all finished before this page was built -- the ordinary case after a
resume -- showed an empty archive and could not reach them.

What is asserted:

* the tab built for a session lists that session's OWN rested delegate
  (red on the base tree);
* it still does not list another session's (green on both -- the control
  that this change did not widen the scope to fix the first);
* a session that owns nothing lists nothing (green on both).

The tab is the real one: ``delfin.dashboard.tab_agent.create_tab``, the
function the dashboard calls, on a context built as the sessions sidebar
builds it. No browser and no node: what is read is the dropdown the page
renders.
"""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from delfin.agent import subagents as SA
from delfin.agent.api_client import KitToolPermissions

_ARCHIVE_MARK = "Subagent archive"


# ---------------------------------------------------------------------------
# Fixtures: an isolated store, a scripted model, and a real tab
# ---------------------------------------------------------------------------

@pytest.fixture
def stores(monkeypatch, tmp_path):
    """The sub-agent stores, redirected into the test's own directory.

    The model settings are written where the product reads them
    (``user_settings.get_settings_path()``) with a key present, and the
    scheduler is this test's own: without them the tab builds a CLI engine
    whose binary the tripwire refuses, the restore aborts before the
    archive is drawn, and the test would fail for a reason that has nothing
    to do with what it names.
    """
    from delfin import user_settings
    from delfin.agent import scheduler as S

    monkeypatch.setattr(SA, "_RUNNING_DIR", tmp_path / "running")
    monkeypatch.setattr(SA, "_SESSIONS_DIR", tmp_path / "sessions")
    monkeypatch.setattr(SA, "_PENDING_DIR", tmp_path / "pending")
    monkeypatch.setattr(SA, "_TELEMETRY_PATH", tmp_path / "telemetry.jsonl")
    SA._FINISHED_SUMMARY_CACHE.clear()
    SA._NOTES.clear()
    (tmp_path / "running").mkdir(parents=True, exist_ok=True)
    monkeypatch.setattr(S, "_GLOBAL", S.Scheduler(path=tmp_path / "cron.json"))
    monkeypatch.setenv("KIT_TOOLBOX_API_KEY", "sk-test")
    path = user_settings.get_settings_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"agent": {
        "provider": "kit", "model": "kit.deepseek-v4-flash",
        "permission_profile": "all_free", "probe_endpoint": False}}))
    yield tmp_path
    SA._NOTES.clear()


def _event(**kw):
    return SimpleNamespace(
        type=kw.get("type", "text_delta"), text=kw.get("text", ""),
        tool_name=kw.get("tool_name", ""), tool_input=kw.get("tool_input", ""),
        tool_output=kw.get("tool_output", ""),
        input_tokens=0, output_tokens=0)


def _client(text: str):
    client = MagicMock()

    def _stream(**_kw):
        yield _event(type="tool_use", tool_name="read_file", tool_input="{}")
        yield _event(type="text_delta", text=text)

    client.stream_message = MagicMock(side_effect=_stream)
    client._permissions = None
    client.set_permissions = MagicMock()
    return client


def _delegate(tmp_path, session_id: str, description: str) -> str:
    """One finished delegate, written by the real writer."""
    perms = KitToolPermissions(workspace=tmp_path, mode="plan")
    perms.task_session_id = session_id
    result = SA.run_subagent(
        subagent_type="explore", description=description,
        prompt="look at " + description,
        parent_client=_client("done: " + description), parent_perms=perms)
    assert result.error == ""
    return result.sa_id


def _save_dashboard_session(tmp_path, session_id: str) -> None:
    """A saved dashboard session, so the tab opens ON it as a resume does."""
    from delfin.agent import session_store as ss
    ss.save_session(
        session_id, mode="solo",
        chat_messages=[{"role": "user", "content": "delegate something",
                        "role_label": ""},
                       {"role": "assistant", "content": "done",
                        "role_label": ""}],
        engine_messages=[{"role": "user", "content": "delegate something"},
                         {"role": "assistant", "content": "done"}],
        provider="kit", model="kit.deepseek-v4-flash", perm_profile="all_free",
        workspace=str(tmp_path / "ws"))
    ss.release_session_lock(session_id)


def _walk(widget, seen=None):
    seen = [] if seen is None else seen
    seen.append(widget)
    for child in (getattr(widget, "children", ()) or ()):
        _walk(child, seen)
    return seen


def _archive_rows(tab) -> list[tuple[str, str]]:
    """The ``(label, sa_id)`` rows of the tab's Subagent-archive dropdown."""
    import ipywidgets as widgets
    for node in _walk(tab):
        if not isinstance(node, widgets.Dropdown):
            continue
        opts = list(getattr(node, "options", ()) or ())
        labels = [str(o[0]) if isinstance(o, (tuple, list)) else str(o)
                  for o in opts]
        if any(_ARCHIVE_MARK in label for label in labels):
            return [(label, (o[1] if isinstance(o, (tuple, list)) else o))
                    for label, o in zip(labels, opts)]
    return []


def _open_tab(tmp_path, session_id: str = ""):
    """Build the real agent tab, as the sessions sidebar builds it."""
    from delfin.agent import scheduler as S
    from delfin.dashboard import tab_agent
    from delfin.dashboard.context import DashboardContext

    for name in ("calc", "archive", "office", "ws"):
        (tmp_path / name).mkdir(exist_ok=True)
    S._GLOBAL = S.Scheduler(path=tmp_path / "cron.json")
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda script, **kw: None
    ctx.agent_workspace = str(tmp_path / "ws")
    ctx.initial_session_id = session_id
    ctx.sessions_managed = True
    ctx.presence_key = (session_id[:8] or "test")
    tab, refs = tab_agent.create_tab(ctx)
    return tab, refs


def _ids(rows) -> list[str]:
    return [str(value) for _label, value in rows if value]


# ---------------------------------------------------------------------------
# The archive at the moment the tab opens
# ---------------------------------------------------------------------------

def test_the_tab_lists_its_own_rested_delegate_when_it_opens(stores):
    """The run finished before this page existed -- the resume case."""
    mine = _delegate(stores, "sess-A", "read the pack")
    _save_dashboard_session(stores, "sess-A")
    tab, refs = _open_tab(stores, "sess-A")
    try:
        rows = _archive_rows(tab)
        assert mine in _ids(rows), (
            "the session's own finished delegate is not in its archive as the"
            f" tab opens; the archive holds {rows!r}")
        assert any("read the pack" in label for label, _v in rows), rows
    finally:
        refs["shutdown"](save=False)


def test_the_open_tab_does_not_list_another_sessions_rested_delegate(stores):
    """The control: the fix must not be "list everything"."""
    theirs = _delegate(stores, "sess-B", "work of B")
    mine = _delegate(stores, "sess-A", "read the pack")
    _save_dashboard_session(stores, "sess-A")
    tab, refs = _open_tab(stores, "sess-A")
    try:
        ids = _ids(_archive_rows(tab))
        assert mine in ids, f"its own delegate is missing: {ids!r}"
        assert theirs not in ids, (
            f"another session's delegate is listed in this session's archive:"
            f" {ids!r}")
    finally:
        refs["shutdown"](save=False)


def test_a_session_that_owns_no_delegate_lists_none(stores):
    """A new session owns nothing and must not fall back to the account."""
    _delegate(stores, "sess-B", "work of B")
    _save_dashboard_session(stores, "sess-A")
    tab, refs = _open_tab(stores, "sess-A")
    try:
        assert _ids(_archive_rows(tab)) == [], _archive_rows(tab)
    finally:
        refs["shutdown"](save=False)


def test_the_archive_refresh_is_reached_without_a_live_delegate(monkeypatch,
                                                                stores):
    """Why the two above failed: nothing refreshed the archive on open.

    Named as its own case so a future change that fixes the archive by
    another route (a widget observer, a page hook) does not have to keep
    this call -- the three tests above are the contract; this one records
    the mechanism that was missing.
    """
    from delfin.dashboard import tab_agent
    calls = {"n": 0}
    real = tab_agent.archive_rows

    def _counting(*a, **kw):
        calls["n"] += 1
        return real(*a, **kw)

    monkeypatch.setattr(tab_agent, "archive_rows", _counting)
    _delegate(stores, "sess-A", "read the pack")
    _save_dashboard_session(stores, "sess-A")
    tab, refs = _open_tab(stores, "sess-A")
    try:
        assert calls["n"] >= 1, (
            "the archive was never built while the tab opened, so a session"
            " with only finished delegates shows an empty archive")
    finally:
        refs["shutdown"](save=False)


def test_the_store_scopes_the_query_by_the_owning_session(stores):
    """The store the archive reads is keyed by the owning session.

    Read at the level the archive uses it, so the value passed to
    ``list_finished`` is shown to be the session id the delegate recorded --
    not a timestamp window or a process stamp.
    """
    mine = _delegate(stores, "sess-A", "read the pack")
    theirs = _delegate(stores, "sess-B", "work of B")
    scoped = SA.list_finished(session="sess-A")
    assert [r["sa_id"] for r in scoped] == [mine]
    assert theirs not in [r["sa_id"] for r in scoped]
    # ...and the record says whose it is, not merely "not mine".
    stored = json.loads(
        (Path(SA._SESSIONS_DIR) / f"{theirs}.json").read_text())
    assert stored["owner_session"] == "sess-B"
