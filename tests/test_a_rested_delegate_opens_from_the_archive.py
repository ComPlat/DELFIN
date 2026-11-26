"""A RESTED delegate listed in a session's archive can be opened from there.

The package asks for two things, and the two other files cover only the first:
the archive must list this session's own finished subagents and nothing else.
This one covers the second -- that a listed one is *reachable*, not merely
visible. ``archive_rows`` producing a ``(label, sa_id)`` pair proves a row was
built; it says nothing about whether selecting that row opens the transcript.

What is asserted here:

* selecting the row a session's own rested delegate was listed under enters
  that delegate's transcript (red on the base tree: there is no such row to
  select, because the archive was never drawn while the tab opened);
* selecting the placeholder changes nothing -- the control that passes on both
  trees, so a green result above is not just "any selection succeeds".

The tab is the real ``delfin.dashboard.tab_agent.create_tab``, the function the
dashboard calls, driven through the dropdown's own ``value`` trait the way the
browser drives it. No browser and no node.
"""

from __future__ import annotations

import json
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from delfin.agent import subagents as SA
from delfin.agent.api_client import KitToolPermissions

_ARCHIVE_MARK = "Subagent archive"


@pytest.fixture
def stores(monkeypatch, tmp_path):
    """The sub-agent stores and model settings, redirected into this test."""
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


def _archive_dropdown(tab):
    import ipywidgets as widgets
    for node in _walk(tab):
        if not isinstance(node, widgets.Dropdown):
            continue
        opts = list(getattr(node, "options", ()) or ())
        labels = [str(o[0]) if isinstance(o, (tuple, list)) else str(o)
                  for o in opts]
        if any(_ARCHIVE_MARK in label for label in labels):
            return node
    return None


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


# ---------------------------------------------------------------------------
# Opening a rested delegate from its archive row
# ---------------------------------------------------------------------------

def test_a_rested_delegate_opens_from_its_archive_row(stores):
    """Selecting the row enters that delegate's transcript, not Main."""
    mine = _delegate(stores, "sess-A", "read the pack")
    _save_dashboard_session(stores, "sess-A")
    tab, refs = _open_tab(stores, "sess-A")
    try:
        dropdown = _archive_dropdown(tab)
        assert dropdown is not None, "the tab has no Subagent-archive dropdown"
        values = [o[1] if isinstance(o, (tuple, list)) else o
                  for o in (dropdown.options or ())]
        assert mine in values, (
            "the session's own rested delegate is not selectable in its "
            f"archive as the tab opens; the dropdown holds {values!r}")

        # What the browser does when the row is picked.
        dropdown.value = mine
        assert refs["state"].get("_view_agent") == mine, (
            "picking the archive row did not enter the delegate's transcript; "
            f"the view is {refs['state'].get('_view_agent')!r}")
    finally:
        refs["shutdown"](save=False)


def test_the_placeholder_selection_changes_nothing(stores):
    """The control: the empty row is not an entry, on either tree."""
    _delegate(stores, "sess-A", "read the pack")
    _save_dashboard_session(stores, "sess-A")
    tab, refs = _open_tab(stores, "sess-A")
    try:
        dropdown = _archive_dropdown(tab)
        assert dropdown is not None, "the tab has no Subagent-archive dropdown"
        before = refs["state"].get("_view_agent")
        dropdown.value = ""
        assert refs["state"].get("_view_agent") == before, (
            "selecting the placeholder moved the view to "
            f"{refs['state'].get('_view_agent')!r}")
    finally:
        refs["shutdown"](save=False)
