"""A new session's archive does not list the pre-attribution delegates.

Measured on the base tree (``80a43c7ff``), the real tab built once per
session: every session's archive carried a ``── older, unassigned ──``
fold holding EVERY pre-attribution delegate of the account. A session
that had never delegated anything therefore opened with one row it never
started, and a second session in another window saw the same rows. The
package's rule is that the tab "must not show what is not true for the
session in view", and "this conversation's delegates" is not true of
records that belong to no conversation.

The fold is account-wide BY DESIGN at the store level
(``list_finished(unassigned=True)``) and stays so: it is the one way to
list records written before ownership was recorded. What is wrong is the
dashboard passing it for every session, so the tab shows it.

What is asserted:

* a session that owns no delegate, opened on a store that holds one
  pre-attribution record, lists no such row (RED on the base tree);
* the store's own account-wide query still returns it, so the record is
  not lost and a caller that lists it deliberately still can (the
  control: passes on both trees);
* an owned delegate still appears while the fold does not (both).

No browser and no node: what is read is the dropdown the page renders.
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


@pytest.fixture
def stores(monkeypatch, tmp_path):
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


def _a_record_without_an_owner(tmp_path, description: str) -> str:
    """A pre-attribution finished record: written before ``owner_session``.

    Written through the product's own saver with no owner, which is the
    shape the 2026-11-03 attribution change left on disk, rather than by
    hand-editing a file.
    """
    import time
    sa_id = "old%d" % int(time.time() * 1000 % 100000)
    SA._save_subagent_session(
        sa_id, subagent_type="explore", description=description,
        messages=[{"role": "user", "content": "look at " + description},
                  {"role": "assistant", "content": "done"}],
        interactions=[], error="", owner_session="")
    return sa_id


def _save_dashboard_session(tmp_path, session_id: str) -> None:
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


def test_a_new_session_does_not_list_the_pre_attribution_fold(stores):
    """The leak: a session that never delegated shows old records anyway."""
    _a_record_without_an_owner(stores, "worked on something else entirely")
    _save_dashboard_session(stores, "sess-new")
    tab, refs = _open_tab(stores, "sess-new")
    try:
        rows = _archive_rows(tab)
        labels = [label for label, _v in rows]
        assert not any("unassigned" in label for label in labels), (
            "a new session's archive lists the account's pre-attribution"
            f" delegates: {rows!r}")
        # Nothing selectable but the placeholder: every other row is an
        # sa_id this session could be sent into.
        selectable = [v for _l, v in rows if v]
        assert selectable == [], (
            "a session that delegated nothing has rows to open:"
            f" {rows!r}")
    finally:
        refs["shutdown"](save=False)


def test_the_store_still_lists_the_pre_attribution_records(stores):
    """The control: the fold is not deleted, only not shown per session."""
    old = _a_record_without_an_owner(stores, "worked on something else")
    assert [r["sa_id"] for r in SA.list_finished(unassigned=True)] == [old]
    # ...and it stays out of a scoped listing, as it always was.
    assert SA.list_finished(session="sess-new") == []


def test_its_own_delegate_is_still_listed_while_the_fold_is_not(stores):
    """Taking the fold out must not take the session's own rows with it."""
    mine = _delegate(stores, "sess-A", "read the pack")
    _a_record_without_an_owner(stores, "worked on something else entirely")
    _save_dashboard_session(stores, "sess-A")
    tab, refs = _open_tab(stores, "sess-A")
    try:
        rows = _archive_rows(tab)
        values = [str(v) for _l, v in rows]
        labels = [label for label, _v in rows]
        assert mine in values, f"its own delegate is missing: {rows!r}"
        assert not any("unassigned" in label for label in labels), rows
    finally:
        refs["shutdown"](save=False)
