"""Sub-agent listings are scoped to the session whose turns started them.

Defect: the dashboard's live chip row (running delegates) and the
"Subagent archive" dropdown read ``read_running()`` and ``list_finished()``
unfiltered, so every session of the account listed every delegate of the
account. Only background runs carried ``owner_session``; foreground runs,
orchestration and fan-out workers, and the finished store carried none.

Contract asserted here:

* every running entry and every finished record carries ``owner_session``,
  the TOP-LEVEL session id (``perms.task_session_id``); a nested delegate
  inherits it through the derived permissions;
* ``read_running(session=...)`` / ``list_finished(session=...)`` return that
  session's records only; ``list_finished(unassigned=True)`` returns the
  records written before ownership existed; no keyword keeps the old
  account-wide answer for callers that look runs up by id;
* the dashboard archive lists the session's own rows, then legacy rows under
  a separate, non-selectable fold header -- never another session's rows;
* foreground runs stay out of the Background panel (as before: they had no
  owner and the panel is strict), so the panel's height does not change.
"""

from __future__ import annotations

import inspect
import json
import os
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from delfin.agent import subagents as SA
from delfin.agent.api_client import KitToolPermissions


@pytest.fixture
def stores(monkeypatch, tmp_path):
    monkeypatch.setattr(SA, "_RUNNING_DIR", tmp_path / "running")
    monkeypatch.setattr(SA, "_SESSIONS_DIR", tmp_path / "sessions")
    monkeypatch.setattr(SA, "_PENDING_DIR", tmp_path / "pending")
    monkeypatch.setattr(SA, "_TELEMETRY_PATH", tmp_path / "telemetry.jsonl")
    SA._FINISHED_SUMMARY_CACHE.clear()
    SA._NOTES.clear()
    yield tmp_path
    SA._NOTES.clear()


def _event(**kw):
    return SimpleNamespace(
        type=kw.get("type", "text_delta"), text=kw.get("text", ""),
        tool_name=kw.get("tool_name", ""), tool_input=kw.get("tool_input", ""),
        tool_output=kw.get("tool_output", ""),
        input_tokens=0, output_tokens=0)


def _client(text="done", during=None):
    """A scripted client; ``during`` runs while the delegate is live."""
    client = MagicMock()

    def _stream(**_kw):
        yield _event(type="tool_use", tool_name="read_file", tool_input="{}")
        if during is not None:
            during()
        yield _event(type="text_delta", text=text)

    client.stream_message = MagicMock(side_effect=_stream)
    client._permissions = None
    client.set_permissions = MagicMock()
    return client


def _perms(tmp_path, sid):
    p = KitToolPermissions(workspace=tmp_path, mode="plan")
    p.task_session_id = sid
    return p


def _run(tmp_path, sid, desc, during=None, perms=None):
    return SA.run_subagent(
        subagent_type="explore", description=desc, prompt="look " + desc,
        parent_client=_client(during=during),
        parent_perms=perms if perms is not None else _perms(tmp_path, sid))


# ---------------------------------------------------------------------------
# The stores record the owner
# ---------------------------------------------------------------------------

def test_two_sessions_each_see_only_their_own_finished_delegate(stores):
    a = _run(stores, "sess-A", "task of A")
    b = _run(stores, "sess-B", "task of B")
    assert a.error == "" and b.error == ""
    own_a = [r["description"] for r in SA.list_finished(session="sess-A")]
    own_b = [r["description"] for r in SA.list_finished(session="sess-B")]
    assert own_a == ["task of A"]
    assert own_b == ["task of B"]
    # Unfiltered keeps the account-wide answer for the by-id callers.
    assert {r["description"] for r in SA.list_finished()} == {
        "task of A", "task of B"}


def test_a_running_foreground_delegate_is_listed_for_its_session_only(stores):
    seen: dict = {}

    def _peek():
        seen["A"] = SA.read_running(session="sess-A")
        seen["B"] = SA.read_running(session="sess-B")
        seen["all"] = SA.read_running()

    _run(stores, "sess-A", "live in A", during=_peek)
    assert [e["description"] for e in seen["A"].values()] == ["live in A"]
    assert seen["B"] == {}
    assert len(seen["all"]) == 1
    entry = next(iter(seen["A"].values()))
    assert entry["owner_session"] == "sess-A"
    # A foreground run is marked so the Background panel keeps leaving it
    # out, exactly as it did while such runs carried no owner at all.
    assert entry.get("foreground") is True


def test_a_nested_delegate_inherits_the_top_level_owner(stores):
    top = _perms(stores, "sess-top")
    child = SA._derive_perms(top, "plan", preset_tools=(), preset_name="explore")
    assert child.task_session_id == "sess-top"
    _run(stores, "", "nested", perms=child)
    assert [r["description"] for r in SA.list_finished(session="sess-top")] == [
        "nested"]


def test_a_background_reservation_owner_wins_over_a_later_perms_change(stores):
    """The background thread reads the client's CURRENT permissions; a tab
    switched to another session in between must not re-own the run."""
    SA.reserve_running("bg1", subagent_type="explore", description="bg",
                       owner_session="sess-A")
    SA.run_subagent(subagent_type="explore", description="bg", prompt="p",
                    parent_client=_client(), parent_perms=_perms(stores, "sess-B"),
                    sa_id="bg1")
    rec = json.loads((stores / "sessions" / "bg1.json").read_text())
    assert rec["owner_session"] == "sess-A"
    assert SA.list_finished(session="sess-B") == []


def test_a_run_without_a_session_owns_nothing(stores):
    _run(stores, "", "no session", perms=None)
    rec = SA.list_finished()
    assert rec and rec[0]["owner_session"] == ""
    assert SA.list_finished(session="sess-A") == []
    # ...and it is filed with the unassigned ones, not shown everywhere.
    assert [r["description"] for r in SA.list_finished(unassigned=True)] == [
        "no session"]


def test_legacy_records_are_unassigned_not_everyones(stores):
    d = stores / "sessions"
    d.mkdir(parents=True)
    (d / "old1.json").write_text(json.dumps({
        "sa_id": "old1", "subagent_type": "explore", "description": "legacy",
        "finished_at": 1.0}))
    _run(stores, "sess-A", "mine")
    assert [r["sa_id"] for r in SA.list_finished(unassigned=True)] == ["old1"]
    assert "old1" not in {r["sa_id"] for r in SA.list_finished(session="sess-A")}
    (stores / "running").mkdir(parents=True, exist_ok=True)
    (stores / "running" / "old2.json").write_text(json.dumps({
        "type": "explore", "description": "legacy live",
        "owner_pid": os.getpid()}))
    assert "old2" in SA.read_running()
    assert SA.read_running(session="sess-A") == {}


def test_the_limit_counts_the_sessions_own_records(stores):
    """Filtering after the cut would let another session's newer runs push
    this session's out of its own archive."""
    _run(stores, "sess-A", "early A")
    for i in range(5):
        _run(stores, "sess-B", f"B{i}")
    assert [r["description"] for r in SA.list_finished(3, session="sess-A")] == [
        "early A"]


# ---------------------------------------------------------------------------
# The dashboard
# ---------------------------------------------------------------------------

def test_the_archive_lists_own_rows_then_a_legacy_fold():
    from delfin.dashboard import tab_agent as T
    own = [{"sa_id": "a1", "subagent_type": "explore", "description": "mine"}]
    legacy = [{"sa_id": "o1", "subagent_type": "explore", "description": "old"}]
    rows = T.archive_rows({}, own, legacy)
    assert [v for _, v in rows] == ["a1", T._ARCHIVE_UNASSIGNED_HEADER, "o1"]
    assert T.archive_rows({}, own) == T.archive_rows({}, own, [])
    assert T._ARCHIVE_UNASSIGNED_HEADER not in [
        v for _, v in T.archive_rows({}, own, [])]


def test_the_dashboard_reads_scoped_lists():
    from delfin.dashboard import tab_agent as T
    src = inspect.getsource(T)
    body = src[src.index("def _refresh_agent_view_chips():"):
               src.index("agent_archive_dropdown.observe")]
    assert "read_running(session=_owner)" in body
    assert "list_finished(_ARCHIVE_LIMIT, session=_owner)" in body
    assert "unassigned=True" in body
    # The fold header is a label, never a view a pick can open.
    pick = src.split("def _on_archive_pick")[1].split("\n\n")[0]
    assert "_ARCHIVE_UNASSIGNED_HEADER" in pick
    agents_cmd = src[src.index('if cmd == "/agents":'):][:3000]
    assert "list_finished(last_n=5, session=_background_owner())" in agents_cmd


def test_the_background_panel_keeps_leaving_foreground_runs_out(monkeypatch, tmp_path):
    from delfin.agent import background_view as BV
    from delfin.agent import bash_jobs as BJ
    from delfin.agent import job_monitor as JM
    from delfin.agent import scheduler as SCH
    monkeypatch.setattr(JM, "_AGENT_WATCH_INDEX_PATH", tmp_path / "index.json")
    monkeypatch.setattr(BJ, "live_jobs", lambda: [])
    monkeypatch.setattr(SCH, "get_scheduler",
                        lambda: SimpleNamespace(list_entries=lambda: []))
    monkeypatch.setattr(SA, "read_running", lambda **_: {
        "fg": {"type": "explore", "description": "fg", "owner_session": "s1",
               "foreground": True},
        "bg": {"type": "explore", "description": "bg", "owner_session": "s1"}})
    monkeypatch.setattr(SA, "_entry_owned_by_us", lambda entry: True)
    agents = BV.collect(tmp_path, session_id="s1")["agents"]
    assert [a["id"] for a in agents] == ["bg"]


# ---------------------------------------------------------------------------
# Other listings
# ---------------------------------------------------------------------------

def test_the_terminal_agents_command_lists_its_own_session(stores, monkeypatch):
    from delfin.agent import repl_commands as RC
    _run(stores, "sess-A", "A's run")
    _run(stores, "sess-B", "B's run")
    ctx = SimpleNamespace(engine=SimpleNamespace(session_id="sess-A",
                                                 kit_permissions=None))
    out = RC._agents(ctx, "").output
    assert "A's run" in out and "B's run" not in out


def test_subagent_message_reaches_only_this_sessions_delegates(stores, monkeypatch):
    from delfin.agent import api_client as AC
    SA.reserve_running("dA", subagent_type="explore", description="a",
                       owner_session="sess-A")
    SA.reserve_running("dB", subagent_type="explore", description="b",
                       owner_session="sess-B")
    perms = _perms(stores, "sess-A")
    listed = json.loads(AC._doc_executor._execute_subagent_message({}, perms))
    assert [d["sa_id"] for d in listed["delegates"]] == ["dA"]
    sent = json.loads(AC._doc_executor._execute_subagent_message(
        {"to": "all", "message": "hi"}, perms))
    assert sent["to"] == ["dA"]


def test_subagent_message_cannot_direct_another_sessions_delegate(stores):
    from delfin.agent import api_client as AC
    SA.reserve_running("dB", subagent_type="explore", description="b",
                       owner_session="sess-B")
    out = json.loads(AC._doc_executor._execute_subagent_message(
        {"to": "dB", "message": "push it"}, _perms(stores, "sess-A")))
    assert "error" in out and SA.notes_waiting("dB") == 0
