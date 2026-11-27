"""A caller with no session id is not told what it cannot know.

Regression test for the defect a real run measured. The branch's dead-resume-id
hint is scoped to the asking session (``subagents.py:4080``), and that scope is
right -- but the limb for "no session id at all" was written as if it had
looked and found the store empty:

    "No subagent this session started is resting either, so there is
     nothing to resume -- run the task as a fresh delegate."

An unnamed caller never looks. ``task_session_id`` is set only by the dashboard
(``tab_agent.py``) and by ``AgentEngine._sync_task_session``; the solo/CLI route
builds its permissions in ``cli_gate.py``, which has no ``task_session_id`` at
all (measured: ``grep -n task_session_id delfin/agent/cli_gate.py`` is empty).
So every tab-less caller is told there is nothing to resume, over a store that
may hold exactly the run it wants -- a false statement where the base made a
weak one.

The real run is the evidence: same prompt and model, ``resume_id='deadbeef'``
after one finished delegate. The base answered with the plain hint; the branch
answered with the "nothing is resting either" claim, on a route where no id
was ever set. ``NAMED_A_REAL_ID=no`` on both sides, and the branch's answer was
also wrong about the store.

Contracts asserted here:

* an unnamed caller gets the wording that asserts nothing about the store;
* the finished record is still there to be found -- the control that makes the
  first assertion mean something;
* a NAMED caller still gets the naming hint (the arm's own behaviour, so the
  fix did not switch the feature off).

Baseline note: the first test is RED on the merged HEAD before this fix and
GREEN after it, and GREEN on ``origin/main`` -- its subject is a wording claim
the merge introduced, so the merge is its baseline, not main. The third test is
the arm's own, unchanged.
"""

from __future__ import annotations

import json
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from delfin.agent import subagents as SA
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor


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


def _client(text="done"):
    client = MagicMock()
    client.stream_message = MagicMock(return_value=iter([
        _event(type="text_delta", text=text),
        _event(type="message_delta", input_tokens=10, output_tokens=5)]))
    client._permissions = None
    client.set_permissions = MagicMock()
    return client


def _perms(tmp_path, sid):
    """Permissions as a route carries them. ``sid=None`` is the solo/CLI case:
    the attribute is simply never set, so it reads back as ``""``."""
    p = KitToolPermissions(workspace=tmp_path, mode="default")
    if sid is not None:
        p.task_session_id = sid
    return p


def _finish_one(perms, text="the config sets retries=3"):
    """One real finished delegate, so the store holds a resumable run."""
    return SA.run_subagent(
        subagent_type="explore", description="read the config",
        prompt="read the config file and report what it sets",
        parent_client=_client(text), parent_perms=perms)


def _dead_id_error(perms):
    """The dead id through the real ``subagent`` tool dispatch."""
    perms.subagent_runner = _runner(perms)
    out = _DocToolExecutor().execute(
        "subagent",
        {"resume_id": "deadbeef", "subagent_type": "explore",
         "description": "follow-up",
         "prompt": "and which component reads that count?"},
        permissions=perms)
    return json.loads(out)["error"]


def _runner(perms):
    def _run(*, subagent_type, description, prompt, isolation="", **kw):
        return SA.run_subagent(
            subagent_type=subagent_type, description=description,
            prompt=prompt, parent_client=_client("ok"), parent_perms=perms,
            isolation=isolation, **kw).to_payload()
    return _run


def test_an_unnamed_caller_is_not_told_the_store_is_empty(stores):
    """RED before the fix: the branch claimed nothing was resting.

    Fails here on the merged HEAD before the fix with
    "No subagent this session started is resting either" in the error.
    """
    perms = _perms(stores, sid=None)
    _finish_one(_perms(stores, "sess-A"))
    err = _dead_id_error(perms)
    assert "nothing to resume" not in err, (
        "an unnamed caller was told the store is empty without having "
        f"looked: {err}")
    assert "no stored subagent session found" in err, err


def test_the_unnamed_caller_control_the_run_was_really_there(stores):
    """Control: the record exists, so the assertion above is not vacuous.

    Without this, the first test would pass on an empty store for the wrong
    reason -- it would measure that nothing was stored, not that nothing was
    claimed.
    """
    _finish_one(_perms(stores, "sess-A"))
    resting = SA.list_finished(5, session="sess-A")
    assert len(resting) == 1, resting
    assert resting[0]["description"] == "read the config"


def test_a_named_caller_still_gets_the_ids(stores):
    """The arm's own behaviour survives: naming is unchanged."""
    perms = _perms(stores, "sess-A")
    _finish_one(perms)
    err = _dead_id_error(perms)
    resting = SA.list_finished(5, session="sess-A")
    assert resting, "the fixture stored nothing to find"
    assert resting[0]["sa_id"] in err, err
