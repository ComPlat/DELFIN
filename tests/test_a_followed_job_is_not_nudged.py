"""A turn that hands its job to DELFIN's watcher is not told nothing runs.

Seen while recording the cluster video (2026-11-12): cluster_submit queued
the batch and told the agent to end its turn; the agent closed on a wait
("Let me wait for the job to finish"), and the window showed "System note: this turn
ended by announcing an action" and "Nothing is running now ... Send
'weiter'". The batch was being followed; the nudge only knew about
wait_registry, and the watch is job_monitor's.

Driven through the real AgentEngine.stream_response.
"""

from __future__ import annotations

from types import SimpleNamespace
from unittest.mock import patch

from delfin.agent import engine as E
from delfin.agent import job_monitor
from delfin.agent.api_client import StreamEvent

CLOSING = "The batch is queued for hpc-a. Let me wait for the job to finish."


def _client(text):
    class _C:
        def stream_message(self, *_a, **_k):
            yield StreamEvent(type="text_delta", text=text)
            yield StreamEvent(type="message_delta", stop_reason="end_turn",
                              input_tokens=1, output_tokens=1)
    return _C()


def _engine(tmp_path, session="sess-1"):
    with patch("delfin.agent.engine.create_client", return_value=_client(CLOSING)):
        eng = E.AgentEngine(repo_dir=str(tmp_path), backend="api", provider="kit",
                            model="kit.deepseek-v4-flash", mode="solo")
    eng._open_task_count = lambda: 0
    # The permissions a KIT client carries: the workspace and this session.
    eng.client._permissions = SimpleNamespace(workspace=str(tmp_path), task_session_id=session)
    return eng


def _watch(eng, session, state=""):
    job_monitor.register_agent_job(
        eng.kit_permissions.workspace, "hpc-a:batch-1", description="1 job(s) for hpc-a",
        extra={"kind": "cluster", "cluster": "hpc-a", "remote_id": "batch-1",
               "session_id": session, "last_state": state})


def _turn(eng):
    shown = []
    eng.stream_response("Submit benzene to hpc-a", on_token=shown.append)
    return "".join(shown)


def test_a_followed_batch_is_neither_nudged_nor_called_nothing(tmp_path):
    eng = _engine(tmp_path)
    _watch(eng, "sess-1")
    shown = _turn(eng)
    assert not eng.pending_turn_continuation
    assert "Nothing is running now" not in shown


def test_without_a_watch_the_same_close_is_still_nudged(tmp_path):
    eng = _engine(tmp_path)
    _turn(eng)
    assert eng.pending_turn_continuation


def test_another_sessions_watch_or_an_ended_job_does_not_count(tmp_path):
    eng = _engine(tmp_path)
    _watch(eng, "sess-2")
    _turn(eng)
    assert eng.pending_turn_continuation
    eng.clear_turn_continuation()
    _watch(eng, "sess-1", state="COMPLETED")
    _turn(eng)
    assert eng.pending_turn_continuation
