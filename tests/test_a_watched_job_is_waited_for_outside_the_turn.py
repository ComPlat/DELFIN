"""A watched job is waited for outside the turn.

Report 20260915-125310: the agent started the test suite in the background,
watched it, and said it would end its turn -- and auto-continue sent it
straight back in, because tasks were still open. The only thing left to do
was wait, so it waited inside the turn: six bash_status calls blocking 300 s
each, on a model where every round costs minutes. While a watched job is
still out, the turn ends; the job wakes the agent when it is done.
"""

from __future__ import annotations

import ast
import inspect
import pathlib
import time
from types import SimpleNamespace

import pytest

from delfin.agent import api_client as A
from delfin.agent import job_monitor as jm


@pytest.fixture
def ws(tmp_path, monkeypatch):
    monkeypatch.setattr(jm, "_AGENT_WATCH_INDEX_PATH", tmp_path / "index.json")
    return tmp_path


def _waiting(ws) -> bool:
    client = SimpleNamespace(_permissions=SimpleNamespace(workspace=ws))
    return A.OpenAIClient._waiting_on_watched_jobs(client)


def test_nothing_watched_is_nothing_to_wait_for(ws):
    assert _waiting(ws) is False


def test_a_watched_cluster_job_is_waited_for(ws):
    jm.register_agent_job(ws, "4976064", "opt freq")
    assert _waiting(ws) is True


def test_a_watched_ci_run_is_waited_for(ws):
    jm.register_ci_watch(ws, "ComPlat/DELFIN", "d4ecf365", "main")
    assert _waiting(ws) is True


def test_a_background_job_counts_only_while_it_runs(ws, monkeypatch):
    from delfin.agent import bash_jobs as BJ
    monkeypatch.setattr(BJ, "_INDEX_PATH", ws / "bash_jobs_index.json")
    job = BJ.get_registry().start("sleep 30", cwd=str(ws), workspace=ws)
    try:
        jm.register_agent_job(ws, job.job_id, "suite")
        assert _waiting(ws) is True
    finally:
        BJ.get_registry().kill(job.job_id)
    deadline = time.monotonic() + 10
    while job.poll() is None and time.monotonic() < deadline:
        time.sleep(0.05)
    assert _waiting(ws) is False


def test_auto_continue_asks_before_it_sends_the_agent_back_in():
    src = inspect.getsource(A.OpenAIClient.stream_message)
    i = src.index("_did_tools_since_cont and _auto_cont_count < _AUTO_CONT_CAP")
    assert "not self._waiting_on_watched_jobs()" in src[i:i + 300]
