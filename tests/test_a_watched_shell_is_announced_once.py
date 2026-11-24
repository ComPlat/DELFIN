"""One finished shell is told to the agent once, with one outcome.

The product's own path, end to end:

1. the model's calls go through the tool dispatch
   (``KitToolClient.execute`` -- what a tool call reaches),
2. the job reaches a terminal state,
3. the idle prompt's wake (``TerminalAgent._wake_text``, the method the
   terminal calls between turns) is read, and its text is the notice the
   agent's next turn opens on.

Measured 2026-11-23 through that path, on a job started with
``bash_background`` and then watched with ``watch_job`` -- the sequence the
``watch_job`` answer itself prescribes ("You will be notified in a future
turn"). The notice named the same job twice:

    [watch -- a system event, not the user] A job you were watching has
    finished:
    - shell 8048adee [exit 126] echo hello; exit 0
    - bash 8048adee [FINISHED]

The ``shell`` line is the job's real outcome. The ``bash`` line says
FINISHED for a job whose exit code is 126: it reports that the job reached
a terminal state, not what the state was, and it reads as success beside a
failure. One finished job, two rows, one of them wrong.

``background_view.collect`` already decided this question the other way:
a bash job gets ONE row, and the watch entry for it is skipped
(background_view.py: the ``kind == "bash"`` branch). The notice and the
panel are two renderings of the same state and only one of them counted.

The de-duplication lives here, not in each caller, because the callers are
four (the terminal prompt, the headless ``--job-wait`` loop, the
dashboard's wake tick, and the panel) and the headless one had already
grown its own copy of this fix. A rule applied at one caller is a rule the
next caller breaks.
"""

from __future__ import annotations

import json
import time
from types import SimpleNamespace

import pytest

from delfin.agent import repl as R
from delfin.agent.api_client import KitToolPermissions, _doc_executor


SESSION = "f2-session"


@pytest.fixture
def job_registry(tmp_path, monkeypatch):
    """The job index in tmp, never the developer's own (bash_jobs fixes its
    index path at import time, so this is the registry object, not HOME)."""
    from delfin.agent import bash_jobs as bj
    monkeypatch.setattr(bj, "_INDEX_PATH", tmp_path / "bash_jobs_index.json")
    return bj


@pytest.fixture
def workspace(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    return ws


def _run_the_calls(ws):
    """The model's two calls, through the dispatch, and the finished job."""
    perms = KitToolPermissions(workspace=ws, mode="default",
                               task_session_id=SESSION)
    started = json.loads(_doc_executor.execute(
        "bash_background",
        {"command": "echo nothing-is-run-in-this-test; exit 3",
         "description": "the work that was waiting"},
        perms))
    assert started.get("status") == "started", started
    job_id = started["job_id"]

    for _ in range(200):
        status = json.loads(_doc_executor.execute(
            "bash_status", {"job_id": job_id}, perms))
        if not status.get("running"):
            break
        time.sleep(0.02)
    assert status.get("running") is False, status

    watched = json.loads(_doc_executor.execute(
        "watch_job", {"job_id": job_id}, perms))
    assert watched.get("status") == "watching", watched
    return job_id, status["exit_code"]


def _notice(ws):
    """What the terminal's idle prompt hands the agent as its next turn."""
    agent = R.TerminalAgent.__new__(R.TerminalAgent)
    agent.engine = SimpleNamespace(session_id=SESSION)
    agent.opts = SimpleNamespace(cwd=str(ws))
    agent._wake_last_look = 0.0
    agent._watch_last_look = -1e9
    return str(agent._wake_text("") or "")


def _lines_for(notice, job_id):
    """The rows of the notice that name ``job_id``."""
    return [ln for ln in notice.splitlines()
            if ln.lstrip().startswith("- ") and job_id in ln]


def test_one_finished_watched_shell_is_named_once(job_registry, workspace):
    job_id, exit_code = _run_the_calls(workspace)

    notice = _notice(workspace)

    assert notice, "the finished job produced no wake-up at all"
    rows = _lines_for(notice, job_id)
    assert len(rows) == 1, (
        "one finished job was named " + str(len(rows)) + " times: "
        + repr(rows))


def test_the_one_row_carries_the_real_outcome(job_registry, workspace):
    """The control for the row above: the surviving row is the one that
    says how the job actually ended, not "FINISHED" for a failed job."""
    job_id, exit_code = _run_the_calls(workspace)

    rows = _lines_for(_notice(workspace), job_id)

    assert len(rows) == 1, rows
    assert f"[exit {exit_code}]" in rows[0], (
        "the row kept for a job that exited "
        + str(exit_code) + " does not say so: " + repr(rows[0]))


def test_a_shell_nobody_watched_is_still_told(job_registry, workspace):
    """The control that must pass before and after: the fix removes a
    duplicate, never the notice. A finished shell with no ``watch_job``
    still wakes the session that started it."""
    perms = KitToolPermissions(workspace=workspace, mode="default",
                               task_session_id=SESSION)
    started = json.loads(_doc_executor.execute(
        "bash_background",
        {"command": "echo unwatched; exit 0", "description": "unwatched work"},
        perms))
    job_id = started["job_id"]
    for _ in range(200):
        status = json.loads(_doc_executor.execute(
            "bash_status", {"job_id": job_id}, perms))
        if not status.get("running"):
            break
        time.sleep(0.02)

    rows = _lines_for(_notice(workspace), job_id)

    assert len(rows) == 1, rows
    # The shell row carries the command (job_wake.finished_shells), so that
    # is what identifies it here -- whether the command got to run or was
    # refused by the sandbox, the notice is about the same job either way.
    assert "echo unwatched" in rows[0], rows[0]


def _model_turn_that_sees(notice, ws):
    """Run one turn with the real streaming loop over a stub backend, and
    return the messages the model was shown a second time.

    The tool loop is ``OpenAIClient.stream_message``; only the backend is a
    double, the same shape ``test_memory_context_nudge_integration`` uses.
    This is the half the wake path is FOR: the notice is worth nothing if
    the model's context never receives it.
    """
    import threading
    import types

    from delfin.agent import api_client as ac

    state = {"round": 0, "requests": []}

    class _Stream:
        def __init__(self, events):
            self._events = events

        def __iter__(self):
            for ev in self._events:
                yield ev

        def close(self):
            pass

    def _answer(text):
        return [
            types.SimpleNamespace(
                usage=None,
                choices=[types.SimpleNamespace(
                    delta=types.SimpleNamespace(
                        content=text, tool_calls=None, reasoning_content=None),
                    finish_reason=None)]),
            types.SimpleNamespace(
                usage=types.SimpleNamespace(prompt_tokens=5, completion_tokens=3),
                choices=[types.SimpleNamespace(
                    delta=types.SimpleNamespace(
                        content=None, tool_calls=None, reasoning_content=None),
                    finish_reason="stop")]),
        ]

    class _Stub:
        chat = types.SimpleNamespace(
            completions=types.SimpleNamespace(create=None))

    def _create(**kw):
        msgs = kw.get("messages") or []
        state["requests"].append(msgs)
        # One answer per user message (the wake-up): no tool round needed,
        # this test is about what the model was READ, not what it did.
        return _Stream(_answer("I see the job's outcome."))

    _Stub.chat.completions.create = _create

    client = ac.OpenAIClient.__new__(ac.OpenAIClient)
    client.client = _Stub()
    client.model = "stub-model"
    client._provider = "kit"
    client._base_url = "https://example.invalid/api/v1"
    client._api_key = "x"
    client.effort = ""
    client._permissions = KitToolPermissions(workspace=ws, mode="default",
                                             task_session_id=SESSION)
    client.on_model_switched = None
    client._steer_lock = threading.Lock()
    client._steer_msgs = []
    client._run_notes = []
    client._test_evidence = []
    client._observed_files_session = set()
    client._red_test_files = []
    client._thrash = None
    client.max_tool_rounds = 2

    list(client.stream_message(
        messages=[{"role": "user", "content": notice}],
        system="You are a test.", max_tokens=64))
    assert state["requests"], "the model was never asked anything"
    return state["requests"][-1]


def test_the_notice_reaches_the_model_in_the_words_that_were_rendered(
        job_registry, workspace):
    """The notice is what the next turn is given, and what the model is
    shown is that notice -- one job, named once."""
    job_id, exit_code = _run_the_calls(workspace)
    notice = _notice(workspace)
    assert notice, "the finished job produced no wake-up at all"

    shown = _model_turn_that_sees(notice, workspace)

    text = "\n".join(str(m.get("content") or "") for m in shown
                     if isinstance(m, dict))
    rows_in_context = [ln for ln in text.splitlines()
                       if ln.lstrip().startswith("- ") and job_id in ln]
    assert rows_in_context, (
        "the model's context does not carry the job's notice row: "
        + repr(text[:400]))
    assert len(rows_in_context) == 1, (
        "the model's context carries " + str(len(rows_in_context))
        + " rows for one job: " + repr(rows_in_context))

