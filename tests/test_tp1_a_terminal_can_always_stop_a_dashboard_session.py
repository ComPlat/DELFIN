"""The terminal can always stop a session the dashboard hosts.

Input: a dashboard-hosted session whose answer ended in "?" (a greeting),
and a request that queues for a slot of the endpoint.

Before (2026-11-16, live): the dashboard's guess that the greeting was a
question went to the terminal as an open question; it took Enter and Esc,
so the running turn could not be stopped from there. Stop did not end a
request still waiting for a slot of the endpoint, and "waiting for a slot
... 1 in use" was said although a slot had been free at once.
Now a guessed question without choices stays in the dashboard; a wait for
a slot ends when the caller says stop; the notice is said only after a
real wait.
"""
from __future__ import annotations

import json
import threading
import time

import pytest

from delfin.agent import endpoint_gate as G
from delfin.dashboard.hosted_sessions import shared_questions

URL = "https://endpoint.example/v1"


@pytest.fixture
def runtime(monkeypatch, tmp_path):
    """A private runtime directory, as rundir resolves it."""
    rt = tmp_path / "rt"
    rt.mkdir()
    monkeypatch.setenv("DELFIN_RUNTIME_DIR", str(rt))
    monkeypatch.delenv("DELFIN_ENDPOINT_MAX_CONCURRENCY", raising=False)
    monkeypatch.delenv("DELFIN_ENDPOINT_PER_MODEL_LIMIT", raising=False)
    monkeypatch.setattr("delfin.user_settings.load_settings", lambda: {})
    return rt


def _gate():
    g = G.EndpointGate(URL, "m")
    if g._dir() is None:
        pytest.skip("no gate directory here")
    return g


def test_a_guessed_question_without_choices_stays_in_the_dashboard():
    pending = [
        {"requestId": "a", "kind": "ask", "question": "Hallo! Wie kann ich helfen?",
         "choices": [], "detected": True},
        {"requestId": "b", "kind": "ask", "question": "Which one?", "choices": ["1: x", "2: y"],
         "detected": True},
        {"requestId": "c", "kind": "ask", "question": "Real question", "choices": []},
        {"requestId": "d", "kind": "confirm", "tool": "bash"},
    ]
    assert [p["requestId"] for p in shared_questions(pending)] == ["b", "c", "d"]
    assert shared_questions([]) == [] and shared_questions(None) == []


def test_stop_ends_the_wait_for_a_slot(runtime):
    g = _gate()
    (g._dir() / "state.json").write_text(json.dumps({"limit": 1}))
    first = G.EndpointGate(URL, "m")
    assert first.acquire() == 0
    try:
        second = G.EndpointGate(URL, "m")
        stop = threading.Event()
        threading.Timer(0.3, stop.set).start()
        t0 = time.monotonic()
        with pytest.raises(G.SlotWaitCancelled):
            second.acquire(timeout_s=30.0, should_stop=stop.is_set)
        assert time.monotonic() - t0 < 5.0, "ended by the stop, not the deadline"
        assert getattr(second, "_held_fd", None) is None
        assert getattr(second, "_cap_fd", None) is None
    finally:
        first.release()


def test_a_slot_free_at_once_is_no_wait_and_says_nothing(runtime):
    """The subject is the notice, not which slot was handed out.

    This asserted slot 0 for the first caller and slot 1 for the second.
    That was true until #694 reserved slot-0 for a person: a headless
    request now scans 1..N-1 and never takes slot-0, so the first
    caller here gets 1 and the second 2. Pinning the index made this
    test fail on a deliberate change elsewhere while saying nothing
    about what it is for -- that a slot which was free at once must not
    be reported as a wait.
    """
    g = _gate()
    (g._dir() / "state.json").write_text(json.dumps({"limit": 4}))
    first = G.EndpointGate(URL, "m")
    first_slot = first.acquire()
    assert first_slot is not None and first_slot >= 1, (
        "slot-0 is the person's; a headless caller must not take it")
    try:
        second = G.EndpointGate(URL, "m")
        second_slot = second.acquire(timeout_s=1.0)
        assert second_slot is not None and second_slot != first_slot
        assert second.wait_notice == "", second.wait_notice
        second.release()
    finally:
        first.release()


def test_a_real_wait_is_still_said(runtime):
    g = _gate()
    (g._dir() / "state.json").write_text(json.dumps({"limit": 1}))
    first = G.EndpointGate(URL, "m")
    assert first.acquire() == 0
    second = G.EndpointGate(URL, "m")
    threading.Timer(0.3, first.release).start()
    assert second.acquire(timeout_s=10.0) == 0
    assert "waiting for a slot" in second.wait_notice
    second.release()


# -- the same suggestions in every window --------------------------------------------

def test_the_worker_offers_its_pending_tasks_as_suggestions():
    from types import SimpleNamespace
    from delfin.agent.service import protocol
    from delfin.agent.service.worker import Worker
    stub = SimpleNamespace(_SUGGESTIONS_MAX=Worker._SUGGESTIONS_MAX, _tasks=[
        {"num": 1, "status": "in_progress", "label": "Doing this"},
        {"num": 2, "status": "pending", "label": "Run the tests"},
        {"num": 3, "status": "completed", "label": "Done thing"},
        {"num": 4, "status": "pending", "label": "Write the summary"},
        {"num": 5, "status": "pending", "label": "  "},
    ])
    assert Worker._suggestions(stub) == ["Run the tests", "Write the summary"]
    assert "suggestions" in protocol.EVENTS["session.info"]
    assert "suggestions" in protocol.EVENTS["session.updated"]


def test_a_hosted_tab_offers_what_its_box_offers():
    from types import SimpleNamespace
    from delfin.dashboard.hosted_sessions import TabHost
    box = SimpleNamespace(placeholder="Run the tests", _dom_classes=("delfin-agent-input",
                                                                      "delfin-agent-input-suggests"))
    host = SimpleNamespace(adapter={"input": box})
    assert TabHost._suggestions(host) == ["Run the tests"]
    box._dom_classes = ("delfin-agent-input",)          # the resting hint
    assert TabHost._suggestions(host) == []
    assert TabHost._suggestions(SimpleNamespace(adapter={})) == []
