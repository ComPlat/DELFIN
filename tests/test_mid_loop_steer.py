"""Mid-loop steering: the user can inject a message INTO the running tool loop
so the model reacts within the SAME turn (user 2026-06-26: "können wir das
Modell während des Prozesses noch Sachen schreiben und beeinflussen")."""
import threading
import pytest
from delfin.agent.api_client import OpenAIClient


def _client():
    # Build without network: skip the real OpenAI client construction.
    c = OpenAIClient.__new__(OpenAIClient)
    c._steer_lock = threading.Lock()
    c._steer_msgs = []
    return c


def test_push_and_drain_roundtrip():
    c = _client()
    assert c._drain_steer() == []
    c.push_steer("build it differently")
    c.push_steer("  ")          # blank ignored
    c.push_steer("also add tests")
    assert c._drain_steer() == ["build it differently", "also add tests"]
    assert c._drain_steer() == []   # drained once


def test_push_is_thread_safe():
    c = _client()
    def worker(n):
        for i in range(50):
            c.push_steer(f"{n}-{i}")
    ts = [threading.Thread(target=worker, args=(n,)) for n in range(4)]
    for t in ts: t.start()
    for t in ts: t.join()
    assert len(c._drain_steer()) == 200    # nothing dropped


def test_engine_steer_delegates():
    from delfin.agent.engine import AgentEngine
    eng = AgentEngine.__new__(AgentEngine)
    eng.client = _client()
    assert eng.steer("go left") is True
    assert eng.client._drain_steer() == ["go left"]
    # An engine whose client lacks push_steer (e.g. a CLI backend) returns False
    class _Dummy: ...
    eng.client = _Dummy()
    assert eng.steer("x") is False
