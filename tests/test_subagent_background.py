"""Background subagents: sa_id round-trip + result collection.

subagent(background=true) now returns the reserved sa_id and forwards it to the
run, and subagent_result(sa_id) collects status/result — so a backgrounded
subagent is referenceable and its output retrievable (not just visible in the
panel).
"""

from __future__ import annotations

import json
import threading

import pytest

from delfin.agent import subagents as S
from delfin.agent import api_client as A


@pytest.fixture(autouse=True)
def _iso(monkeypatch, tmp_path):
    monkeypatch.setattr(S, "_RUNNING_PATH", tmp_path / "running.json")
    monkeypatch.setattr(S, "_SESSIONS_DIR", tmp_path / "sessions")
    monkeypatch.setattr(S, "_TELEMETRY_PATH", tmp_path / "telemetry.jsonl")
    yield


# ---------------------------------------------------------------------------
# get_subagent_result
# ---------------------------------------------------------------------------


def test_result_unknown():
    assert S.get_subagent_result("nope")["status"] == "unknown"


def test_result_requires_id():
    assert "error" in S.get_subagent_result("")


def test_result_running(monkeypatch):
    monkeypatch.setattr(S, "read_running",
                        lambda: {"abc": {"type": "explore",
                                         "description": "probe", "started_at": 1}})
    out = S.get_subagent_result("abc")
    assert out["status"] == "running"
    assert out["subagent_type"] == "explore"


def test_result_finished(monkeypatch):
    monkeypatch.setattr(S, "read_running", lambda: {})
    monkeypatch.setattr(S, "load_subagent_session", lambda sa: {
        "subagent_type": "explore", "description": "probe",
        "messages": [{"role": "user", "content": "go"},
                     {"role": "assistant", "content": "FOUND IT"}],
        "error": ""})
    out = S.get_subagent_result("xyz")
    assert out["status"] == "finished"
    assert out["final_text"] == "FOUND IT"


# ---------------------------------------------------------------------------
# subagent_result tool dispatch
# ---------------------------------------------------------------------------


def test_subagent_result_tool(monkeypatch):
    monkeypatch.setattr(S, "get_subagent_result",
                        lambda sa: {"sa_id": sa, "status": "finished",
                                    "final_text": "REPORT"})
    out = json.loads(A._doc_executor._execute_subagent_result({"sa_id": "k"}))
    assert out["final_text"] == "REPORT"


# ---------------------------------------------------------------------------
# Background launch returns + forwards the sa_id
# ---------------------------------------------------------------------------


def test_background_returns_and_forwards_sa_id():
    received: dict = {}
    done = threading.Event()

    def _runner(**kw):
        received.update(kw)
        done.set()
        return {"ok": True}

    class _Perms:
        subagent_runner = staticmethod(_runner)

    out = json.loads(A._doc_executor._execute_subagent(
        {"subagent_type": "explore", "description": "probe the thing",
         "prompt": "look around carefully for the bug in module x",
         "background": True},
        _Perms()))

    assert out["status"] == "started_in_background"
    assert out["sa_id"]
    assert done.wait(timeout=3)
    assert received.get("sa_id") == out["sa_id"]
