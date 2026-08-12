"""Durable attention inbox (delfin.agent.attention) + kit_confirm hooks.

Covers: durable emit that survives a simulated restart, list_pending
filtering, resolve + drain-exactly-once semantics, the notify_command
transport (invoked detached, failures swallowed), and the confirm-broker
integration (enqueue emits; a timeout PARKS the event instead of
resolving it).
"""

from __future__ import annotations

import importlib
import json
import threading
import time
from pathlib import Path

import pytest

import delfin.agent.attention as attention
from delfin.agent.kit_confirm import KitConfirmBroker


@pytest.fixture(autouse=True)
def _fake_home(tmp_path, monkeypatch):
    """Isolate ~ (inbox + settings) and mute the desktop/webhook transports."""
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setattr(Path, "home", lambda: home)
    # Desktop notifications must never fire from the test suite.
    import delfin.agent.notify as notify
    monkeypatch.setattr(notify, "send_notification", lambda *a, **k: True)
    monkeypatch.setattr(
        notify, "send_remote_trigger",
        lambda *a, **k: notify.TriggerResult(sent=False, error="test"))
    return home


def _inbox(home: Path) -> Path:
    return home / ".delfin" / "attention_inbox.jsonl"


# ---------------------------------------------------------------------------
# emit / persistence
# ---------------------------------------------------------------------------

def test_emit_persists_and_survives_restart(_fake_home):
    eid = attention.emit_attention(
        "question_pending",
        session_id="s1",
        title="Which basis set?",
        detail="def2-SVP or def2-TZVP",
        options=["def2-SVP", "def2-TZVP"],
        workspace="/tmp/ws",
    )
    assert eid.startswith("att-")
    assert _inbox(_fake_home).exists()

    # Simulated restart: a fresh module (no in-process state) still sees it.
    fresh = importlib.reload(attention)
    pending = fresh.list_pending()
    assert len(pending) == 1
    ev = pending[0]
    assert ev["id"] == eid
    assert ev["kind"] == "question_pending"
    assert ev["session_id"] == "s1"
    assert ev["title"] == "Which basis set?"
    assert ev["options"] == ["def2-SVP", "def2-TZVP"]
    assert ev["status"] == "pending"
    assert ev["acknowledged"] is False


def test_emit_rejects_unknown_kind(_fake_home):
    with pytest.raises(ValueError):
        attention.emit_attention("nonsense", title="x")
    assert attention.list_pending() == []


def test_load_tolerates_torn_lines(_fake_home):
    eid = attention.emit_attention("run_finished", title="done")
    p = _inbox(_fake_home)
    with p.open("a", encoding="utf-8") as fh:
        fh.write('{"torn": \n')            # crash artifact / foreign line
    pending = attention.list_pending()
    assert [ev["id"] for ev in pending] == [eid]


def test_list_pending_filters_by_kind(_fake_home):
    a = attention.emit_attention("confirm_pending", title="confirm rm")
    b = attention.emit_attention("run_finished", title="scheduled run done")
    assert {ev["id"] for ev in attention.list_pending()} == {a, b}
    only = attention.list_pending("run_finished")
    assert [ev["id"] for ev in only] == [b]
    assert attention.list_pending("plan_pending") == []


# ---------------------------------------------------------------------------
# resolve + drain-once
# ---------------------------------------------------------------------------

def test_resolve_and_drain_exactly_once_across_restart(_fake_home):
    eid = attention.emit_attention(
        "question_pending", session_id="s1", title="Pick one",
        options=["A", "B"])
    assert attention.resolve(eid, answer="A") is True
    assert attention.resolve(eid, answer="A") is False     # no longer pending
    assert attention.list_pending() == []

    # Restart between resolve and drain: acknowledged flag is file-backed.
    fresh = importlib.reload(attention)
    drained = fresh.drain_resolved("s1")
    assert len(drained) == 1
    assert drained[0]["id"] == eid
    assert drained[0]["answer"] == "A"
    assert fresh.drain_resolved("s1") == []                # exactly once
    fresh2 = importlib.reload(fresh)
    assert fresh2.drain_resolved("s1") == []               # still once


def test_drain_ignores_the_session_id_the_event_was_parked_with(_fake_home):
    """The session id is not stable, so it must not gate delivery.

    It is re-minted on a new cycle, cleared when the CLI backend restarts
    its process, and overwritten mid-turn by the backend's own id. Exact
    matching meant the user answered, the agent never saw it, and the
    record was invisible to every surface."""
    parked = attention.emit_attention(
        "question_pending", session_id="session-of-the-day", title="parked")
    unrouted = attention.emit_attention(
        "confirm_pending", title="unrouted confirm")
    attention.resolve(parked, answer="yes")
    attention.resolve(unrouted, answer="approved")
    got = attention.drain_resolved("a-completely-different-id")
    assert [ev["id"] for ev in got] == [parked, unrouted]
    assert attention.drain_resolved("a-completely-different-id") == []


def test_drain_can_be_scoped_by_workspace_and_kind(_fake_home):
    mine = attention.emit_attention(
        "question_pending", title="mine", workspace="/ws/a")
    theirs = attention.emit_attention(
        "question_pending", title="theirs", workspace="/ws/b")
    homeless = attention.emit_attention("confirm_pending", title="no ws")
    for eid in (mine, theirs, homeless):
        attention.resolve(eid, answer="ok")
    # An event without a workspace belongs to whoever drains first.
    got = attention.drain_resolved(workspace="/ws/a")
    assert [ev["id"] for ev in got] == [mine, homeless]
    got = attention.drain_resolved(workspace="/ws/b",
                                   kinds=("run_finished",))
    assert got == []                                       # kind filtered
    got = attention.drain_resolved(workspace="/ws/b")
    assert [ev["id"] for ev in got] == [theirs]


def test_drain_skips_pending_and_resolve_missing_is_false(_fake_home):
    attention.emit_attention("plan_pending", session_id="s1", title="plan")
    assert attention.drain_resolved("s1") == []            # pending ≠ answered
    assert attention.resolve("att-does-not-exist", answer="x") is False


def test_resolve_with_acknowledged_skips_drain(_fake_home):
    eid = attention.emit_attention("confirm_pending", title="live click")
    assert attention.resolve(eid, answer="approved", acknowledged=True)
    # Decision was delivered in-band — never replayed to the model.
    assert attention.drain_resolved("") == []


# ---------------------------------------------------------------------------
# notify_command transport
# ---------------------------------------------------------------------------

def _write_settings(home: Path, notify_command: str) -> None:
    (home / ".delfin_settings.json").write_text(json.dumps({
        "agent": {"attention": {"notify_command": notify_command}},
    }), encoding="utf-8")


def test_notify_command_invoked_with_title_argv(_fake_home, monkeypatch):
    _write_settings(_fake_home, "/usr/bin/my-notifier --urgent")
    calls: list[tuple] = []

    def _fake_run(argv, **kwargs):
        # Runs on the detached worker thread (exceptions are swallowed
        # there by design) — record everything, assert on the main thread.
        calls.append((list(argv), dict(kwargs)))
        return None

    monkeypatch.setattr(attention.subprocess, "run", _fake_run)
    attention.emit_attention("budget_warning", title="cost at 80%")
    t = attention._LAST_NOTIFY_THREAD
    assert t is not None
    t.join(timeout=5)
    assert len(calls) == 1
    argv, kwargs = calls[0]
    assert argv == ["/usr/bin/my-notifier", "--urgent", "cost at 80%"]
    assert kwargs.get("timeout") == 5
    assert "shell" not in kwargs                           # never shell=True
    assert kwargs.get("start_new_session") is True


def test_notify_command_failure_never_raises(_fake_home, monkeypatch):
    _write_settings(_fake_home, "/does/not/exist")

    def _boom(*_a, **_k):
        raise OSError("spawn failed")

    monkeypatch.setattr(attention.subprocess, "run", _boom)
    eid = attention.emit_attention("run_failed", title="entry x failed")
    t = attention._LAST_NOTIFY_THREAD
    if t is not None:
        t.join(timeout=5)
    # Event was still durably recorded despite the broken hook.
    assert [ev["id"] for ev in attention.list_pending()] == [eid]


def test_unparseable_notify_command_is_ignored(_fake_home, monkeypatch):
    _write_settings(_fake_home, 'broken "quote')

    def _never(*_a, **_k):                                 # must not be hit
        raise AssertionError("subprocess.run called for unparseable command")

    monkeypatch.setattr(attention.subprocess, "run", _never)
    attention.emit_attention("run_finished", title="ok")
    assert len(attention.list_pending()) == 1


# ---------------------------------------------------------------------------
# kit_confirm integration
# ---------------------------------------------------------------------------

def test_confirm_enqueue_emits_and_click_resolves(_fake_home):
    broker = KitConfirmBroker(default_timeout_s=10.0)
    out = {}

    def _worker():
        out["ok"] = broker.callback(
            "bash", {"command": "rm -rf build"}, "preview text")

    t = threading.Thread(target=_worker, daemon=True)
    t.start()
    for _ in range(200):                                   # wait for enqueue
        if attention.list_pending("confirm_pending"):
            break
        time.sleep(0.01)
    pending = attention.list_pending("confirm_pending")
    assert len(pending) == 1
    assert "bash" in pending[0]["title"]
    assert "preview text" in pending[0]["detail"]

    req = list(broker._pending)[0]
    with broker._lock:
        req.decision = True
    req.event.set()
    t.join(timeout=5)
    assert out["ok"] is True
    assert broker.last_timed_out is False
    # A real click resolves the event (acknowledged: answered in-band).
    assert attention.list_pending("confirm_pending") == []
    assert attention.drain_resolved("") == []
    all_events = attention._load_events()
    assert all_events[0]["status"] == "resolved"
    assert all_events[0]["answer"] == "approved"
    assert all_events[0]["acknowledged"] is True


def test_confirm_timeout_parks_event_pending(_fake_home):
    broker = KitConfirmBroker(default_timeout_s=0.05)
    ok = broker.callback("bash", {"command": "sbatch job.sh"}, "")
    assert ok is False
    assert broker.last_timed_out is True
    # Timeout must NOT resolve the event — it stays pending (parked) so
    # the user can still act on it later from the attention inbox.
    parked = attention.list_pending("confirm_pending")
    assert len(parked) == 1
    assert "sbatch job.sh" in parked[0]["detail"]
    assert attention.drain_resolved("") == []


def test_confirm_flow_survives_broken_inbox(_fake_home, monkeypatch):
    """emit_attention failing must never break the confirm dialog."""
    import delfin.agent.kit_confirm  # noqa: F401  (module under test)

    def _boom(*_a, **_k):
        raise OSError("disk full")

    monkeypatch.setattr(attention, "emit_attention", _boom)
    broker = KitConfirmBroker(default_timeout_s=0.05)
    ok = broker.callback("bash", {"command": "true"}, "")
    assert ok is False                                     # timeout, no crash
