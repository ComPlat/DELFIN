"""Keep session for a session in the dashboard's kernel
(delfin/dashboard/kept_agent_sessions.py).

Switched on, the session is resumed in the agent service once its kernel
has ended, kept there; switched off, it ends with the kernel. Judged by
the kernel's process (pid and start time), never by a widget. An
emergency stop given after the switch wins. A failed hand-over keeps the
record, with the reason, for the next caller.
"""
from __future__ import annotations

import os
import subprocess
import sys
import time

import pytest

from delfin.dashboard import kept_agent_sessions as K


@pytest.fixture(autouse=True)
def records_here(tmp_path, monkeypatch):
    monkeypatch.setenv(K.DIR_ENV, str(tmp_path / "kept"))
    monkeypatch.setattr(K, "_stopped_after", lambda armed_at: False)
    return tmp_path / "kept"


class _Service:
    """What hand_over needs of a LiveClient."""

    def __init__(self, running=(), fail: str = ""):
        self.created: list[tuple] = []
        self.running = list(running)
        self.fail = fail
        self.closed = False

    def sessions(self):
        return {"sessions": [{"sessionId": s} for s in self.running]}

    def create(self, cwd, **settings):
        if self.fail:
            raise RuntimeError(self.fail)
        self.created.append((cwd, settings))
        return {"sessionId": settings.get("resume")}

    def close(self):
        self.closed = True


def _ended_process() -> int:
    proc = subprocess.Popen([sys.executable, "-c", "pass"])
    proc.wait()
    return proc.pid


def test_arm_writes_an_owner_only_record_and_disarm_removes_it(records_here):
    path = K.arm("abc123", "/work", {"model": "m", "provider": "", "bogus": 1})
    assert K.is_armed("abc123")
    assert oct(path.stat().st_mode & 0o777) == "0o600"
    assert oct(records_here.stat().st_mode & 0o777) == "0o700"
    (row,) = K.records()
    assert (row["session_id"], row["workspace"], row["pid"]) == ("abc123", "/work", os.getpid())
    assert row["settings"] == {"model": "m"}, "only start settings, no empty ones"
    assert K.disarm("abc123") and not K.is_armed("abc123")


def test_a_session_whose_kernel_runs_is_held():
    K.arm("live", "/work")                      # held by this very process
    service = _Service()
    assert K.hand_over(connect=lambda: service) == [
        {"session_id": "live", "result": "held"}]
    assert service.created == [] and K.is_armed("live")


def test_once_its_kernel_ended_it_is_resumed_kept_in_the_service():
    K.arm("gone1", "/work", {"permissionMode": "acceptEdits"}, pid=_ended_process())
    service = _Service()
    rows = K.hand_over(connect=lambda: service)
    assert rows == [{"session_id": "gone1", "result": "resumed"}]
    assert service.created == [("/work", {"resume": "gone1", "keep": True,
                                          "permissionMode": "acceptEdits"})]
    assert not K.is_armed("gone1") and service.closed


def test_one_the_service_runs_already_is_not_started_twice():
    K.arm("twice", "/work", pid=_ended_process())
    service = _Service(running=["twice"])
    assert K.hand_over(connect=lambda: service)[0]["result"] == "resumed"
    assert service.created == []


def test_a_failure_keeps_the_record_with_its_reason():
    K.arm("later", "/work", pid=_ended_process())
    rows = K.hand_over(connect=lambda: _Service(fail="no service here"))
    assert rows[0]["result"].startswith("failed: RuntimeError: no service here")
    (row,) = K.records()
    assert row["session_id"] == "later" and "no service here" in row["last_error"]


def test_an_emergency_stop_after_the_switch_drops_it(monkeypatch):
    monkeypatch.setattr(K, "_stopped_after", lambda armed_at: True)
    K.arm("stopped", "/work", pid=_ended_process())
    service = _Service()
    rows = K.hand_over(connect=lambda: service)
    assert rows[0]["result"].startswith("dropped")
    assert service.created == [] and K.records() == []


def test_a_record_of_another_host_is_left_alone(monkeypatch):
    K.arm("far", "/work", pid=_ended_process())
    monkeypatch.setattr(K, "_host", lambda: "another-host")
    assert K.hand_over(connect=lambda: _Service())[0]["result"] == "held"


def test_a_reused_pid_is_not_the_kernel():
    """Same pid, other start time: the kernel is gone."""
    path = K.arm("reused", "/work")
    import json
    data = json.loads(path.read_text())
    data["start_ticks"] = (data["start_ticks"] or 0) + 12345
    path.write_text(json.dumps(data))
    assert K.holder_gone(K.records()[0])


def test_the_hand_over_waits_for_the_kernel_it_was_started_from():
    proc = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(1.5)"])
    K.arm("waits", "/work", pid=proc.pid)
    service = _Service()
    t0 = time.monotonic()
    rows = K.hand_over(wait_pid=proc.pid, wait_s=30, connect=lambda: service)
    proc.wait()
    assert time.monotonic() - t0 >= 1.0
    assert rows == [{"session_id": "waits", "result": "resumed"}]


def test_a_kernel_ending_spawns_the_hand_over_only_for_its_own(monkeypatch):
    spawned = []
    monkeypatch.setattr(K, "spawn_hand_over", lambda **kw: spawned.append(kw) or "proc")
    assert K.on_kernel_end() is None, "nothing armed here"
    K.arm("mine", "/work")
    assert K.on_kernel_end(emergency=True) is None, "an emergency stop hands nothing over"
    assert K.on_kernel_end() == "proc"
    assert spawned == [{"wait_pid": os.getpid()}]
