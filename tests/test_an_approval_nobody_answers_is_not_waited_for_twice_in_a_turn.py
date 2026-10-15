"""An approval nobody answers is not waited for twice in one turn.

Measured on the shared bug archive, 2026-09-16 .. 10-14: 21 approval
dialogs timed out at 300 s each (85 min in 6 sessions). Three of the full
waits came after another approval of the same kind had already expired
in the same turn; five times the very next call was the identical one,
although the timeout text said not to try it again now. The change the
agent meant to make was lost with the refusal: the user, back, had to
ask for it again.

  first approval of a kind expires     as before: the full window
  the next one of that kind, same turn "not now" at once, nobody waits
  another kind                         asked as usual
  next turn, or the user writes        asked again
  a write that was not approved        kept as a pending change the user
                                       applies with /approve <id>
  a real refusal                       unchanged: recorded as a denial
"""

from __future__ import annotations

import threading

import pytest

from delfin.agent import api_client as A
from delfin.agent import pending_changes as pc
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor


class _Away:
    """A dialog nobody answers: every request expires."""

    def __init__(self, timed_out=True):
        self.calls: list[str] = []
        self.last_timed_out = False
        self._timed_out = timed_out

    def callback(self, name, args, preview):
        self.calls.append(name)
        self.last_timed_out = self._timed_out
        return False


@pytest.fixture
def gate(tmp_path, monkeypatch):
    monkeypatch.setattr(pc, "_pending_root", lambda: tmp_path / "pending")
    ws = tmp_path / "ws"
    ws.mkdir()
    (ws / "a.txt").write_text("old a\n", encoding="utf-8")
    (ws / "b.txt").write_text("old b\n", encoding="utf-8")
    perms = KitToolPermissions(workspace=ws, mode="default")
    perms.task_session_id = "sess1"
    perms.turn_token = A._next_turn_token()
    broker = _Away()
    perms.confirm_callback = broker.callback
    eng = _DocToolExecutor.__new__(_DocToolExecutor)
    eng._permissions = perms

    def _write(path, content):
        return eng._run_permission_gate(
            "write_file", {"path": path, "content": content}, perms)
    return eng, perms, broker, _write


def test_the_second_write_in_the_turn_is_not_asked(gate):
    eng, perms, broker, write = gate
    first = write("a.txt", "new a\n")
    assert "TIMED OUT" in first
    second = write("b.txt", "new b\n")
    assert "TIMED OUT" in second
    assert broker.calls == ["write_file"], "the second write waited again"


def test_the_identical_retry_is_not_asked(gate):
    eng, perms, broker, write = gate
    write("a.txt", "new a\n")
    write("a.txt", "new a\n")
    assert len(broker.calls) == 1


def test_another_kind_is_still_asked(gate):
    eng, perms, broker, write = gate
    write("a.txt", "new a\n")
    eng._run_permission_gate(
        "bash", {"command": "make build"}, perms)
    assert broker.calls[0] == "write_file"
    assert len(broker.calls) == 2


def test_a_new_turn_asks_again(gate):
    eng, perms, broker, write = gate
    write("a.txt", "new a\n")
    perms.turn_token = A._next_turn_token()
    write("b.txt", "new b\n")
    assert len(broker.calls) == 2


def test_a_message_from_the_user_asks_again(gate):
    eng, perms, broker, write = gate
    write("a.txt", "new a\n")
    c = A.OpenAIClient.__new__(A.OpenAIClient)
    c._steer_lock = threading.Lock()
    c._steer_msgs = []
    c._permissions = perms
    c.push_steer(A._HarnessText("[watch] a job finished"))
    write("b.txt", "new b\n")
    assert len(broker.calls) == 1, "harness text is not the user coming back"
    c.push_steer("I am back, go on")
    write("b.txt", "newer b\n")
    assert len(broker.calls) == 2


def test_the_change_is_kept_for_the_user(gate):
    eng, perms, broker, write = gate
    first = write("a.txt", "new a\n")
    second = write("b.txt", "new b\n")
    pending = pc.list_pending("sess1")
    assert {p["path"].rsplit("/", 1)[-1] for p in pending} == {"a.txt", "b.txt"}
    ids = {p["path"].rsplit("/", 1)[-1]: p["id"] for p in pending}
    assert f"/approve {ids['a.txt']}" in first
    assert f"/approve {ids['b.txt']}" in second
    # Nothing was written: the user decides.
    assert (perms.workspace / "a.txt").read_text() == "old a\n"
    res = pc.approve("sess1", ids["a.txt"], workspace=perms.workspace)
    assert res["status"] == "applied"
    assert (perms.workspace / "a.txt").read_text() == "new a\n"


def test_a_protected_path_is_kept_with_its_warning(gate, monkeypatch):
    eng, perms, broker, write = gate
    perms.mode = "bypassPermissions"
    target = perms.workspace / "delfin" / "agent" / "engine.py"
    target.parent.mkdir(parents=True)
    target.write_text("x = 1\n", encoding="utf-8")
    perms.read_tracker[str(target.resolve())] = target.stat().st_mtime
    out = eng._run_permission_gate(
        "edit_file", {"path": "delfin/agent/engine.py",
                      "old_string": "x = 1", "new_string": "x = 2"}, perms)
    assert "TIMED OUT" in out
    pending = pc.list_pending("sess1")
    assert len(pending) == 1
    assert "SELF-MODIFICATION" in pending[0]["note"]
    assert target.read_text() == "x = 1\n"


def test_a_real_refusal_is_still_a_refusal(gate):
    eng, perms, broker, write = gate
    broker._timed_out = False
    out = write("a.txt", "new a\n")
    assert "denied" in out
    assert pc.list_pending("sess1") == []
    assert not perms.approval_away
