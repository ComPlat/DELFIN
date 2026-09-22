"""A refusal from the supervisor can carry its reason to the agent.

Seen while supervising, 2026-09-21: a session asked to read a file outside
its worktree. Approving would have opened the whole directory for the rest
of the session, so the supervisor refused -- and the agent was told only
"user denied". The reason, and what to do instead, had to travel as a
session message, which arrives at the idle prompt: after the turn that
needed it. An agent refused without a reason guesses, and the likeliest
guess is the same thing spelled differently -- one more dialog, and one
more refusal.

Claude Code lets a refusal say what to do differently, and the model reads
it in the same turn. Here the refusal from outside can carry a reason: it
is written into the answer file, read under the same checks as the answer
itself, kept per thread like the expiry flag the gate already reads, and
appended to the refusal the model receives. A reason never outlives its
refusal: the next dialog, or a command that asked nothing, carries none.
"""

from __future__ import annotations

import json
import threading
import time

import pytest

from delfin.agent import api_client as A
from delfin.agent import file_confirm as fc
from delfin.agent import terminal_confirm as tc

REASON = "Read it with read_file inside your worktree; ~/.delfin stays closed."


@pytest.fixture
def room(tmp_path, monkeypatch):
    monkeypatch.setattr(tc, "_PENDING_DIR", tmp_path / "pending")
    return tmp_path / "pending"


def _ask_in_thread(broker, box, tool="bash", command="rm -r build"):
    def run():
        box["decision"] = broker.callback(tool, {"command": command}, "x")
        box["reason"] = broker.last_refusal_reason
    t = threading.Thread(target=run, daemon=True)   # red must fail, not hang
    t.start()
    return t


def _wait_published(room, n=1, timeout=5.0):
    end = time.monotonic() + timeout
    while time.monotonic() < end:
        found = sorted(room.glob("*.request.json")) if room.exists() else []
        if len(found) >= n:
            return found[-1].name[:-len(".request.json")]
        time.sleep(0.02)
    raise AssertionError("the question was never published")


# ---------------------------------------------------------------------------
# The answer file
# ---------------------------------------------------------------------------

class TestTheAnswerFile:
    def _request(self, room, rid="1-a"):
        room.mkdir(parents=True, exist_ok=True)
        (room / f"{rid}.request.json").write_text(json.dumps({"id": rid}))
        return rid

    def test_a_refusal_keeps_its_reason(self, tmp_path):
        rid = self._request(tmp_path)
        assert fc.answer(rid, False, room=tmp_path, reason=REASON)
        rec = json.loads((tmp_path / f"{rid}.answer.json").read_text())
        assert rec["decision"] == fc.DENY and rec["reason"] == REASON

    def test_an_approval_carries_none(self, tmp_path):
        rid = self._request(tmp_path)
        assert fc.answer(rid, True, room=tmp_path, reason="ignored")
        rec = json.loads((tmp_path / f"{rid}.answer.json").read_text())
        assert "reason" not in rec

    def test_control_characters_do_not_reach_the_model(self, tmp_path):
        rid = self._request(tmp_path)
        fc.answer(rid, False, room=tmp_path,
                  reason="no\x1b[2J\r\nclear\x07 " + "x" * 2000)
        rec = json.loads((tmp_path / f"{rid}.answer.json").read_text())
        assert "\x1b" not in rec["reason"] and "\x07" not in rec["reason"]
        assert "\r" not in rec["reason"] and "\n" not in rec["reason"]
        assert len(rec["reason"]) <= fc.REASON_MAX


# ---------------------------------------------------------------------------
# The broker
# ---------------------------------------------------------------------------

class TestTheBroker:
    def test_the_asking_thread_learns_the_reason(self, room):
        broker = tc.TerminalConfirmBroker(session_key="s", poll_s=0.02, timeout_s=10)
        box: dict = {}
        t = _ask_in_thread(broker, box)
        rid = _wait_published(room)
        assert tc.answer_waiting(rid, False, by="op", reason=REASON)
        t.join(5)
        assert box["decision"] is False
        assert box["reason"] == REASON

    def test_the_next_dialog_does_not_inherit_it(self, room):
        broker = tc.TerminalConfirmBroker(session_key="s", poll_s=0.02, timeout_s=10)
        first: dict = {}
        t = _ask_in_thread(broker, first)
        tc.answer_waiting(_wait_published(room, 1), False, reason=REASON)
        t.join(5)
        second: dict = {}
        t = _ask_in_thread(broker, second, command="rm -r dist")
        tc.answer_waiting(_wait_published(room, 1), True)
        t.join(5)
        assert second["decision"] is True and second["reason"] == ""

    def test_another_thread_does_not_see_it(self, room):
        broker = tc.TerminalConfirmBroker(session_key="s", poll_s=0.02, timeout_s=10)
        box: dict = {}
        t = _ask_in_thread(broker, box)
        tc.answer_waiting(_wait_published(room), False, reason=REASON)
        t.join(5)
        assert broker.last_refusal_reason == ""      # this thread asked nothing


# ---------------------------------------------------------------------------
# What the model reads
# ---------------------------------------------------------------------------

class _Refuser:
    """Stands in for the broker: refuses, and says why."""

    def __init__(self, reason):
        self.last_refusal_reason = ""
        self.last_timed_out = False
        self._reason = reason
        self.asked = 0

    def callback(self, name, args, preview):
        self.asked += 1
        self.last_refusal_reason = self._reason
        return False


def test_the_refusal_the_model_reads_carries_the_reason(tmp_path):
    refuser = _Refuser(REASON)
    perms = A.KitToolPermissions(workspace=str(tmp_path), mode="default",
                                 confirm_callback=refuser.callback)
    err = A._doc_executor._run_permission_gate(
        "bash", {"command": "python3 -c 'import os; os.system(\"x\")'"}, perms)
    assert refuser.asked == 1
    assert err and REASON in err


def test_a_command_that_asked_nothing_carries_no_old_reason(tmp_path):
    refuser = _Refuser(REASON)
    refuser.last_refusal_reason = REASON            # left over from before
    perms = A.KitToolPermissions(workspace=str(tmp_path), mode="default",
                                 confirm_callback=refuser.callback)
    err = A._doc_executor._run_permission_gate(
        "bash", {"command": "rm -rf /"}, perms)       # refused by the deny list
    assert refuser.asked == 0
    assert err and REASON not in err


def test_the_cli_passes_the_reason(monkeypatch):
    from delfin.agent import cli
    seen = {}
    monkeypatch.setattr(fc, "answer", lambda *a, **k: False)
    monkeypatch.setattr(tc, "answer_waiting",
                        lambda rid, d, by="", reason="": seen.update(
                            rid=rid, d=d, reason=reason) or True)
    args = cli.build_parser().parse_args(
        ["approvals", "deny", "1-a", "--reason", REASON])
    assert cli.cmd_approvals(args) == 0
    assert seen == {"rid": "1-a", "d": False, "reason": REASON}
