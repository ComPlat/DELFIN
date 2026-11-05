"""One answer closes a question or approval everywhere, from whichever side.

Field report (Pi, an owner session): an approval and an ask_user_question
answered on the phone stayed open in the session's window; one was still
pending after the phone had answered.

The one authority for a request is the broker of the process that asked
(TerminalConfirmBroker, in the agent service its ProtocolBroker). It has
two doors, and the phone uses the one the session has:

* a session in the agent service: permission.reply through the service --
  the broker resolves it and announces permission.resolved with the same
  request id to every window (TUI, dashboard live view);
* a session outside the service (the terminal REPL, a dashboard tab that
  runs its own engine): the request room ~/.delfin/terminal_confirmations,
  answered by the bridge with the session's one-time code, read by the
  waiting broker's own checked reader. A request of a service session is
  never offered through the room (one door per request).

Output, through the real gateway and bridge and a real broker
(tests/_messenger_people_rig.py), each with its control:

* a question and an approval of a service session answered on the phone:
  resolved, permission.resolved with the same id, nothing pending;
* answered in a window instead: the phone's card says "answered at the
  computer" (an update of the card), and a late phone reply frees nothing;
* a question and an approval of a session outside the service (room door)
  answered on the phone: the waiting broker returns the phone's answer.
"""
from __future__ import annotations

import asyncio
import threading
import time

import pytest

from _messenger_people_rig import Account
from delfin.agent import approval_code, terminal_confirm as tc
from delfin.agent.service import approvals
from delfin.messenger import activation, store


class BrokerService:
    """A person's agent service with one session whose questions a REAL
    ProtocolBroker asks: session.pending is what its permission.asked
    events announced, permission.reply is the worker's answer path."""

    def __init__(self, ws):
        self.events = []
        self.pending = {}
        self.broker = approvals.ProtocolBroker(self._emit, session_id="s-anna")
        self.ws = ws
        self.closed = False
        self.asked = []

    def _emit(self, event, **fields):
        self.events.append((event, fields))
        if event == "permission.asked":
            self.pending[fields["requestId"]] = fields
        elif event == "permission.resolved":
            self.pending.pop(fields["requestId"], None)

    def request(self, method, params=None, timeout=30.0):
        params = params or {}
        self.asked.append((method, params))
        if method == "sessions.list":
            return {"sessions": [{"sessionId": "s-anna", "title": "Butan",
                                  "workspace": self.ws, "state": "waiting"}], "elsewhere": []}
        if method == "session.pending":
            return {"pending": list(self.pending.values())}
        if method == "permission.reply":
            landed, message = self.broker.answer(
                params["requestId"], choice=params.get("choice", ""),
                answers=params.get("answers") or [], text=params.get("text", ""),
                persist=bool(params.get("persist")), client="delfin-messenger-bridge")
            return {"applied": landed, "message": message}
        return {"accepted": True}

    def close(self):
        self.closed = True


@pytest.fixture
def acct(tmp_path, monkeypatch):
    return make_account(tmp_path, monkeypatch)


def make_account(tmp_path, monkeypatch):
    """anna's account with one service session s-anna whose questions a
    real ProtocolBroker asks, enabled for her phone 21 in room #a1."""
    monkeypatch.setenv("DELFIN_MESSENGER_LANG", "en")
    account = Account(tmp_path, monkeypatch)
    anna_home = account.home / "people" / "anna"
    # Every broker of this test is anna's process: its room and codes are hers.
    monkeypatch.setattr(tc, "_PENDING_DIR", anna_home / ".delfin" / "terminal_confirmations")
    monkeypatch.setattr(approval_code, "codes_dir", lambda: anna_home / ".delfin" / "approval_codes")
    service = BrokerService(account.services["anna"].sessions[0]["workspace"])
    account.services["anna"] = service
    account.bridge._person_bridge("anna")._service = service
    store.save("21", person="anna", label="anna phone")
    account.enable("anna", "s-anna")
    store.put_room("#a1", "s-anna", "Butan", person="anna", ready=True)
    account.service = service
    return account


def _ask_in_thread(fn, *args):
    box = {}
    th = threading.Thread(target=lambda: box.update(result=fn(*args)), daemon=True)
    th.start()
    return th, box


def _wait_for(cond, seconds=5.0):
    deadline = time.monotonic() + seconds
    while time.monotonic() < deadline:
        if cond():
            return True
        time.sleep(0.02)
    return False


def _card(acct, chat):
    """The newest card the phone got in chat, and its item id."""
    mids = [k[1] for k in acct.rooms._reply_targets if k[0] == chat]
    return mids[-1] if mids else None


QUESTION = {"question": "Which solvent?", "options": [{"label": "water"}, {"label": "toluene"}],
            "multiSelect": False}


def test_a_phone_answer_closes_a_service_question_everywhere(acct):
    broker = acct.service.broker
    th, box = _ask_in_thread(broker.ask_user, QUESTION)
    assert _wait_for(lambda: acct.service.pending)
    rid = next(iter(acct.service.pending))
    acct.say("21", "#a1", "pending")
    mid = _card(acct, "#a1")
    assert mid
    acct.say("21", "#a1", "2", quote=mid)
    th.join(5)
    assert box["result"]["answers"] == ["toluene"]
    resolved = [f for k, f in acct.service.events if k == "permission.resolved"]
    assert [f["requestId"] for f in resolved] == [rid] and not broker.pending()
    assert not acct.service.pending


def test_a_phone_answer_closes_a_service_approval_everywhere(acct):
    broker = acct.service.broker
    th, box = _ask_in_thread(broker.callback, "write_file", {"path": "notes.txt"}, "notes.txt: +5")
    assert _wait_for(lambda: acct.service.pending)
    rid = next(iter(acct.service.pending))
    acct.say("21", "#a1", "pending")
    acct.say("21", "#a1", "ja", quote=_card(acct, "#a1"))
    th.join(5)
    assert box["result"] is True
    resolved = [f for k, f in acct.service.events if k == "permission.resolved"]
    assert [f["requestId"] for f in resolved] == [rid] and not broker.pending()


def test_answered_at_a_window_the_phone_card_says_so_and_a_late_reply_frees_nothing(acct):
    broker = acct.service.broker
    th, box = _ask_in_thread(broker.callback, "write_file", {"path": "notes.txt"}, "notes.txt: +5")
    assert _wait_for(lambda: acct.service.pending)
    rid = next(iter(acct.service.pending))
    acct.say("21", "#a1", "pending")
    mid = _card(acct, "#a1")
    broker.answer(rid, choice="y", client="delfin-tui")              # the TUI answers
    th.join(5)
    asyncio.run(_refresh(acct))
    assert any(chat == "#a1" and str(item) == str(mid) and "computer" in text
               for chat, item, text in acct.daemon.updated)
    before = len([1 for k, _ in acct.service.events if k == "permission.resolved"])
    acct.say("21", "#a1", "nein", quote=mid)                          # late
    assert len([1 for k, _ in acct.service.events if k == "permission.resolved"]) == before


async def _refresh(acct):
    from delfin.messenger.gateway import _PERSON
    token = _PERSON.set("anna")
    try:
        await acct.rooms.refresh()
        await acct.rooms.approvals_news()
    finally:
        _PERSON.reset(token)


def test_a_session_outside_the_service_is_answered_through_its_room(acct):
    # A terminal REPL session of anna's, enabled for the phone, not in her service.
    activation.configure("s-repl", directory=acct.activation_dir("anna"))
    repl = tc.TerminalConfirmBroker(session_id="s-repl", poll_s=0.02)
    th, box = _ask_in_thread(repl.callback, "write_file", {"path": "x.txt"}, "x.txt: +1")
    assert _wait_for(lambda: list(tc._PENDING_DIR.glob("*.request.json")))
    asyncio.run(_refresh(acct))
    card = [k for k in acct.rooms._reply_targets if k[0] == "@21"]
    assert card, acct.daemon.to("@21")
    acct.say("21", "@21", "ja", quote=card[-1][1])
    th.join(5)
    assert box["result"] is True
    # The question of the same session, answered with an option number.
    th2, box2 = _ask_in_thread(repl.ask_user, QUESTION)
    assert _wait_for(lambda: list(tc._PENDING_DIR.glob("*.request.json")))
    asyncio.run(_refresh(acct))
    card = [k for k in acct.rooms._reply_targets if k[0] == "@21"]
    acct.say("21", "@21", "1", quote=card[-1][1])
    th2.join(5)
    assert box2["result"]["answers"] == ["water"]


def test_a_service_sessions_request_is_never_offered_through_the_room(acct):
    broker = acct.service.broker
    th, _ = _ask_in_thread(broker.callback, "write_file", {"path": "n.txt"}, "n.txt: +1")
    assert _wait_for(lambda: acct.service.pending)
    reply = asyncio.run(acct.bridge._run_command({"cmd": "room", "person": "anna", "action": "pending"}))
    assert reply.get("items") == []
    broker.abort_all()
    th.join(5)
