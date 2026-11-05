"""The phone is a contact channel: it hears questions, approvals, done and
failed, and what the agent deliberately sends it; nothing else.

Defaults (coordinator brief, real use on the Pi):

* no tool-call stream; a "working" message only for a turn the person
  started from the phone;
* ``details`` shows more on request (the current activity);
* push_notification reaches the phone when the session is connected,
  rate-limited (one per minute per session), and never when it is not.

Output, through the real gateway and bridge (the rig) and the real tool
executor, each with its control.
"""
from __future__ import annotations

import asyncio
import json

import pytest

from _messenger_people_rig import Account
from delfin.agent import notify
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor
from delfin.messenger import activation, people, store
from delfin.messenger.gateway import _PERSON


@pytest.fixture
def acct(tmp_path, monkeypatch):
    monkeypatch.setenv("DELFIN_MESSENGER_LANG", "en")
    monkeypatch.setattr(notify, "send_notification", lambda *a, **k: False)
    account = Account(tmp_path, monkeypatch)
    store.save("21", person="anna", label="anna phone")
    account.enable("anna", "s-anna")
    store.put_room("#a1", "s-anna", "Butan", person="anna", ready=True)
    account.service = account.services["anna"]
    account.row = account.service.sessions[0]
    account.tmp = tmp_path
    return account


def refresh(acct):
    async def go():
        token = _PERSON.set("anna")
        try:
            await acct.rooms.refresh()
        finally:
            _PERSON.reset(token)
    asyncio.run(go())


def _turn(acct, answer):
    acct.row["state"] = "running"
    refresh(acct)
    acct.row.update(state="idle", answer=answer)
    refresh(acct)


def test_a_turn_started_elsewhere_sends_no_working_message_but_its_answer(acct):
    refresh(acct)
    _turn(acct, "Barrier: 12.3 kcal/mol.")
    sent = acct.daemon.to("#a1")
    assert not [t for t in sent if "Working" in t], sent
    assert sent[-1] == "Barrier: 12.3 kcal/mol."


def test_a_turn_started_from_the_phone_sends_one_working_message(acct):
    refresh(acct)
    acct.say("21", "#a1", "compute the barrier")
    _turn(acct, "Barrier: 12.3 kcal/mol.")
    assert len([t for t in acct.daemon.to("#a1") if "Working" in t]) == 1


def test_details_shows_the_current_activity(acct):
    acct.row.update(state="running", activity="optimising conformer 2 of 5")
    acct.say("21", "#a1", "details")
    out = acct.daemon.to("#a1")[-1]
    assert "optimising conformer 2 of 5" in out
    assert not acct.service.calls("session.prompt"), "details is not a prompt"


def _push(acct, monkeypatch, sid="s-anna", title="Done", body="Barrier 12.3 kcal/mol"):
    monkeypatch.setattr(activation, "_DIR", people.activation_dir("anna"))
    perms = KitToolPermissions(workspace=str(acct.tmp), mode="acceptEdits")
    perms.task_session_id = sid
    out = _DocToolExecutor().execute("push_notification", {"title": title, "body": body}, perms)
    return json.loads(out)


def test_push_notification_reaches_the_connected_phone_once_a_minute(acct, monkeypatch):
    refresh(acct)
    first = _push(acct, monkeypatch)
    assert first["phone"] is True, first
    refresh(acct)
    note = acct.daemon.to("#a1")[-1]
    assert note == "🔔 *Done*\nBarrier 12.3 kcal/mol", note
    again = _push(acct, monkeypatch, title="Again")
    assert again["phone"] == "rate_limited"
    refresh(acct)
    assert not [t for t in acct.daemon.to("#a1") if "Again" in t]
    refresh(acct)
    assert acct.daemon.to("#a1").count(note) == 1, "a note is sent once"


def test_push_notification_never_reaches_a_phone_the_session_is_not_on(acct, monkeypatch):
    out = _push(acct, monkeypatch, sid="s-anna2")
    assert out["phone"] is False
    refresh(acct)
    assert not [t for chat, t in acct.daemon.sent if "Done" in t]
