"""The messenger gateway: only the paired contact commands, with limits.

* a message from another contact is refused, whatever it says;
* the same message id acts once, its replay is refused;
* the rate limit turns messages away; a message over the size limit is
  refused;
* status/sessions/pending/stop/hilfe reach the bridge and the answer is
  rendered; an approval in pending is named with where to approve;
* a number chooses the question, the next message answers it; free text
  answers the first open question; with nothing open it becomes a task
  for the phone session; @name routes to the named session;
* the bridge being down is a readable reply, not a crash;
* only the paired contact gets replies (Listener over a fake daemon).
"""

from __future__ import annotations

import asyncio

import pytest

from delfin.messenger import gateway as gw
from delfin.messenger.bridge_client import BridgeDown
from delfin.messenger.gateway import Gateway, PhoneSession, Refused


class FakeBridge:
    """Stands in for the bridge socket: records what was asked."""

    def __init__(self, replies=None):
        self.asked: list[dict] = []
        self.replies = replies or {}
        self.down = False

    def ask(self, command: dict) -> dict:
        if self.down:
            raise BridgeDown("socket gone")
        self.asked.append(command)
        cmd = command.get("cmd")
        if cmd in self.replies:
            return self.replies[cmd]
        if cmd == "status":
            return {"ok": True, "sessions": [
                {"id": "s1", "title": "Benzol", "state": "running",
                 "pending": []}]}
        if cmd == "sessions":
            return {"ok": True, "sessions": [
                {"id": "s1", "name": "Benzol", "state": "running"}]}
        if cmd == "pending":
            return {"ok": True, "items": [
                {"id": "q1", "kind": "question", "session": "Benzol",
                 "text": "welches Loesungsmittel?"},
                {"id": "c1", "kind": "confirm", "session": "Benzol",
                 "text": "may run orca"}]}
        if cmd == "stop":
            return {"ok": True, "stopping": True, "sessions": 2}
        if cmd == "prompt":
            return {"ok": True, "session": "phone-1", "accepted": True}
        if cmd == "answer":
            if command.get("id") == "c1":
                return {"error": "approvals are never answered from the phone; "
                                 "approve at the dashboard or terminal"}
            return {"ok": True}
        return {"error": f"no command {cmd!r}"}


def make(contact_id: str = "phone-contact"):
    bridge = FakeBridge()
    g = Gateway(contact_id=contact_id, bridge=bridge)
    return g, bridge


def msg(text: str, cid: str = "phone-contact", mid: str = "m1") -> PhoneSession:
    return PhoneSession(cid, text, mid)


class TestSenderGate:
    def test_a_message_from_another_contact_is_refused(self):
        async def run():
            g, bridge = make()
            with pytest.raises(Refused):
                await g.handle(msg("status", cid="stranger"))
            assert bridge.asked == []
        asyncio.run(run())

    def test_the_paired_contact_is_served(self):
        async def run():
            g, _ = make()
            reply = await g.handle(msg("status"))
            assert "Benzol" in reply
        asyncio.run(run())


class TestLimits:
    def test_the_same_message_id_acts_once(self):
        async def run():
            g, bridge = make()
            await g.handle(msg("stop", mid="m7"))
            with pytest.raises(Refused):
                await g.handle(msg("stop", mid="m7"))
            stops = [a for a in bridge.asked if a.get("cmd") == "stop"]
            assert len(stops) == 1
        asyncio.run(run())

    def test_the_rate_limit_turns_flood_away(self):
        async def run():
            g, _ = make()
            g.now = lambda: 0.0
            for i in range(gw.RATE_LIMIT):
                await g.handle(msg("status", mid=f"r{i}"))
            with pytest.raises(Refused):
                await g.handle(msg("status", mid="flood"))
        asyncio.run(run())

    def test_the_rate_limit_relaxes_after_the_window(self):
        async def run():
            g, _ = make()
            now = [0.0]
            g.now = lambda: now[0]
            for i in range(gw.RATE_LIMIT):
                await g.handle(msg("status", mid=f"w{i}"))
            now[0] = gw.RATE_WINDOW_S + 1.0
            reply = await g.handle(msg("status", mid="later"))
            assert "Benzol" in reply
        asyncio.run(run())

    def test_an_oversized_message_is_refused(self):
        async def run():
            g, bridge = make()
            with pytest.raises(Refused):
                await g.handle(msg("x" * (gw.SIZE_LIMIT + 1), mid="big"))
            assert bridge.asked == []
        asyncio.run(run())


class TestCommands:
    def test_help_says_where_approvals_live(self):
        async def run():
            g, _ = make()
            reply = await g.handle(msg("hilfe"))
            assert "Dashboard" in reply or "dashboard" in reply
            reply = await g.handle(msg("help", mid="m2"))
            assert "Dashboard" in reply or "dashboard" in reply
        asyncio.run(run())

    def test_sessions_in_german_and_english(self):
        async def run():
            g, bridge = make()
            await g.handle(msg("sitzungen", mid="m2"))
            await g.handle(msg("sessions", mid="m3"))
            assert all(a.get("cmd") == "sessions" for a in bridge.asked)
        asyncio.run(run())

    def test_pending_names_the_approval_with_where_to_approve(self):
        async def run():
            g, _ = make()
            reply = await g.handle(msg("pending", mid="m2"))
            assert "1. Frage" in reply
            assert "Freigabe nötig" in reply
            assert "Dashboard/Terminal" in reply
            assert "orca" in reply
        asyncio.run(run())

    def test_stop_reaches_the_bridge(self):
        async def run():
            g, bridge = make()
            reply = await g.handle(msg("stop", mid="m2"))
            assert bridge.asked[-1] == {"cmd": "stop"}
            assert "angehalten" in reply
        asyncio.run(run())

    def test_at_name_routes_to_the_named_session(self):
        async def run():
            g, bridge = make()
            await g.handle(msg("@Benzol starte die Frequenzrechnung", mid="m2"))
            prompt = bridge.asked[-1]
            assert prompt["cmd"] == "prompt"
            assert prompt["session"] == "Benzol"
            assert "Frequenzrechnung" in prompt["text"]
        asyncio.run(run())

    def test_free_text_becomes_a_phone_session_task(self):
        async def run():
            g, bridge = make()
            await g.handle(msg("rechne die Hessian von Benzol", mid="m2"))
            prompt = bridge.asked[-1]
            assert prompt["cmd"] == "prompt"
            assert "session" not in prompt      # the bridge makes the session
            assert "Hessian" in prompt["text"]
        asyncio.run(run())


class TestAnswers:
    def test_a_number_selects_the_question_then_text_answers_it(self):
        async def run():
            g, bridge = make()
            await g.handle(msg("pending", mid="m2"))
            reply = await g.handle(msg("1", mid="m3"))
            assert "Ausgewählt" in reply
            reply = await g.handle(msg("DMF", mid="m4"))
            answer = bridge.asked[-1]
            assert answer["cmd"] == "answer"
            assert answer["id"] == "q1"
            assert answer["text"] == "DMF"
            assert "weitergegeben" in reply
        asyncio.run(run())

    def test_an_answer_to_an_approval_is_refused_verbatim(self):
        async def run():
            g, bridge = make()
            await g.handle(msg("pending", mid="m2"))
            await g.handle(msg("2", mid="m3"))
            reply = await g.handle(msg("yes", mid="m4"))
            assert "never answered from the phone" in reply
            assert "Antwort weitergegeben" not in reply
        asyncio.run(run())

    def test_text_answers_the_first_open_question(self):
        async def run():
            g, bridge = make()
            await g.handle(msg("pending", mid="m2"))
            await g.handle(msg("DMF bitte", mid="m3"))
            answer = bridge.asked[-1]
            assert answer["cmd"] == "answer"
            assert answer["id"] == "q1"
        asyncio.run(run())

    def test_unknown_number_is_a_task_not_an_answer(self):
        async def run():
            g, bridge = make()
            await g.handle(msg("pending", mid="m2"))
            await g.handle(msg("42", mid="m3"))
            assert bridge.asked[-1]["cmd"] == "prompt"
        asyncio.run(run())


class TestBridgeDown:
    def test_a_down_bridge_is_a_readable_reply(self):
        async def run():
            g, bridge = make()
            bridge.down = True
            for text, mid in (("status", "m2"), ("stop", "m3"),
                              ("rechne was", "m4")):
                reply = await g.handle(msg(text, mid=mid))
                assert "nicht erreichbar" in reply
        asyncio.run(run())


class TestListener:
    def test_only_the_paired_contact_gets_replies(self):
        sent: list[tuple[str, str]] = []

        class FakeDaemon:
            async def subscribe(self):
                pass

            async def events(self):
                yield {"chatItems": [{"chatItem": {
                    "contact": {"contactId": "phone-contact"},
                    "msgContent": {"type": "text", "text": "status"}}}]}
                yield {"chatItems": [{"chatItem": {
                    "contact": {"contactId": "stranger"},
                    "msgContent": {"type": "text", "text": "stop"}}}]}
                yield {"chatItems": []}

            async def send_text(self, cid, text):
                sent.append((cid, text))

        async def run():
            g, _ = make()
            listener = gw.Listener(contact_id="phone-contact", gateway=g)
            await listener.run(FakeDaemon())

        asyncio.run(run())
        assert len(sent) == 1
        assert sent[0][0] == "phone-contact"
        assert "Benzol" in sent[0][1]


def test_trim_cuts_with_an_ellipsis():
    out = gw._trim("x" * 5_000, 100)
    assert len(out) == 100
    assert out.endswith("…")
