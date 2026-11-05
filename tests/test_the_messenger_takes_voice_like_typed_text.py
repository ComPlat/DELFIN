"""A voice message is a prompt or an answer, and never approves.

Input: Rooms.voice (delfin/messenger/gateway.py), called with the local
transcript of a voice message, its chat and, when the person used the reply
gesture, the message it answers.

- Every reply starts with one line echoing the transcript
  ("🎙 understood: ..."), so the person sees what arrived.
- Without a reply reference the transcript goes to the session as its next
  prompt, marked as transcribed; while exactly one question waits, it is
  that question's answer.
- As a reply to a question the whole transcript is the answer.
- As a reply to an approval, a spoken yes or an offered option's label
  approves nothing and says to type it; a spoken no denies (refusing is
  never the risky direction). A sentence approves nothing.
- A lasting rule is never confirmed by voice; that stays typed.
"""
from __future__ import annotations

import asyncio

import pytest

from delfin.messenger import store
from delfin.messenger.gateway import Rooms

from test_the_messenger_gateway_only_serves_the_paired_phone import (
    Bridge, Daemon)


@pytest.fixture
def rig(tmp_path, monkeypatch):
    monkeypatch.setenv("DELFIN_MESSENGER_DIR", str(tmp_path))
    store.save("7")
    store.put_room("#10", "s1", "same title", owner="7", ready=True)
    daemon, bridge = Daemon(), Bridge()
    return Rooms(daemon, bridge), daemon, bridge


CONFIRM = {"id": "c1", "sessionId": "s1", "kind": "confirm",
           "text": "run calculation", "phoneApproval": True,
           "options": [{"key": "y", "label": "Once"}, {"key": "n", "label": "No"},
                       {"key": "A", "label": "Always", "detail": "run calc"}]}


def _request(rooms, item):
    async def run():
        await rooms.send_request("#10", store.rooms()["#10"], item)
        return next(mid for (chat, mid), it in rooms._reply_targets.items()
                    if it["id"] == item["id"])
    return asyncio.run(run())


def _voice(rooms, transcript, quote_id=None):
    return asyncio.run(rooms.voice(transcript, store.rooms()["#10"], "#10",
                                   quote_id))


def _decisions(bridge):
    return [c for c in bridge.calls if c.get("action") == "decide"]


def test_a_voice_message_without_a_reply_is_the_next_prompt(rig):
    rooms, _, bridge = rig
    out = _voice(rooms, "Optimiere bitte beide Konformere.")
    assert bridge.calls == [{"cmd": "session", "session": "s1",
                             "action": "prompt",
                             "text": "Voice message from the user (transcript):\n"
                                     "Optimiere bitte beide Konformere."}]
    assert out.splitlines()[0] == "🎙 verstanden: Optimiere bitte beide Konformere."


def test_the_echo_is_english_on_request(rig, monkeypatch):
    monkeypatch.setenv("DELFIN_MESSENGER_LANG", "en")
    rooms, _, _ = rig
    out = _voice(rooms, "Optimise both conformers.")
    assert out.splitlines()[0] == "🎙 understood: Optimise both conformers."


def test_without_a_reply_a_voice_message_answers_the_one_waiting_question(rig):
    rooms, _, bridge = rig
    bridge.pending["s1"] = [{"id": "q1", "sessionId": "s1", "kind": "question",
                             "text": "Which solvent?"}]
    _voice(rooms, "Toluol")
    answers = [c for c in bridge.calls if c.get("action") == "answer"]
    assert [(c["id"], c["text"]) for c in answers] == [("q1", "Toluol")]
    assert not [c for c in bridge.calls if c.get("action") == "prompt"]


@pytest.mark.parametrize("spoken", [" Ja. ", "yes", "okay", "Once."])
def test_a_spoken_yes_in_reply_approves_nothing_and_says_to_type_it(rig, spoken):
    rooms, _, bridge = rig
    mid = _request(rooms, CONFIRM)
    out = _voice(rooms, spoken, quote_id=mid)
    assert _decisions(bridge) == []
    assert out.startswith("🎙 ") and "tippen" in out.lower(), out


def test_a_spoken_no_in_reply_denies(rig):
    rooms, _, bridge = rig
    mid = _request(rooms, CONFIRM)
    _voice(rooms, "Nein!", quote_id=mid)
    assert [(c["id"], c["choice"]) for c in _decisions(bridge)] == [("c1", "n")]


def test_a_sentence_in_reply_to_an_approval_frees_nothing(rig):
    rooms, _, bridge = rig
    mid = _request(rooms, CONFIRM)
    for spoken in ("ja aber nicht jetzt", "Ja, mach das, und lösche danach alles.",
                   "", "hmm"):
        out = _voice(rooms, spoken, quote_id=mid)
        assert "nichts freigegeben" in out.lower(), (spoken, out)
    assert _decisions(bridge) == []


def test_a_lasting_rule_is_never_confirmed_by_voice(rig):
    rooms, _, bridge = rig
    mid = _request(rooms, {**CONFIRM, "confirmChoice": "A"})
    out = _voice(rooms, "Ja.", quote_id=mid)
    assert _decisions(bridge) == []
    assert "text" in out.lower()


def test_a_spoken_answer_to_a_question_is_the_whole_transcript(rig):
    rooms, _, bridge = rig
    mid = _request(rooms, {"id": "q1", "sessionId": "s1", "kind": "question",
                           "text": "Which solvent?"})
    _voice(rooms, "Wasser, bitte.", quote_id=mid)
    answers = [c for c in bridge.calls if c.get("action") == "answer"]
    assert [(c["id"], c["text"]) for c in answers] == [("q1", "Wasser, bitte.")]
