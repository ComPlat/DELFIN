"""A voice message does what the same words typed would do, and no more.

Input: Rooms.voice (delfin/messenger/gateway.py), called with the local
transcript of a voice message, its chat and, when the person used the reply
gesture, the message it answers.

- Without a reply reference the transcript goes to the session as its next
  prompt, at once (before: it was shown and needed a typed "senden").
- As a reply to an approval, a transcript that is exactly a decision word
  ("Ja.", "nein!", "yes", an offered option's label) answers it like the
  typed word, and the phone is told what was understood.
- Anything else in reply to an approval -- a sentence, a mumble, "ja aber
  nicht jetzt" -- approves nothing and says so.
- A lasting rule is never confirmed by voice; that stays typed.
- As a reply to a question the whole transcript is the answer.
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
                             "text": "Sprachnachricht des Nutzers (Transkript):\n"
                                     "Optimiere bitte beide Konformere."}]
    assert "Optimiere bitte beide Konformere." in out


def test_a_spoken_yes_in_reply_answers_the_approval(rig):
    rooms, _, bridge = rig
    mid = _request(rooms, CONFIRM)
    out = _voice(rooms, " Ja. ", quote_id=mid)
    assert [(c["id"], c["choice"]) for c in _decisions(bridge)] == [("c1", "y")]
    assert "ja" in out.lower()


def test_a_spoken_no_in_reply_denies(rig):
    rooms, _, bridge = rig
    mid = _request(rooms, CONFIRM)
    _voice(rooms, "Nein!", quote_id=mid)
    assert [(c["id"], c["choice"]) for c in _decisions(bridge)] == [("c1", "n")]


def test_an_offered_label_spoken_in_reply_is_that_option(rig):
    rooms, _, bridge = rig
    mid = _request(rooms, CONFIRM)
    _voice(rooms, "Once.", quote_id=mid)
    assert [(c["id"], c["choice"]) for c in _decisions(bridge)] == [("c1", "y")]


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
