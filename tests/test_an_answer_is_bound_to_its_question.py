"""An answer on the phone is bound to its question; a gone question is never a prompt.

Field report (Pi, 2026-11-12): in a session that had asked nothing, the
model received "Antwort auf deine Rückfrage: Eins, zwei, drei." and
"Antwort auf deine Rückfrage: Erzähle mir einen Witz." — the person's
answers to question cards of an EARLIER session/chat.

Where the text is produced: delfin/dashboard/tab_agent.py
(answer_detected_question, wired as DashboardSession's ``question=``)
stuffs "Antwort auf deine Rückfrage: " + text into the input box and
sends it — whenever ``state['_pending_question']`` is set. That state is
NOT bound to the question the card was for:

* it is a question DETECTED in a turn's last text (detect_question), and
  it survives the turn, an answer given at the computer and a session
  restart in the same tab (the tab's ``state`` is one dict for every
  session it shows);
* ``DashboardSession.pending()`` publishes it as a pending question with
  no session binding, so the gateway sees a question waiting in the NEW
  session too, and a plain phone line is taken as its answer
  (answer_open_question answers "the one question" — chosen by count);
* the reply is then injected into the session that is live NOW — which
  had asked nothing.

Output:

* a question detected in another session's turn (or before a restart) is
  never published as this session's pending question, and an answer to
  its card is refused: nothing is injected into any session, and the
  phone is told the question is gone;
* after the question was answered at the computer (its state cleared),
  a late reply is refused with the same explanation;
* end to end through the gateway: a plain line in the chat of the NEW
  session reaches no session as a prompt while the only published
  question belongs to the old one.
"""
from __future__ import annotations

import asyncio

import pytest

from delfin.messenger import store
from delfin.messenger.dashboard import DashboardSession
from delfin.messenger.gateway import Rooms


DETECTED = {'type': 'open', 'options': []}


def adapter(sid, state, **callbacks):
    callbacks = {**dict(send=lambda text: True, stop=lambda: None,
                        answer=lambda text: False, plan=lambda yes: None), **callbacks}
    return DashboardSession(sid, state, **callbacks)


def test_an_answer_reaches_only_the_session_whose_turn_detected_the_question():
    """The detected question of an EARLIER session is not this session's
    question: its card's answer is refused, nothing is injected."""
    state = {'active_session_id': 'one', '_pending_question': dict(DETECTED),
             '_pending_question_sid': 'old', '_last_answer_text': 'First question?'}
    injected = []
    ui = adapter('one', state, question=lambda text: injected.append(text) or True)
    published = ui.pending()
    assert published == [], 'another session\'s detected question is not pending here'
    reply = ui({'method': 'permission.reply', 'requestId': 'whatever', 'text': 'Eins, zwei, drei.'})
    assert not reply['applied'] and injected == []
    assert 'no longer' in reply.get('error', '') or 'gone' in reply.get('error', '')


def test_a_question_answered_at_the_computer_is_gone_for_the_phone():
    """The person answered the detected question at the computer (the
    question UI hid it and the state was cleared): a late phone reply is
    refused and says so — it is not a prompt for anything."""
    state = {'active_session_id': 'one', '_pending_question': dict(DETECTED),
             '_pending_question_sid': 'one', '_last_answer_text': 'Which solvent?'}
    injected = []
    ui = adapter('one', state, question=lambda text: injected.append(text) or True)
    item = ui.pending()[0]
    state.pop('_pending_question')            # answered at the computer / turn ended
    state.pop('_pending_question_sid')
    ui.pending()                              # the refresh the next phone round does
    reply = ui({'method': 'permission.reply', 'requestId': item['requestId'],
                'text': 'Eins, zwei, drei.'})
    assert not reply['applied'] and injected == []
    assert 'no longer' in reply.get('error', '') or 'gone' in reply.get('error', '')


def test_the_live_question_still_takes_its_answer():
    """Control: the session's own, still-published detected question
    takes the phone's answer through the same door."""
    state = {'active_session_id': 'one', '_pending_question': dict(DETECTED),
             '_pending_question_sid': 'one', '_last_answer_text': 'Which solvent?'}
    injected = []
    ui = adapter('one', state, question=lambda text: injected.append(text) or True)
    item = ui.pending()[0]
    reply = ui({'method': 'permission.reply', 'requestId': item['requestId'], 'text': 'toluene'})
    assert reply['applied'] and injected == ['toluene']


# -- end to end: the gateway rig, the answer reaching the dashboard door ----

class Daemon:
    def __init__(self):
        self.sent = []

    async def send_chat(self, chat, text):
        self.sent.append((chat, text))
        return {"chatItems": [{"chatItem": {"meta": {"itemId": len(self.sent)}}}]}

    async def update_message(self, chat, mid, text):
        self.sent.append(('update', chat, mid, text))


class Bridge:
    """The bridge's session door the way bridge._cmd_session answers a
    question: a dashboard-backed session goes to its adapter (the
    bridge's kind translation ask->question and its ok/error mapping
    included), and a plain line while a question waits is taken as its
    answer."""

    def __init__(self, adapter_of):
        self.adapter_of = adapter_of
        self.prompts = []

    def ask(self, command):
        if command['cmd'] == 'sessions':
            return {'sessions': [{'id': sid, 'name': 'same title', 'state': 'idle',
                                  'messengerEnabled': True} for sid in self.adapter_of]}
        sid = command['session']
        ui = self.adapter_of[sid]
        action = command['action']
        if action == 'status':
            return {'ok': True, 'title': 'same title', 'state': 'idle',
                    'messengerEnabled': True}
        if action == 'pending':
            items = [{'kind': 'question' if i['kind'] == 'ask' else i['kind'],
                      'id': i['requestId'], 'sessionId': sid, **{
                          k: v for k, v in i.items()
                          if k not in {'kind', 'requestId'}}} for i in ui.pending()]
            return {'ok': True, 'items': items}
        if action == 'answer':
            reply = ui({'method': 'permission.reply', 'requestId': command['id'],
                        'text': command['text']})
            return {'ok': bool(reply.get('applied')), 'error': reply.get('error')}
        if action == 'prompt':
            self.prompts.append((sid, command['text']))
            return {'ok': True, 'accepted': True}
        return {'ok': True}


@pytest.fixture
def rig(tmp_path, monkeypatch):
    monkeypatch.setenv('DELFIN_MESSENGER_DIR', str(tmp_path))
    store.save('7')
    store.put_room('#10', 's-new', 'same title', owner='7', ready=True)
    return store


def _message(chat='#10', text='calculate', sender='7', mid='1'):
    return {'chat': chat, 'text': text, 'sender': sender, 'id': mid}


def test_the_answer_to_an_old_card_is_never_a_prompt_of_the_new_session(rig):
    """The card was for the EARLIER session's detected question; the chat
    is bound to the NEW session now. The plain line answers nothing and
    reaches no session as a prompt."""
    state = {'active_session_id': 's-new', '_pending_question': dict(DETECTED),
             '_pending_question_sid': 's-old', '_last_answer_text': 'First question?'}
    injected = []
    ui = adapter('s-new', state, question=lambda text: injected.append(text) or True)
    bridge = Bridge({'s-new': ui})
    rooms = Rooms(Daemon(), bridge)
    asyncio.run(rooms.handle(_message(text='Eins, zwei, drei.', mid='a')))
    assert injected == [], 'the answer went to the new session as a prompt'
    assert bridge.prompts == [], 'the answer reached a session as a prompt'
    reply = rooms.daemon.sent[-1][1] if rooms.daemon.sent else ''
    assert 'nicht mehr zugeordnet' in reply or 'no longer' in reply


def test_a_line_of_the_new_session_itself_still_reaches_it(rig):
    """Control: with no stale question published, the same line is the
    new session's prompt as always."""
    state = {'active_session_id': 's-new'}
    ui = adapter('s-new', state, question=lambda text: True)
    bridge = Bridge({'s-new': ui})
    rooms = Rooms(Daemon(), bridge)
    asyncio.run(rooms.handle(_message(text='Erzähle mir einen Witz.', mid='b')))
    assert bridge.prompts == [('s-new', 'Erzähle mir einen Witz.')]
