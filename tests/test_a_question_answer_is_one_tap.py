"""A question card's options are tappable: /'N TOKEN' answers that option.

Like the approval taps (/ja TOKEN, /nein TOKEN), a question tap is a new
message -- the card's random token is what binds it to this one card in
this one chat, through the same checks as a reply to it, once, and not
after the reply window. A question with multiSelect keeps typed answers
only (no number taps). A cut answer's ``more`` hint is a tappable
/'more' (English) or /'mehr' (German).
"""
from __future__ import annotations

import asyncio
import re

import pytest

from delfin.messenger import render
from test_the_messenger_gateway_only_serves_the_paired_phone import message
from test_the_messenger_gateway_only_serves_the_paired_phone import rig as _rig_fixture

rig = pytest.fixture(name='rig')(_rig_fixture.__wrapped__) if hasattr(_rig_fixture, '__wrapped__') else _rig_fixture

QCARD = {'id': 'q1', 'sessionId': 's1', 'kind': 'question',
         'question': 'What color?', 'choices': ['red', 'blue']}
MULTI = {'id': 'q2', 'sessionId': 's1', 'kind': 'question',
         'question': 'Which colors?', 'choices': ['red', 'blue'], 'multiSelect': True}


def _token(daemon, chat='#10'):
    text = next(t for c, t in daemon.sent if c == chat and 'What color' in t)
    return re.search(r"/'\d ([0-9a-f]{4})'", text).group(1)


def test_a_tap_answers_the_question_like_a_reply(rig):
    rooms, daemon, bridge = rig
    bridge.pending['s1'] = [dict(QCARD)]
    asyncio.run(rooms.refresh())
    token = _token(daemon)
    asyncio.run(rooms.handle(message(text=f'/2 {token}', mid='t1')))
    assert any(c.get('action') == 'answer' and c.get('text') == '2' for c in bridge.calls), bridge.calls


def test_a_tap_works_once_and_only_in_its_chat(rig):
    rooms, daemon, bridge = rig
    bridge.pending['s1'] = [dict(QCARD)]
    asyncio.run(rooms.refresh())
    token = _token(daemon)
    asyncio.run(rooms.handle(message('#20', f'/1 {token}', mid='x1')))
    assert not any(c.get('id') == 'q1' for c in bridge.calls)
    asyncio.run(rooms.handle(message(text=f'/1 {token}', mid='t1')))
    n = len(bridge.calls)
    asyncio.run(rooms.handle(message(text=f'/1 {token}', mid='t2')))
    assert not any(c.get('id') == 'q1' for c in bridge.calls[n:])


def test_an_expired_tap_answers_nothing(rig):
    rooms, daemon, bridge = rig
    bridge.pending['s1'] = [dict(QCARD)]
    asyncio.run(rooms.refresh())
    token = _token(daemon)
    mid, _when = rooms._taps[('#10', token)]
    rooms._taps[('#10', token)] = (mid, -1e9)
    asyncio.run(rooms.handle(message(text=f'/2 {token}', mid='t1')))
    assert not any(c.get('id') == 'q1' for c in bridge.calls)


def test_a_multiselect_question_has_no_tap(rig):
    rooms, daemon, bridge = rig
    bridge.pending['s1'] = [dict(MULTI)]
    asyncio.run(rooms.refresh())
    text = next(t for c, t in daemon.sent if c == '#10' and 'Which colors' in t)
    assert "/'" not in text


def test_a_question_with_no_choices_has_no_tap(rig):
    rooms, daemon, bridge = rig
    bridge.pending['s1'] = [dict(QCARD, choices=[])]
    asyncio.run(rooms.refresh())
    text = next(t for c, t in daemon.sent if c == '#10' and 'What color' in t)
    assert "/'" not in text


def test_more_hint_is_tappable_in_english(monkeypatch):
    monkeypatch.setenv('DELFIN_MESSENGER_LANG', 'en')
    assert "/'more'" in render.more_hint()


def test_more_hint_is_tappable_in_german(monkeypatch):
    monkeypatch.setenv('DELFIN_MESSENGER_LANG', 'de')
    assert "/'mehr'" in render.more_hint()
