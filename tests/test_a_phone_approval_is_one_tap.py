"""An approval card on the phone carries two tappable commands.

SimpleX (6.4.3+) highlights /'...' in a bot's message and sends it on a tap.
A tap is a new message, not a reply to the card, so the card's random token
is what binds it: /ja TOKEN answers exactly that card in that chat, through
the same checks as a reply, once, and not after the reply window.
"""
from __future__ import annotations

import asyncio
import re

import pytest

from test_the_messenger_gateway_only_serves_the_paired_phone import message
from test_the_messenger_gateway_only_serves_the_paired_phone import rig as _rig_fixture

rig = pytest.fixture(name='rig')(_rig_fixture.__wrapped__) if hasattr(_rig_fixture, '__wrapped__') else _rig_fixture

CARD = {'id': 'abc', 'sessionId': 's1', 'kind': 'confirm', 'text': 'run calculation',
        'phoneApproval': True, 'options': [{'key': 'y', 'label': 'yes'}, {'key': 'n', 'label': 'no'}]}


def _token(daemon, chat='#10'):
    text = next(t for c, t in daemon.sent if c == chat and 'run calculation' in t)
    return re.search(r"/'(?:ja|yes) ([0-9a-f]{4})'", text).group(1)


def test_a_tap_answers_the_card_like_a_reply(rig):
    rooms, daemon, bridge = rig
    bridge.pending['s1'] = [dict(CARD)]
    asyncio.run(rooms.refresh())
    token = _token(daemon)
    asyncio.run(rooms.handle(message(text=f'/ja {token}', mid='t1')))
    tapped = list(bridge.calls)
    assert tapped and all(c['session'] == 's1' for c in tapped)
    assert any(c.get('id') == 'abc' for c in tapped), tapped


def test_a_tap_works_once_and_only_in_its_chat(rig):
    rooms, daemon, bridge = rig
    bridge.pending['s1'] = [dict(CARD)]
    asyncio.run(rooms.refresh())
    token = _token(daemon)
    asyncio.run(rooms.handle(message('#20', f'/ja {token}', mid='x1')))
    assert not any(c.get('id') == 'abc' for c in bridge.calls)
    asyncio.run(rooms.handle(message(text=f'/ja {token}', mid='t1')))
    n = len(bridge.calls)
    asyncio.run(rooms.handle(message(text=f'/ja {token}', mid='t2')))
    assert not any(c.get('id') == 'abc' for c in bridge.calls[n:])


def test_an_expired_tap_approves_nothing(rig):
    rooms, daemon, bridge = rig
    bridge.pending['s1'] = [dict(CARD)]
    asyncio.run(rooms.refresh())
    token = _token(daemon)
    mid, _when = rooms._taps[('#10', token)]
    rooms._taps[('#10', token)] = (mid, -1e9)
    asyncio.run(rooms.handle(message(text=f'/ja {token}', mid='t1')))
    assert not any(c.get('id') == 'abc' for c in bridge.calls)


def test_a_card_the_phone_may_not_approve_has_no_tap(rig):
    rooms, daemon, bridge = rig
    bridge.pending['s1'] = [{**CARD, 'phoneApproval': False}]
    asyncio.run(rooms.refresh())
    text = next(t for c, t in daemon.sent if c == '#10' and 'run calculation' in t)
    assert "/'ja" not in text and "/'yes" not in text
