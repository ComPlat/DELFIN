"""A session's SimpleX chat takes the session's name.

Seen on the phone on 2026-11-12: a chat made while its session was untitled
was named after the head of the session id ("only a number") and kept it;
the session's later title never reached the chat.
"""
from __future__ import annotations

import asyncio

from delfin.messenger import gateway, store
from test_the_messenger_gateway_only_serves_the_paired_phone import message  # noqa: F401
from test_the_messenger_gateway_only_serves_the_paired_phone import rig  # noqa: F401


def _renames(rooms):
    seen = []

    async def rename_group(chat, title, commands=None):
        seen.append((chat, title))
    rooms.daemon.rename_group = rename_group
    return seen


def test_an_untitled_session_chat_is_not_a_number():
    name = gateway.group_name('', 'a1b2c3d4e5f6', '')
    assert 'a1b2c3d4' not in name


def test_the_chat_follows_a_later_title(rig):  # noqa: F811
    rooms, daemon, bridge = rig
    seen = _renames(rooms)
    asyncio.run(rooms._follow_title('#10', store.rooms()['#10'], {'title': 'Spin states'}))
    assert seen and seen[-1][0] == '#10' and seen[-1][1].endswith('Spin states')
    assert store.rooms()['#10']['title'] == 'Spin states'
    seen.clear()
    asyncio.run(rooms._follow_title('#10', store.rooms()['#10'], {'title': 'Spin states'}))
    assert not seen                                 # unchanged: nothing sent


def test_a_command_title_or_an_error_is_not_taken(rig):  # noqa: F811
    rooms, daemon, bridge = rig
    seen = _renames(rooms)
    asyncio.run(rooms._follow_title('#10', store.rooms()['#10'], {'title': '/sync'}))
    asyncio.run(rooms._follow_title('#10', store.rooms()['#10'], {'error': 'x', 'title': 'B'}))
    asyncio.run(rooms._follow_title('#10', store.rooms()['#10'], {'title': ''}))
    assert not seen
