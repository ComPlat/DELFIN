"""Renaming from the phone works: typed, tapped and after a restart.

Seen live (2026-11-12): ``name Test`` in a session's SimpleX group did not
rename the session. The phone's menu shows ``name`` as a tappable command
(SimpleX 6.4.3+ highlights ``/'...'`` and sends it on a tap), so the line
can arrive either as the person typed it or as the tapped form.
"""

from __future__ import annotations

import asyncio

import pytest

from delfin.messenger import store
from delfin.messenger.gateway import Rooms, group_name


class Daemon:
    def __init__(self):
        self.renamed = []

    async def create_group(self, title, commands=None):
        return "#501"

    async def rename_group(self, chat, title, commands=None):
        self.renamed.append((chat, title))

    async def invite(self, chat, cid):
        pass

    async def send_chat(self, chat, text, **kw):
        pass


class Bridge:
    """The service's answers; the session's title lives here."""

    def __init__(self, titles=None):
        self.calls = []
        self.titles = dict(titles or {})

    def ask(self, command):
        self.calls.append(command)
        if command['cmd'] == 'sessions':
            return {'sessions': [{'id': sid, 'name': self.titles.get(sid, ''),
                                  'state': 'idle', 'messengerEnabled': True}
                                 for sid in ['s1']]}
        sid = command['session']
        if command['action'] == 'status':
            return {'ok': True, 'title': self.titles.get(sid, ''),
                    'state': 'idle', 'messengerEnabled': True}
        if command['action'] == 'rename':
            title = str(command.get('text') or '')
            self.titles[sid] = title
            return {'ok': True, 'title': title}
        return {'ok': True}


def message(chat='#10', text='calculate', sender='7', mid='1'):
    return {'chat': chat, 'text': text, 'sender': sender, 'id': mid}


@pytest.fixture
def rig(tmp_path, monkeypatch):
    monkeypatch.setenv('DELFIN_MESSENGER_DIR', str(tmp_path))
    store.save('7')
    store.put_room('#10', 's1', '/sync', owner='7', ready=True)
    daemon, bridge = Daemon(), Bridge(titles={'s1': '/sync'})
    return Rooms(daemon, bridge), daemon, bridge


def _renamed(rig):
    return rig[1].renamed


def _rename_calls(rig):
    return [c for c in rig[2].calls if c.get('action') == 'rename']


def test_a_typed_name_renames_the_session_and_the_group(rig):
    rooms, daemon, _ = rig
    result = asyncio.run(rooms.handle(message(text='name Test', mid='2')))
    assert _rename_calls(rig) == [
        {'cmd': 'session', 'session': 's1', 'action': 'rename', 'text': 'Test'}]
    assert _renamed(rig) == [('#10', group_name('Test', 's1'))]
    assert store.rooms()['#10']['title'] == 'Test'


@pytest.mark.parametrize("tapped", ["/'name Test'", "/name Test"])
def test_a_tapped_name_is_a_name_not_a_task(rig, tapped):
    # The menu's tappable form must do what the typed word does, not become
    # a message to the agent (the old behaviour made "name Test" a prompt).
    rooms, daemon, _ = rig
    result = asyncio.run(rooms.handle(message(text=tapped, mid='3')))
    assert _rename_calls(rig), result
    assert _renamed(rig) == [('#10', group_name('Test', 's1'))]
    assert store.rooms()['#10']['title'] == 'Test'


def test_a_session_titled_a_slash_command_renames_from_the_phone(rig):
    # The title the first line left ('/sync'): a name command still lands,
    # and the new name is what the chat carries afterwards.
    rooms, daemon, _ = rig
    result = asyncio.run(rooms.handle(message(text='name Butan', mid='4')))
    assert store.rooms()['#10']['title'] == 'Butan'
    assert _renamed(rig) == [('#10', group_name('Butan', 's1'))]
