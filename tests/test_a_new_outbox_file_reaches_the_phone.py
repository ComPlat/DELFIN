"""A file the session puts in its outbox reaches the phone without typing.

Live report (2026-11-12): the person wrote "schick mir ein Bild"; the
agent rendered benzene.png into the outbox and told them to fetch it
with "datei benzene.png" -- a whole extra step of typing for something
the phone could just show. Since then the gateway offers what appears
in a session's outbox on its own:

* an image under 1 MB is pushed directly (send_file) -- the phone shows
  it, no text needed;
* every other file sends one message with its name and a tappable
  SimpleX command link /'datei NAME' -- one tap, any size;
* only files that appear after the gateway has seen the chat once are
  offered (at start the existing outbox is the baseline, no replay);
* each file is offered once; the outbox itself stays the only source --
  nothing is offered that is not in it.
"""
import asyncio
import base64
import hashlib
import os

import pytest

from delfin.messenger import store
from delfin.messenger.gateway import Rooms

SID = 's1'
OWN = hashlib.sha256(SID.encode()).hexdigest()[:24]


class Daemon:
    def __init__(self):
        self.sent = []
        self.files = []
    async def send_chat(self, chat, text):
        self.sent.append((chat, text))
        return {"chatItems": [{"chatItem": {"meta": {"itemId": len(self.sent)}}}]}
    async def send_file(self, chat, path, name):
        self.files.append((chat, name, os.path.basename(path)))
    async def create_group(self, title):
        return '#' + title
    async def invite(self, chat, owner):
        pass
    async def rename_group(self, chat, title):
        return {'chatItems': []}


class Bridge:
    def __init__(self, outbox):
        self.outbox = outbox
        self.active = {SID}
        self.chunk_name = ''
    def ask(self, command):
        if command['cmd'] == 'sessions':
            return {'sessions': [{'id': SID, 'name': 'title', 'state': 'idle',
                                  'messengerEnabled': True}]}
        sid, action = command['session'], command.get('action')
        if action == 'status':
            return {'ok': True, 'state': 'idle', 'answer': '', 'messengerEnabled': True}
        if action == 'pending':
            return {'ok': True, 'items': []}
        if action == 'files':
            files = [{'name': n, 'size': os.path.getsize(os.path.join(self.outbox, n))}
                     for n in sorted(os.listdir(self.outbox))]
            return {'ok': True, 'files': files, 'outbox': self.outbox}
        if action == 'download.begin':
            name = command.get('text', '')
            data = open(os.path.join(self.outbox, name), 'rb').read()
            self.chunk_name = name
            return {'ok': True, 'id': 'x', 'name': name, 'sha256': hashlib.sha256(data).hexdigest()}
        if action == 'download.chunk':
            return {'ok': True, 'data': base64.b64encode(
                open(os.path.join(self.outbox, self.chunk_name), 'rb').read()).decode(),
                'done': True}
        return {'ok': True}


@pytest.fixture
def rig(tmp_path, monkeypatch):
    monkeypatch.setenv('DELFIN_MESSENGER_DIR', str(tmp_path))
    store.save('7')
    store.put_room('#10', SID, 'title', owner='7', ready=True)
    outbox = tmp_path / '.delfin-messenger' / OWN / 'outbox'
    outbox.mkdir(parents=True)
    daemon, bridge = Daemon(), Bridge(str(outbox))
    rooms = Rooms(daemon, bridge)
    rooms._active = {SID}
    return rooms, daemon, bridge, outbox


def test_an_image_under_a_megabyte_is_pushed_to_the_phone(rig):
    rooms, daemon, _, outbox = rig
    asyncio.run(rooms.refresh())                       # baseline
    (outbox / 'benzene.png').write_bytes(b'\x89PNG' + b'x' * 100)
    asyncio.run(rooms.refresh())
    assert [f[:2] for f in daemon.files] == [('#10', 'benzene.png')]
    assert not daemon.sent                             # the image is the message


def test_a_document_sends_one_message_with_a_tappable_fetch(rig):
    rooms, daemon, _, outbox = rig
    asyncio.run(rooms.refresh())
    (outbox / 'report.pdf').write_bytes(b'x' * 10)
    asyncio.run(rooms.refresh())
    assert not daemon.files
    assert len(daemon.sent) == 1
    text = daemon.sent[0][1]
    assert 'report.pdf' in text
    assert "/'datei report.pdf'" in text


def test_a_large_image_is_offered_not_pushed(rig):
    rooms, daemon, _, outbox = rig
    asyncio.run(rooms.refresh())
    (outbox / 'big.png').write_bytes(b'x' * (1024 * 1024 + 1))
    asyncio.run(rooms.refresh())
    assert not daemon.files
    assert daemon.sent and 'big.png' in daemon.sent[0][1]


def test_the_baseline_at_start_is_not_replayed_and_nothing_is_offered_twice(rig):
    rooms, daemon, _, outbox = rig
    (outbox / 'old.txt').write_text('x')
    asyncio.run(rooms.refresh())
    assert not daemon.sent and not daemon.files        # pre-existing: baseline
    (outbox / 'old.txt').write_text('x')               # unchanged: still once
    asyncio.run(rooms.refresh())
    (outbox / 'new.txt').write_text('y')
    asyncio.run(rooms.refresh())
    assert not daemon.files and 'new.txt' in daemon.sent[0][1]
    asyncio.run(rooms.refresh())
    assert len(daemon.sent) == 1 and not daemon.files


def test_only_the_session_s_own_outbox_is_the_source(rig):
    rooms, daemon, bridge, outbox = rig
    asyncio.run(rooms.refresh())
    # A path the session did not report is never the source: the bridge's
    # 'outbox' is where files come from; anything else must not be read.
    assert bridge.outbox.endswith(os.path.join('.delfin-messenger', OWN, 'outbox'))
    (outbox / 'img.jpg').write_bytes(b'x' * 5)
    asyncio.run(rooms.refresh())
    assert [(c, n) for c, n, _ in daemon.files] == [('#10', 'img.jpg')]
