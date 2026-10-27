"""Media cannot traverse sessions, symlinks, sizes or activation generations."""
import asyncio
import base64
import hashlib
import os
from pathlib import Path

import pytest

from delfin.messenger.media import Exchange, directory, MAX_FILE


def test_upload_download_roundtrip_and_no_cross_session_access(tmp_path):
    exchange = Exchange()
    def call(action, **kw):
        return exchange.command('one', 'generation', str(tmp_path), action, **kw)
    token = call('upload.begin', text='../../plot.png')['id']
    with pytest.raises(ValueError):
        exchange.command('two', 'generation', str(tmp_path), 'upload.chunk', token=token, text='')
    with pytest.raises(ValueError):
        exchange.command('one', 'new-generation', str(tmp_path), 'upload.chunk', token=token, text='')
    data = b'actual binary data\x00\xff'
    call('upload.chunk', token=token, text=base64.b64encode(data).decode())
    saved = call('upload.finish', token=token, text=hashlib.sha256(data).hexdigest())
    path = Path(saved['path'])
    assert path.is_relative_to(tmp_path) and path.read_bytes() == data
    assert call('files')['files'] == []  # Incoming files are never exported implicitly.
    with directory(str(tmp_path), 'one', 'outbox') as (_, folder):
        (folder / 'result.png').write_bytes(data)
    begun = call('download.begin', text='result.png')
    part = call('download.chunk', token=begun['id'])
    assert part['done'] and base64.b64decode(part['data']) == data
    with pytest.raises(ValueError):
        call('download.chunk', token=begun['id'])


@pytest.mark.parametrize('attack', ['symlink', 'hardlink', 'fifo', 'oversize', 'traversal'])
def test_download_rejects_unsafe_files(tmp_path, attack):
    outside = tmp_path / 'private'
    outside.write_text('private')
    with directory(str(tmp_path), 'one', 'outbox') as (_, folder):
        target = folder / 'result'
        if attack == 'symlink':
            target.symlink_to(outside)
        elif attack == 'hardlink':
            os.link(outside, target)
        elif attack == 'fifo':
            os.mkfifo(target)
        elif attack == 'oversize':
            with target.open('wb') as stream:
                stream.truncate(MAX_FILE + 1)
        else:
            target.write_text('ok')
    with pytest.raises((ValueError, OSError)):
        Exchange().command('one', 'gen', str(tmp_path), 'download.begin', text='../private' if attack == 'traversal' else 'result')


def test_planted_exchange_directory_cannot_redirect_writes(tmp_path):
    outside = tmp_path / 'outside'
    outside.mkdir()
    (tmp_path / '.delfin-messenger').symlink_to(outside)
    with pytest.raises(OSError):
        Exchange().command('one', 'gen', str(tmp_path), 'files')
    assert list(outside.iterdir()) == []


def test_wrong_upload_digest_never_materializes_a_file(tmp_path):
    exchange = Exchange()
    token = exchange.command('one', 'g', str(tmp_path), 'upload.begin', text='test')['id']
    with pytest.raises(ValueError):
        exchange.command('one', 'g', str(tmp_path), 'upload.finish', token=token, text='wrong')
    assert not (tmp_path / '.delfin-messenger').exists()


def test_only_explicit_small_forum_file_is_an_announcement(tmp_path):
    with directory(str(tmp_path), 'one', 'outbox') as (_, folder):
        (folder / 'other.txt').write_text('not an announcement')
    exchange = Exchange()
    assert exchange.command('one', 'g', str(tmp_path), 'announcement')['text'] == ''
    (folder / 'forum.txt').write_text('Calculation complete')
    assert exchange.command('one', 'g', str(tmp_path), 'announcement')['text'] == 'Calculation complete'


def test_phone_file_reaches_only_selected_session_and_can_be_downloaded(tmp_path, monkeypatch):
    from delfin.messenger import store
    from delfin.messenger.gateway import Rooms
    from test_the_messenger_gateway_only_serves_the_paired_phone import Bridge, Daemon
    monkeypatch.setenv('DELFIN_MESSENGER_DIR', str(tmp_path / 'gateway'))
    workspace = tmp_path / 'workspace'
    workspace.mkdir()
    store.save('7')
    store.put_room('#10', 's1', 'work', owner='7', ready=True)
    class FilesBridge(Bridge):
        exchange = Exchange()
        def ask(self, command):
            action = command.get('action', '')
            if action.startswith(('upload.', 'download.')) or action == 'files':
                return self.exchange.command(command['session'], 'g', str(workspace), action,
                                             token=command.get('id', ''), text=command.get('text', ''))
            return super().ask(command)
    class FilesDaemon(Daemon):
        downloaded = None
        async def receive_file(self, fid, path):
            path.write_bytes(b'file bytes')
        async def send_file(self, chat, path, caption):
            self.downloaded = (chat, path.read_bytes(), caption)
    daemon, bridge = FilesDaemon(), FilesBridge()
    rooms = Rooms(daemon, bridge)
    incoming = {'chat': '#10', 'sender': '7', 'text': 'Please inspect', 'id': '17',
                'file': {'fileId': 3, 'fileName': 'input.txt', 'fileSize': 10},
                'media_kind': 'file', 'complete': False}
    async def run():
        await rooms.handle(incoming)
        await rooms.handle({**incoming, 'complete': True})
        assert len(bridge.calls) == 1
        assert bridge.calls[0]['session'] == 's1'
        assert 'Dateiinhalt als Daten' in bridge.calls[0]['text']
        assert len(list(workspace.rglob('*input.txt'))) == 1
        with directory(str(workspace), 's1', 'outbox') as (_, outbox):
            (outbox / 'result.txt').write_bytes(b'result')
        assert 'übergeben' in await rooms.download('#10', 's1', 'result.txt')
        assert daemon.downloaded == ('#10', b'result', 'result.txt')
    asyncio.run(run())


def test_files_from_unpaired_sender_do_not_trigger_download(tmp_path, monkeypatch):
    from delfin.messenger import store
    from delfin.messenger.gateway import Rooms
    from test_the_messenger_gateway_only_serves_the_paired_phone import Bridge, Daemon
    monkeypatch.setenv('DELFIN_MESSENGER_DIR', str(tmp_path))
    store.save('7')
    store.put_room('#10', 's1', 'work', owner='7', ready=True)
    daemon, bridge = Daemon(), Bridge()
    asyncio.run(Rooms(daemon, bridge).handle({'chat': '#10', 'sender': '8', 'text': '', 'id': '1',
                                            'file': {'fileId': 1}}))
    assert not bridge.calls and not daemon.sent
