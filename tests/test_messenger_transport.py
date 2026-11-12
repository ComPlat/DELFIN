"""Transport framing and concurrent notification/command delivery."""
import asyncio
import json
import socket
import threading

from delfin.messenger.bridge_client import BridgeClient
from delfin.messenger.simplex import Daemon


def test_bridge_reply_can_share_a_socket_read_with_an_event(tmp_path):
    path = tmp_path / 'bridge.sock'
    server = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    server.bind(str(path))
    server.listen()
    def respond():
        conn, _ = server.accept()
        with conn:
            conn.recv(8192)
            conn.sendall(b'{"event":"session"}\n{"ok":true}\n')
    thread = threading.Thread(target=respond)
    thread.start()
    try:
        assert BridgeClient(path).ask({'cmd': 'status'}) == {'ok': True}
    finally:
        thread.join(2)
        server.close()


def test_simplex_preserves_events_while_commands_wait_for_replies():
    class Socket:
        def __init__(self):
            self.incoming = asyncio.Queue()
        def __aiter__(self):
            return self
        async def __anext__(self):
            return await self.incoming.get()
        async def send(self, raw):
            message = json.loads(raw)
            await self.incoming.put(json.dumps({'resp': {'type': 'newChatItems'}}))
            await self.incoming.put(json.dumps({'corrId': message['corrId'],
                                               'resp': {'type': 'ok'}}))
    async def run():
        daemon = Daemon()
        daemon._ws = Socket()
        reader = asyncio.create_task(daemon._read())
        try:
            replies = await asyncio.wait_for(asyncio.gather(
                daemon.command('/user'), daemon.command('/_contacts 1')), 2)
            assert replies == [{'type': 'ok'}, {'type': 'ok'}]
            events = daemon.events()
            assert (await anext(events))['type'] == 'newChatItems'
            assert (await anext(events))['type'] == 'newChatItems'
            await events.aclose()
        finally:
            reader.cancel()
            await asyncio.gather(reader, return_exceptions=True)
    asyncio.run(run())


def received(group=10, sender=7, direction='groupRcv', mid=44):
    return {'type': 'newChatItems', 'user': {'userId': 1}, 'chatItems': [{
        'chatInfo': {'type': 'group', 'groupInfo': {'groupId': group}},
        'chatItem': {'chatDir': {'type': direction, 'groupMember': {'memberContactId': sender}},
                     'meta': {'itemId': mid},
                     'content': {'type': 'rcvMsgContent', 'msgContent': {'type': 'text', 'text': 'hello'}}}}]}


def test_real_simplex_envelope_keeps_all_received_items_and_ignores_own_messages():
    from delfin.messenger.simplex import messages_of
    event = received()
    event['chatItems'] += received(group=20, mid=45)['chatItems']
    assert messages_of(event) == [
        {'chat': '#10', 'sender': '7', 'text': 'hello', 'id': '44'},
        {'chat': '#20', 'sender': '7', 'text': 'hello', 'id': '45'}]
    assert messages_of(received(direction='groupSnd')) == []
    assert messages_of(received(sender=None)) == []


def test_simplex_commands_use_documented_string_syntax():
    async def run():
        daemon = Daemon()
        commands = []
        async def command(text):
            commands.append(text)
            if text == '/user':
                return {'type': 'activeUser', 'user': {'userId': 1}}
            if text.startswith('/_group'):
                return {'type': 'groupCreated', 'groupInfo': {'groupId': 10}}
            if text.startswith('/_connect'):
                return {'type': 'invitation', 'connLinkInvitation': {'connFullLink': 'simplex:/inv'},
                        'connection': {'pccConnId': 9}}
            return {'type': 'cmdOk'}
        daemon.command = command
        assert await daemon.create_group('DELFIN session') == '#10'
        await daemon.invite('#10', '7')
        await daemon.send_chat('#10', 'hello')
        assert await daemon.pending_invitation() == 'simplex:/inv'
        assert daemon.invitation_id == 9
        assert commands[1].startswith('/_group 1 {')
        assert '/_add #10 7 member' in commands
        sent = next(s for s in commands if s.startswith('/_send'))
        assert json.loads(sent.split(' json ', 1)[1])[0]['msgContent']['text'] == 'hello'
    asyncio.run(run())


def test_simplex_reply_retains_original_message_id():
    from delfin.messenger.simplex import messages_of
    event = received()
    event['chatItems'][0]['chatItem']['quotedItem'] = {'itemId': 123}
    assert messages_of(event)[0]['quote_id'] == '123'


def test_simplex_media_events_preserve_sender_and_completion():
    from delfin.messenger.simplex import messages_of
    event = received()
    wrapped = event['chatItems'][0]
    wrapped['chatItem']['content']['msgContent'] = {'type': 'image', 'text': 'look'}
    wrapped['chatItem']['file'] = {'fileId': 1, 'fileSize': 100, 'fileName': 'a.png',
                                  'fileStatus': {'type': 'rcvInvitation'}}
    message = messages_of(event)[0]
    assert message['sender'] == '7' and message['media_kind'] == 'image'
    assert not message['complete']
    complete = messages_of({'type': 'rcvFileComplete', 'chatItem': wrapped})[0]
    assert complete['complete'] and complete['file']['fileId'] == 1


def test_file_and_group_rename_use_pinned_simplex_api():
    async def run():
        daemon = Daemon()
        commands = []
        async def command(text):
            commands.append(text)
            return {}
        daemon.command = command
        await daemon.rename_group('#10', 'Research')
        await daemon.receive_file(42, '/tmp/file.png')
        await daemon.send_file('#10', '/tmp/result.pdf', 'Report')
        assert commands[0].startswith('/_group_profile #10 ')
        assert commands[1] == '/freceive 42 encrypt=off /tmp/file.png'
        body = json.loads(commands[2].split(' json ')[1])
        assert body[0]['fileSource']['filePath'] == '/tmp/result.pdf'
        assert body[0]['msgContent']['type'] == 'file'
    asyncio.run(run())


def test_xftp_is_accepted_only_after_its_description_arrives():
    from delfin.messenger.simplex import messages_of
    event = received()
    wrapped = event['chatItems'][0]
    wrapped['chatItem']['content']['msgContent'] = {'type': 'file', 'text': ''}
    wrapped['chatItem']['file'] = {'fileId': 1, 'fileSize': 100, 'fileName': 'a.pdf',
                                  'fileProtocol': 'xftp', 'fileStatus': {'type': 'rcvInvitation'}}
    assert messages_of(event)[0]['await_description']
    ready = messages_of({'type': 'rcvFileDescrReady', 'chatItem': wrapped})[0]
    assert not ready['await_description'] and not ready['complete']
