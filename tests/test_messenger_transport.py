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
                daemon.command({'type': 'one'}), daemon.command({'type': 'two'})), 2)
            assert replies == [{'type': 'ok'}, {'type': 'ok'}]
            events = daemon.events()
            assert (await anext(events))['type'] == 'newChatItems'
            assert (await anext(events))['type'] == 'newChatItems'
            await events.aclose()
        finally:
            reader.cancel()
            await asyncio.gather(reader, return_exceptions=True)
    asyncio.run(run())
