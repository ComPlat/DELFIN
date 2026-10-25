"""Chat IDs, never titles or the last active agent, determine recipients."""
import asyncio

import pytest

from delfin.messenger import store
from delfin.messenger.gateway import Rooms, UPDATE


class Daemon:
    def __init__(self):
        self.sent = []
        self.groups = []
        self.invited = []
    async def send_chat(self, chat, text):
        self.sent.append((chat, text))
    async def create_group(self, title):
        self.groups.append(title)
        return '#' + str(100 + len(self.groups))
    async def invite(self, chat, owner):
        self.invited.append((chat, owner))


class Bridge:
    def __init__(self):
        self.calls = []
        self.missing = set()
        self.pending = {}
    def ask(self, command):
        self.calls.append(command)
        if command['cmd'] == 'sessions':
            return {'sessions': [{'id': sid, 'name': 'same title', 'state': 'idle'} for sid in ['s1', 's2']]}
        sid = command['session']
        if sid in self.missing:
            return {'error': 'Sitzung nicht verfügbar'}
        action = command['action']
        if action == 'status':
            return {'ok': True, 'title': 'same title', 'state': 'idle', 'answer': 'short result'}
        if action == 'pending':
            return {'ok': True, 'items': self.pending.get(sid, [])}
        return {'ok': True, 'accepted': True}


@pytest.fixture
def rig(tmp_path, monkeypatch):
    monkeypatch.setenv('DELFIN_MESSENGER_DIR', str(tmp_path))
    store.save('7')
    store.put_room('#10', 's1', 'same title', owner='7', ready=True)
    store.put_room('#20', 's2', 'same title', owner='7', ready=True)
    store.put_room('#30', None, 'Forum', owner='7', ready=True)
    daemon, bridge = Daemon(), Bridge()
    return Rooms(daemon, bridge), daemon, bridge


def message(chat='#10', text='calculate', sender='7', mid='1'):
    return {'chat': chat, 'text': text, 'sender': sender, 'id': mid}


def test_chat_routes_to_exact_session_even_when_titles_match(rig):
    rooms, _, bridge = rig
    asyncio.run(rooms.handle(message()))
    assert bridge.calls == [{'cmd': 'session', 'session': 's1', 'action': 'prompt', 'text': 'calculate'}]


def test_forum_delivers_once_to_each_connected_session(rig):
    rooms, daemon, bridge = rig
    asyncio.run(rooms.handle(message('#30', 'report progress')))
    assert [r['session'] for r in bridge.calls] == ['s1', 's2']
    assert 's1' in daemon.sent[-1][1] and 's2' in daemon.sent[-1][1]


def test_replay_stays_refused_after_restart(rig):
    rooms, daemon, bridge = rig
    asyncio.run(rooms.handle(message('#30')))
    asyncio.run(Rooms(daemon, bridge).handle(message('#30')))
    assert len(bridge.calls) == 2


@pytest.mark.parametrize('chat,sender', [('#99', '7'), ('#10', '8'), ('@8', '8')])
def test_unknown_group_or_sender_is_ignored(rig, chat, sender):
    rooms, daemon, bridge = rig
    asyncio.run(rooms.handle(message(chat, sender=sender)))
    assert not bridge.calls and not daemon.sent


def test_old_direct_chat_no_longer_dispatches_tasks(rig):
    rooms, daemon, bridge = rig
    asyncio.run(rooms.handle(message('@7')))
    assert not bridge.calls
    assert 'Einrichtung' in daemon.sent[-1][1]


def test_revoke_stops_running_gateway_and_notifications(rig):
    rooms, daemon, bridge = rig
    store.revoke()
    asyncio.run(rooms.handle(message()))
    asyncio.run(rooms.refresh())
    assert not bridge.calls and not daemon.sent


def test_status_and_update_stay_in_the_room(rig):
    rooms, _, bridge = rig
    asyncio.run(rooms.handle(message(text='/status')))
    asyncio.run(rooms.handle(message(text='/update', mid='2')))
    assert bridge.calls == [
        {'cmd': 'session', 'session': 's1', 'action': 'status'},
        {'cmd': 'session', 'session': 's1', 'action': 'update', 'text': UPDATE}]


def test_forum_update_and_stop_only_target_connected_sessions(rig):
    rooms, _, bridge = rig
    asyncio.run(rooms.handle(message('#30', '/update')))
    asyncio.run(rooms.handle(message('#30', '/stop', mid='2')))
    assert [r['session'] for r in bridge.calls] == ['s1', 's2', 's1', 's2']
    assert [r['action'] for r in bridge.calls] == ['update', 'update', 'stop', 'stop']
    assert all(r['cmd'] == 'session' for r in bridge.calls)


def test_no_cross_session_redirect_or_global_approval(rig):
    rooms, daemon, bridge = rig
    asyncio.run(rooms.handle(message(text='@s2 do something')))
    asyncio.run(rooms.handle(message('#30', '/approve abc', mid='2')))
    assert not bridge.calls
    assert 'Einzelchat' in daemon.sent[-1][1]


def test_approval_answer_and_denial_carry_fixed_session(rig):
    rooms, _, bridge = rig
    for i, text in enumerate(['/approve abc', '/deny def', '/answer ghi DMF']):
        asyncio.run(rooms.handle(message(text=text, mid=str(i))))
    assert [r['action'] for r in bridge.calls] == ['approve', 'deny', 'answer']
    assert all(r['session'] == 's1' for r in bridge.calls)
    assert bridge.calls[-1]['text'] == 'DMF'


def test_missing_session_never_creates_or_uses_another(rig):
    rooms, daemon, bridge = rig
    bridge.missing.add('s1')
    asyncio.run(rooms.handle(message()))
    assert len(bridge.calls) == 1
    assert 'nicht verfügbar' in daemon.sent[-1][1]


def test_notifications_go_to_own_room_and_forum_only(rig):
    rooms, daemon, bridge = rig
    bridge.pending['s1'] = [{'id': 'abc', 'kind': 'confirm', 'text': 'run calculation', 'phoneApproval': True}]
    asyncio.run(rooms.refresh())
    s1 = [(chat, text) for chat, text in daemon.sent if '· s1]' in text]
    assert {chat for chat, _ in s1} == {'#10', '#30'}
    assert any('/approve abc' in text for _, text in s1)
    count = len(daemon.sent)
    asyncio.run(rooms.refresh())
    assert len(daemon.sent) == count
    assert all(r['cmd'] == 'session' and r['action'] in {'status', 'pending'} for r in bridge.calls)


def test_setup_is_idempotent_and_new_room_is_fixed(rig):
    rooms, daemon, _ = rig
    asyncio.run(rooms.sync())
    assert not daemon.groups
    asyncio.run(rooms.ensure_room('s3'))
    asyncio.run(rooms.ensure_room('s3'))
    assert len(daemon.groups) == len(daemon.invited) == 1
    with pytest.raises(ValueError, match='rebound'):
        store.put_room('#101', 's4', 'changed', owner='7')


def test_partial_broadcast_failure_is_reported_without_resending(rig):
    rooms, daemon, bridge = rig
    bridge.missing.add('s1')
    asyncio.run(rooms.handle(message('#30')))
    assert [r['session'] for r in bridge.calls] == ['s1', 's2']
    assert 'nicht verfügbar' in daemon.sent[-1][1] and 'Nachricht übergeben' in daemon.sent[-1][1]


def test_limits_and_no_automatic_answer_to_a_pending_question(rig):
    rooms, _, bridge = rig
    asyncio.run(rooms.handle(message(text='x' * 4097)))
    assert not bridge.calls
    for i in range(20):
        asyncio.run(rooms.handle(message(text='yes', mid=str(i))))
    assert len(bridge.calls) == 15
    assert all(r['action'] == 'prompt' for r in bridge.calls)


def test_disconnected_recipient_does_not_hide_other_broadcast_results(rig):
    rooms, daemon, bridge = rig
    original = bridge.ask
    def ask(command):
        if command.get('session') == 's1':
            bridge.calls.append(command)
            raise ConnectionError('private exception detail')
        return original(command)
    bridge.ask = ask
    asyncio.run(rooms.handle(message('#30')))
    assert [r['session'] for r in bridge.calls] == ['s1', 's2']
    assert 'Zustellung unklar' in daemon.sent[-1][1]
    assert 'private exception' not in daemon.sent[-1][1]
    assert 'Nachricht übergeben' in daemon.sent[-1][1]


def test_failed_notification_is_retried(rig):
    rooms, daemon, _ = rig
    original = daemon.send_chat
    async def fail(chat, text):
        raise ConnectionError('offline')
    daemon.send_chat = fail
    with pytest.raises(ConnectionError):
        asyncio.run(rooms.refresh())
    daemon.send_chat = original
    asyncio.run(rooms.refresh())
    assert {chat for chat, _ in daemon.sent} == {'#10', '#20', '#30'}


def test_pairing_rejects_an_unidentified_invitation(tmp_path, monkeypatch):
    from delfin.messenger import __main__ as cli
    monkeypatch.setenv('DELFIN_MESSENGER_DIR', str(tmp_path))
    class InvalidDaemon:
        async def __aenter__(self):
            return self
        async def __aexit__(self, *args):
            pass
        async def subscribe(self):
            pass
        async def pending_invitation(self):
            return 'simplex:/inv'
    monkeypatch.setattr(cli.simplex, 'Daemon', InvalidDaemon)
    assert cli._pair() == 1
    assert store.paired_contact() is None
