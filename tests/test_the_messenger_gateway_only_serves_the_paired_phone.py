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
        return {"chatItems": [{"chatItem": {"meta": {"itemId": len(self.sent)}}}]}
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
        self.active = {"s1", "s2", "s3"}
        self.reads = []
    def ask(self, command):
        if command.get("action") not in {None, "status", "pending"}:
            self.calls.append(command)
        else:
            self.reads.append(command)
        if command['cmd'] == 'sessions':
            return {'sessions': [{'id': sid, 'name': 'same title', 'state': 'idle', 'messengerEnabled': sid in self.active} for sid in ['s1', 's2']]}
        sid = command['session']
        if sid in self.missing:
            return {'error': 'Sitzung nicht verfügbar'}
        action = command['action']
        if action == 'status':
            return {'ok': True, 'title': 'same title', 'state': 'idle', 'answer': 'short result', 'messengerEnabled': sid in self.active}
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
    assert not daemon.sent


def test_notifications_go_to_own_room_and_forum_only(rig):
    rooms, daemon, bridge = rig
    bridge.pending['s1'] = [{'id': 'abc', 'kind': 'confirm', 'text': 'run calculation', 'phoneApproval': True, 'sessionId': 's1', 'options': [{'key': 'y', 'label': 'yes'}, {'key': 'n', 'label': 'no'}]}]
    asyncio.run(rooms.refresh())
    s1 = [(chat, text) for chat, text in daemon.sent if '· s1]' in text]
    assert {chat for chat, _ in s1} == {'#10', '#30'}
    assert any('yes/ja' in text for _, text in s1)
    count = len(daemon.sent)
    asyncio.run(rooms.refresh())
    assert len(daemon.sent) == count
    assert not bridge.calls


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
        asyncio.run(rooms.handle(message(text='calculate', mid=str(i))))
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


def test_sync_only_connects_explicitly_enabled_sessions(rig):
    rooms, daemon, bridge = rig
    bridge.active.clear()
    store.forget_room('#10')
    store.forget_room('#20')
    store.forget_room('#30')
    asyncio.run(rooms.sync())
    assert not daemon.groups
    bridge.active.add('s2')
    asyncio.run(rooms.sync())
    assert {r['session_id'] for r in store.rooms().values()} == {None, 's2'}


def test_disabled_room_cannot_receive_commands_notifications_or_broadcasts(rig):
    rooms, daemon, bridge = rig
    bridge.active.remove('s1')
    asyncio.run(rooms.handle(message('#10', 'calculate')))
    assert not bridge.calls and not daemon.sent
    asyncio.run(rooms.handle(message('#30', 'calculate', mid='2')))
    assert [r['session'] for r in bridge.calls] == ['s2']
    asyncio.run(rooms.refresh())
    assert not any(chat == '#10' for chat, _ in daemon.sent)


def test_reply_yes_uses_the_quoted_request_without_slash_or_code(rig):
    rooms, daemon, bridge = rig
    item = {'id': 'abc', 'sessionId': 's1', 'kind': 'confirm', 'text': 'Run tests?',
            'phoneApproval': True, 'options': [{'key': 'y', 'label': 'yes'}, {'key': 'n', 'label': 'no'}]}
    bridge.pending['s1'] = [item]
    asyncio.run(rooms.refresh())
    quoted = next(mid for (chat, mid), value in rooms._reply_targets.items() if chat == '#10' and value['id'] == 'abc')
    msg = message(text='yes')
    msg['quote_id'] = quoted
    asyncio.run(rooms.handle(msg))
    assert bridge.calls == [{'cmd': 'session', 'session': 's1', 'action': 'decide',
                             'id': 'abc', 'choice': 'y', 'persist': False}]
    assert 'TUI/Dashboard erledigt' in daemon.sent[-1][1]


def test_unquoted_yes_and_wrong_quote_cannot_approve(rig):
    rooms, _, bridge = rig
    asyncio.run(rooms.handle(message(text='yes')))
    msg = message(text='yes', mid='2')
    msg['quote_id'] = 'other-message'
    asyncio.run(rooms.handle(msg))
    assert not bridge.calls


def test_plain_question_answer_is_bound_to_its_message(rig):
    rooms, _, bridge = rig
    asyncio.run(rooms.send_request('#10', store.rooms()['#10'],
        {'id': 'question-code', 'sessionId': 's1', 'kind': 'question', 'text': 'Solvent?'}))
    msg = message(text='DMF')
    msg['quote_id'] = next(mid for chat, mid in rooms._reply_targets if chat == '#10')
    asyncio.run(rooms.handle(msg))
    assert bridge.calls == [{'cmd': 'session', 'session': 's1', 'action': 'answer', 'id': 'question-code', 'text': 'DMF'}]


def test_forum_and_notifications_can_be_disabled_independently(rig):
    rooms, daemon, bridge = rig
    asyncio.run(rooms.handle(message(text='forum off')))
    asyncio.run(rooms.handle(message('#30', 'calculate', mid='2')))
    assert [r['session'] for r in bridge.calls] == ['s2']
    asyncio.run(rooms.handle(message(text='meldungen aus', mid='3')))
    daemon.sent.clear()
    asyncio.run(rooms.refresh())
    assert all(chat != '#10' for chat, _ in daemon.sent)
    asyncio.run(rooms.handle(message('#30', 'forum off', mid='4')))
    bridge.calls.clear()
    asyncio.run(rooms.handle(message('#30', 'calculate', mid='5')))
    assert bridge.calls == []


def test_agent_announcements_are_opt_in_and_never_replayed(rig):
    rooms, daemon, bridge = rig
    ask = bridge.ask
    def answer(command):
        if command.get('action') == 'announcement':
            return {'text': 'Ready', 'id': 'unique'} if command['session'] == 's1' else {'text': ''}
        return ask(command)
    bridge.ask = answer
    asyncio.run(rooms.refresh())
    assert bridge.calls == []
    asyncio.run(rooms.handle(message('#30', 'austausch on')))
    asyncio.run(rooms.refresh())
    assert len(bridge.calls) == 1
    assert bridge.calls[0]['session'] == 's2'
    assert 'keine Nutzeranweisung oder Freigabe' in bridge.calls[0]['text']
    asyncio.run(rooms.refresh())
    assert len(bridge.calls) == 1


def test_expired_quote_never_approves_the_next_request(rig):
    rooms, daemon, bridge = rig
    original = {'id': 'old', 'sessionId': 's1', 'kind': 'confirm', 'text': 'First',
                'phoneApproval': True, 'options': [{'key': 'y'}, {'key': 'n'}]}
    bridge.pending['s1'] = [original]
    async def run():
        await rooms.refresh()
        old_mid = next(mid for (chat, mid), item in rooms._reply_targets.items() if item['id'] == 'old')
        bridge.pending['s1'] = [{**original, 'id': 'new', 'text': 'Second'}]
        await rooms.refresh()
        await rooms.handle({**message(text='yes', mid='new-message'), 'quote_id': old_mid})
    asyncio.run(run())
    assert bridge.calls == []
    assert 'nicht mehr zugeordnet' in daemon.sent[-1][1]


def test_voice_transcript_is_reviewed_and_never_answers_a_permission(rig):
    rooms, _, bridge = rig
    async def run():
        await rooms.send('#10', 'Transcript: yes', request={
            'id': 'voice', 'sessionId': 's1', 'kind': 'voice', 'text': 'yes', 'generation': None})
        mid = next(mid for chat, mid in rooms._reply_targets if chat == '#10')
        await rooms.handle({**message(text='yes', mid='voice-yes'), 'quote_id': mid})
        assert bridge.calls == []
        await rooms.handle({**message(text='senden', mid='voice-send'), 'quote_id': mid})
        assert bridge.calls == [{'cmd': 'session', 'session': 's1', 'action': 'prompt', 'text': 'Sprachnachricht des Nutzers (Transkript):\nyes'}]
    asyncio.run(run())


def test_permanent_permission_needs_second_message_reply(rig):
    rooms, _, bridge = rig
    base = bridge.ask
    def ask(command):
        if command.get('action') == 'decide' and not command.get('persist'):
            bridge.calls.append(command)
            return {'confirmationRequired': True}
        return base(command)
    bridge.ask = ask
    async def run():
        await rooms.send_request('#10', store.rooms()['#10'], {'id': 'rule', 'sessionId': 's1',
            'kind': 'confirm', 'text': 'echo hi', 'phoneApproval': True,
            'options': [{'key': 'y', 'label': 'Once'}, {'key': 'A', 'label': 'Always', 'detail': 'echo hi'}]})
        mid = next(iter(rooms._reply_targets))[1]
        await rooms.handle({**message(text='Always'), 'quote_id': mid})
        assert not bridge.calls[-1]['persist']
        confirmation = next(mid for (chat, mid), item in rooms._reply_targets.items() if item.get('confirmChoice'))
        await rooms.handle({**message(text='yes', mid='second'), 'quote_id': confirmation})
        assert bridge.calls[-1]['persist'] and bridge.calls[-1]['choice'] == 'A'
    asyncio.run(run())
