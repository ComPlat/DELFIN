"""Phone decisions are bound to an unchanged, still-open request and session."""
import asyncio
from types import SimpleNamespace

import pytest

from delfin.agent.service.bridge import Bridge
from delfin.agent.service.redact import Redactor
from delfin.agent.service.supervisor import Supervisor


@pytest.fixture
def scoped():
    redactor = Redactor()
    redactor._known = ['SECRET-CREDENTIAL']
    redactor._loaded = float('inf')
    bridge = Bridge(allowed_uid=1, redactor=redactor)
    requests = [{'requestId': 'p123-1', 'kind': 'confirm', 'tool': 'shell',
                 'preview': 'run calculation', 'options': [{'key': k} for k in ['y', 'n', 'A', 'e']]}]
    calls = []
    async def request(method, params=None):
        calls.append((method, params))
        if method == 'sessions.list':
            return {'sessions': [{'sessionId': s, 'title': 'same name'} for s in ['s1', 's2']]}
        if method == 'session.pending':
            return {'pending': list(requests)}
        if method == 'permission.reply':
            requests.clear()
            return {'applied': True}
        return {'ok': True}
    bridge._request = request
    return bridge, requests, calls


async def code_for(bridge):
    reply = await bridge._cmd_session({'session': 's1', 'action': 'pending'})
    return reply['items'][0]['id']


def test_approve_once_uses_the_existing_permission_gate(scoped):
    bridge, _, calls = scoped
    async def run():
        code = await code_for(bridge)
        assert (await bridge._cmd_session({'session': 's1', 'action': 'approve', 'id': code}))['ok']
        assert 'error' in await bridge._cmd_session({'session': 's1', 'action': 'approve', 'id': code})
    asyncio.run(run())
    decisions = [p for m, p in calls if m == 'permission.reply']
    assert decisions == [{'sessionId': 's1', 'requestId': 'p123-1', 'choice': 'y', 'persist': False}]


@pytest.mark.parametrize('change', ['session', 'changed', 'resolved', 'expired', 'text_instead'])
def test_wrong_stale_or_changed_request_cannot_be_approved(scoped, change):
    bridge, requests, calls = scoped
    async def run():
        code = await code_for(bridge)
        params = {'session': 's1', 'action': 'approve', 'id': code}
        if change == 'session':
            params['session'] = 's2'
        elif change == 'changed':
            requests[0]['preview'] = 'different action'
        elif change == 'resolved':
            requests.clear()
        elif change == 'expired':
            issued = bridge._issued[code]
            bridge._issued[code] = (*issued[:3], issued[3] - 601)
        else:
            params.update(action='answer', text='yes')
        assert 'error' in await bridge._cmd_session(params)
    asyncio.run(run())
    assert not any(m == 'permission.reply' for m, _ in calls)


def test_long_preview_cannot_be_approved_but_can_be_denied(scoped):
    bridge, requests, calls = scoped
    requests[0]['preview'] = 'x' * 1500
    async def run():
        reply = await bridge._cmd_session({'session': 's1', 'action': 'pending'})
        item = reply['items'][0]
        assert not item['phoneApproval']
        assert 'error' in await bridge._cmd_session({'session': 's1', 'action': 'approve', 'id': item['id']})
        assert (await bridge._cmd_session({'session': 's1', 'action': 'deny', 'id': item['id']}))['ok']
    asyncio.run(run())
    assert [p['choice'] for m, p in calls if m == 'permission.reply'] == ['n']


def test_plan_approval_cannot_enable_accept_edits(scoped):
    bridge, requests, calls = scoped
    requests[0].update(kind='plan', options=[{'key': k} for k in ['d', 'e', 'n']])
    async def run():
        code = await code_for(bridge)
        assert (await bridge._cmd_session({'session': 's1', 'action': 'approve', 'id': code}))['ok']
    asyncio.run(run())
    assert [p['choice'] for m, p in calls if m == 'permission.reply'] == ['d']


def test_question_text_is_not_a_permission_choice(scoped):
    bridge, requests, calls = scoped
    requests[0].update(kind='ask', question='Which solvent?')
    async def run():
        code = await code_for(bridge)
        assert 'error' in await bridge._cmd_session({'session': 's1', 'action': 'approve', 'id': code})
        assert (await bridge._cmd_session({'session': 's1', 'action': 'answer', 'id': code, 'text': 'DMF'}))['ok']
    asyncio.run(run())
    assert [p for m, p in calls if m == 'permission.reply'] == [{'sessionId': 's1', 'requestId': 'p123-1', 'text': 'DMF'}]


def test_phone_preview_redacts_credentials(scoped):
    bridge, requests, _ = scoped
    requests[0]['preview'] = 'command SECRET-CREDENTIAL'
    reply = asyncio.run(bridge._cmd_session({'session': 's1', 'action': 'pending'}))
    assert 'SECRET-CREDENTIAL' not in str(reply)
    assert '[redacted:credential]' in str(reply)


def test_unknown_session_does_not_fall_back(scoped):
    bridge, _, calls = scoped
    reply = asyncio.run(bridge._cmd_session({'session': 'same name', 'action': 'prompt', 'text': 'run'}))
    assert 'error' in reply
    assert [m for m, _ in calls] == ['sessions.list']


def test_pending_read_does_not_attach_wake_or_replay():
    service = SimpleNamespace(sessions={'s1': SimpleNamespace(pending={'r': {'requestId': 'r'}})})
    assert asyncio.run(Supervisor._m_session_pending(service, None, {'sessionId': 's1'})) == {'pending': [{'requestId': 'r'}]}
    assert asyncio.run(Supervisor._m_session_pending(service, None, {'sessionId': 'sleeping'})) == {'pending': []}


@pytest.mark.parametrize('prompt,action,text', [('bash', 'approve', ''), ('bash', 'deny', ''), ('ask', 'answer', 'water')])
def test_phone_decision_resolves_a_real_worker_request(tmp_path, prompt, action, text):
    from tui_service import ServiceProcess
    with ServiceProcess(tmp_path) as service:
        client = service.client()
        sid = service.session(client, mode='default')
        client.request('session.prompt', {'sessionId': sid, 'text': prompt})
        client.wait_for(lambda e: e.get('type') == 'permission.asked' and e.get('session') == sid)
        bridge = Bridge(allowed_uid=1, service_path=service.socket)
        async def run():
            pending = await bridge._run_command({'cmd': 'session', 'session': sid, 'action': 'pending'})
            assert pending['items'], pending
            code = pending['items'][0]['id']
            reply = await bridge._run_command({'cmd': 'session', 'session': sid, 'action': action, 'id': code, 'text': text})
            assert reply['ok'], reply
            repeat = await bridge._run_command({'cmd': 'session', 'session': sid, 'action': action, 'id': code, 'text': text})
            assert 'error' in repeat
        try:
            asyncio.run(run())
            event, _ = client.wait_for(lambda e: e.get('type') == 'permission.resolved' and e.get('session') == sid)
            assert event['decision'] == {'approve': 'approved', 'deny': 'denied', 'answer': 'answered'}[action]
            client.wait_for(lambda e: e.get('type') == 'turn.ended' and e.get('session') == sid)
        finally:
            if bridge._service is not None:
                bridge._service.close()
