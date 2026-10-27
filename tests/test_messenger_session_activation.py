"""Opt-in and phone decisions use the same live dashboard gates."""
import threading
from types import SimpleNamespace

import pytest

from delfin.messenger import activation
from delfin.messenger.dashboard import DashboardSession
from delfin.agent.kit_confirm import KitConfirmBroker, _ConfirmRequest


def test_activation_is_per_session_and_off_invalidates_old_requests(tmp_path, monkeypatch):
    monkeypatch.setattr(activation, '_DIR', tmp_path)
    assert not activation.enabled('one')
    activation.configure('one')
    first = activation.record('one')['generation']
    assert activation.enabled('one') and not activation.enabled('two')
    activation.configure('one', 'on')
    assert activation.record('one')['generation'] == first
    activation.configure('one', 'off')
    assert not activation.enabled('one')
    activation.configure('one')
    assert activation.record('one')['generation'] != first


def adapter(state, **callbacks):
    callbacks = {**dict(send=lambda text: True, stop=lambda: None,
                       answer=lambda text: False, plan=lambda yes: None), **callbacks}
    return DashboardSession('one', state, **callbacks)


@pytest.mark.parametrize('choice,expected', [('y', True), ('n', False)])
def test_phone_clears_actual_dashboard_broker_and_refreshes_widget(choice, expected):
    broker = KitConfirmBroker()
    req = _ConfirmRequest(1, 'bash', {'command': 'echo hi'}, 'echo hi')
    broker._pending.append(req)
    refreshed = []
    broker._refresh_panel = lambda: refreshed.append(True)
    ui = adapter({'active_session_id': 'one', '_kit_confirm_broker': broker})
    item = ui({'method': 'session.pending'})['pending'][0]
    reply = {'method': 'permission.reply', 'requestId': item['requestId'], 'choice': choice}
    assert ui(reply)['applied']
    assert req.decision is expected and req.event.is_set()
    assert broker.pending_requests() == [] and refreshed
    assert ui({'method': 'session.pending'}) == {'pending': []}
    assert not ui(reply)['applied']


def test_phone_closes_question_and_plan_using_dashboard_callbacks():
    question, plan = threading.Event(), threading.Event()
    state = {'active_session_id': 'one', '_ask_user_event': question,
             '_ask_user_args': {'question': 'Solvent?', 'options': [{'label': 'Water'}]},
             '_plan_approval_event': plan, '_pending_plan_body': 'Run calculation'}
    answers = []
    def answer(text):
        answers.append(text)
        question.set()
        return True
    ui = adapter(state, answer=answer, plan=lambda yes: plan.set())
    pending = ui.pending()
    assert pending[0]['choices'] == ['Water']
    assert ui({'method': 'permission.reply', 'requestId': pending[0]['requestId'], 'text': 'Water'})['applied']
    assert ui({'method': 'permission.reply', 'requestId': pending[1]['requestId'], 'choice': 'e'})['applied']
    assert answers == ['Water'] and ui.pending() == []


def test_dashboard_action_uses_existing_approve_and_deny_handler():
    state = {'active_session_id': 'one', '_pending_dashboard_action': {'description': 'Cancel job 123'}}
    decisions = []
    def decide(yes):
        decisions.append(yes)
        state['_pending_dashboard_action'] = None
    ui = adapter(state, approval=decide)
    item = ui.pending()[0]
    assert item['preview'] == 'Cancel job 123'
    assert ui({'method': 'permission.reply', 'requestId': item['requestId'], 'choice': 'n'})['applied']
    assert decisions == [False] and ui.pending() == []


def test_switched_dashboard_never_answers_another_session():
    ui = adapter({'active_session_id': 'two'})
    assert 'error' in ui({'method': 'session.pending'})


@pytest.fixture
def short_socket_runtime(monkeypatch):
    import tempfile
    with tempfile.TemporaryDirectory(prefix='dmsg-', dir='/tmp') as path:
        monkeypatch.setenv('DELFIN_RUNTIME_DIR', path)
        yield path


def test_local_dashboard_socket_obeys_activation(tmp_path, monkeypatch, short_socket_runtime):
    from delfin.agent.service.supervisor import Supervisor
    monkeypatch.setattr(activation, '_DIR', tmp_path / 'activation')
    # The real gateway is a separate user-service process, not this dashboard process.
    monkeypatch.setattr(Supervisor, '_from_an_agent', lambda *_: False)
    server = activation.LocalSession('one', lambda command: {'pending': []})
    try:
        activation.configure('one', endpoint=str(server.path))
        assert activation.ask_dashboard('one', 'session.pending', {}) == {'pending': []}
        activation.configure('one', 'off')
        assert 'error' in activation.ask_dashboard('one', 'session.pending', {})
    finally:
        server.close()


def test_phone_resolution_hides_real_dashboard_widget():
    pytest.importorskip('ipywidgets')
    broker = KitConfirmBroker()
    req = _ConfirmRequest(1, 'bash', {'command': 'echo hi'}, 'echo hi')
    broker._pending.append(req)
    panel = broker.build_widget()
    assert panel.layout.display != 'none'
    ui = adapter({'active_session_id': 'one', '_kit_confirm_broker': broker})
    request = ui.pending()[0]
    assert ui({'method': 'permission.reply', 'requestId': request['requestId'], 'choice': 'n'})['applied']
    assert panel.layout.display == 'none' and broker._requests_box.children == ()


def test_replaced_dashboard_question_rejects_old_phone_answer():
    state = {'active_session_id': 'one', '_pending_question': {'type': 'open', 'options': []},
             '_last_answer_text': 'First question?'}
    answers = []
    ui = adapter(state, question=lambda text: answers.append(text) or True)
    original = ui.pending()[0]['requestId']
    state['_pending_question'] = {'type': 'open', 'options': []}
    state['_last_answer_text'] = 'Second question?'
    assert not ui({'method': 'permission.reply', 'requestId': original, 'text': 'yes'})['applied']
    assert answers == []


def test_activation_store_is_hidden_from_agent_files_and_shells(tmp_path, monkeypatch):
    # Masked as a home secret in every cage (like the trust store), not as
    # a runtime directory: a workspace that contains ~/.delfin keeps its
    # write root under Landlock.
    from delfin.agent import boundaries, sandbox, api_client
    monkeypatch.setattr(activation, '_DIR', tmp_path / 'messenger_sessions')
    assert '.delfin/messenger_sessions' in sandbox._HOME_SECRET_DIRS
    assert '.delfin/messenger_sessions' in boundaries.MASKED_HOME_ENTRIES
    assert str(activation._DIR.resolve()) not in sandbox.runtime_masked_dirs()
    assert not activation._DIR.exists()
    import fnmatch
    assert any(fnmatch.fnmatch('messenger_sessions/a.json', glob)
               for glob in api_client._DEFAULT_PATH_DENY_GLOBS)


def test_resuming_in_tui_replaces_a_previous_dashboard_endpoint(tmp_path, monkeypatch):
    monkeypatch.setattr(activation, '_DIR', tmp_path)
    activation.configure('one', endpoint='/tmp/previous-dashboard.sock')
    generation = activation.record('one')['generation']
    activation.configure('one', endpoint='')
    assert activation.enabled('one') and 'endpoint' not in activation.record('one')
    assert activation.record('one')['generation'] != generation


@pytest.mark.parametrize('kind', ['question', 'plan'])
def test_elapsed_dashboard_deadline_rejects_phone_decision(kind):
    import time
    state = {'active_session_id': 'one'}
    event = threading.Event()
    if kind == 'question':
        state.update(_ask_user_event=event, _ask_user_args={'question': 'Proceed?'})
        deadline = '_ask_user_deadline'
    else:
        state.update(_plan_approval_event=event, _pending_plan_body='First plan')
        deadline = '_plan_approval_deadline'
    calls = []
    ui = adapter(state, answer=lambda text: calls.append(text), plan=lambda yes: calls.append(yes))
    original = ui.pending()[0]['requestId']
    state[deadline] = time.monotonic() - 1
    assert not ui({'method': 'permission.reply', 'requestId': original, 'text': 'yes', 'choice': 'e'})['applied']
    assert calls == [] and not event.is_set()


def test_dashboard_gate_cannot_be_replaced_during_phone_resolution():
    from delfin.messenger.dashboard import serialized_gate
    state = {'active_session_id': 'one', '_pending_dashboard_action': {'description': 'First'}}
    started, replaced = threading.Event(), threading.Event()
    @serialized_gate(state)
    def replace():
        started.set()
        state['_pending_dashboard_action'] = {'description': 'Second'}
        replaced.set()
    threads = []
    def decide(yes):
        thread = threading.Thread(target=replace)
        threads.append(thread)
        thread.start()
        assert not replaced.wait(0.05)
        assert state['_pending_dashboard_action']['description'] == 'First'
        state['_pending_dashboard_action'] = None
    ui = adapter(state, approval=decide)
    item = ui.pending()[0]
    assert ui({'method': 'permission.reply', 'requestId': item['requestId'], 'choice': 'y'})['applied']
    threads[0].join(2)
    assert replaced.is_set() and state['_pending_dashboard_action']['description'] == 'Second'
