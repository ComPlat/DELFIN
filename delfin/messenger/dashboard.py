"""Adapt a live dashboard tab to the bridge, using its existing UI callbacks."""
from __future__ import annotations
import secrets
import threading
import time
from functools import wraps


def gate(state):
    return state.setdefault('_messenger_gate_lock', threading.RLock())


def serialized_gate(state):
    """Use the same lock for local answers, phone answers and gate replacement."""
    def decorate(fn):
        @wraps(fn)
        def guarded(*args, **kwargs):
            with gate(state):
                return fn(*args, **kwargs)
        return guarded
    return decorate


class DashboardSession:
    def __init__(self, sid, state, *, send, stop, answer, plan, approval=None, question=None):
        self.sid, self.state = sid, state
        self.send, self.stop, self.answer, self.plan = send, stop, answer, plan
        self.approval = approval
        self.question = question
        self._requests = {}

    def _current(self):
        engine = self.state.get('engine')
        return str(getattr(engine, 'session_id', '') or self.state.get('active_session_id') or '')

    def pending(self):
        current = []
        if self.approval:
            action = self.state.get('_pending_dashboard_action')
            blocked = self.state.get('_pending_approval')
            obj = action or blocked
            if obj:
                preview = (action.get('description', '') if action else
                           str(blocked.get('detail', '')) + '\n' + str(blocked.get('tool', '')))
                current.append((obj, {'kind': 'confirm', 'tool': 'dashboard', 'preview': preview,
                                      'options': [{'key': 'y', 'label': 'Approve'},
                                                  {'key': 'n', 'label': 'Deny'}]}))
        broker = self.state.get('_kit_confirm_broker')
        if broker:
            for req in broker.pending_requests():
                if req.event.is_set():
                    continue
                pattern, kind, _ = broker.persistence_for(req)
                once = 'Grant permanently' if req.preview.startswith('⚠') else 'Allow (once)'
                options = [{'key': 'y', 'label': once, 'detail': 'Diese konkrete Aktion ausführen.'},
                           {'key': 'n', 'label': 'Deny', 'detail': 'Aktion ablehnen.'}]
                if pattern and broker._persist_callback is not None:
                    options.append({'key': 'p', 'label': 'Allow + Permanent',
                                    'detail': f'{kind}: {pattern}; gilt auch für spätere Sitzungen.'})
                current.append((req, {'kind': 'confirm', 'tool': req.tool_name,
                                      'preview': req.preview, 'options': options}))
        ev = self.state.get('_ask_user_event')
        detected = self.state.get('_pending_question')
        if self.question and detected and not (ev and not ev.is_set()):
            choices = [f'{number}: {label}' for number, label in detected.get('options', [])]
            current.append((detected, {'kind': 'ask', 'question': self.state.get('_last_answer_text', '')[-1400:],
                                       'choices': choices, 'detected': True, 'options': []}))
        if ev and not ev.is_set() and time.monotonic() < self.state.get('_ask_user_deadline', float('inf')):
            args = self.state.get('_ask_user_args') or {}
            choices = [o.get('label', '') if isinstance(o, dict) else str(o) for o in args.get('options', [])]
            current.append((ev, {'kind': 'ask', 'question': args.get('question', ''), 'choices': choices,
                                 'choiceDetails': [str(o.get('description') or o.get('preview') or '') if isinstance(o, dict) else ''
                                                   for o in args.get('options', [])],
                                 'preview': '', 'options': []}))
        ev = self.state.get('_plan_approval_event')
        if (ev and not ev.is_set() and not (self.state.get('_plan_approval_result') or {}).get('timed_out')
                and time.monotonic() < self.state.get('_plan_approval_deadline', float('inf'))):
            current.append((ev, {'kind': 'plan', 'preview': self.state.get('_pending_plan_body', ''),
                                 'options': [{'key': 'e', 'label': 'Accept plan & execute',
                                              'detail': 'Plan ausführen; Accept Edits aktivieren.'},
                                             {'key': 'n', 'label': 'Reject', 'detail': 'Plan ablehnen.'}]}))
        self._requests = {key: value for key, value in self._requests.items()
                          if any(value[0] is obj for obj, _ in current)}
        result = []
        for obj, fields in current:
            key = next((k for k, v in self._requests.items() if v[0] is obj), None)
            if key is None:
                key = secrets.token_hex(12)
            self._requests[key] = (obj, fields)
            result.append({'requestId': key, **fields})
        return result

    def __call__(self, command):
        with gate(self.state):
            return self._dispatch(command)

    def _dispatch(self, command):
        if self.state.get('_closed') or self._current() != self.sid:
            return {'error': 'Diese Dashboard-Sitzung ist nicht mehr geöffnet.'}
        method = command['method']
        if method == 'session.setTitle':
            from delfin.agent.service.protocol import clean_title
            from delfin.agent.session_store import rename_session
            title = clean_title(command.get('title'))
            if not title:
                return {'error': 'Bitte einen Namen angeben.'}
            rename_session(self.sid, title)
            self.state.setdefault('_messenger_titles', {})[self.sid] = title
            return {'title': title}
        if method == 'session.status':
            return {'sessionId': self.sid, 'title': self.state.get('_messenger_titles', {}).get(self.sid) or self.sid,
                    'state': 'running' if self.state.get('streaming') else 'idle',
                    'workspace': str(getattr(getattr(self.state.get('engine'), 'kit_permissions', None), 'workspace', '') or ''),
                    'answer': str(self.state.get('_last_answer_text') or '')[:1000], 'pending': len(self.pending())}
        if method == 'session.pending':
            return {'pending': self.pending()}
        if method == 'session.cancel':
            self.stop()
            return {'ok': True}
        if method == 'session.prompt':
            accepted = self.send(command.get('text', ''))
            return {'accepted': bool(accepted), **({} if accepted else {'error': 'Dashboard hat eine Eingabe offen; Nachricht nicht gesendet.'})}
        if method == 'permission.reply':
            self.pending()
            entry = self._requests.get(command.get('requestId'))
            if not entry:
                return {'applied': False}
            obj, fields = entry
            choice = command.get('choice')
            if fields['kind'] == 'ask':
                if fields.get('detected'):
                    return {'applied': bool(self.question(command.get('text', '')))}
                return {'applied': bool(self.answer(command.get('text', '')))}
            if choice not in {o['key'] for o in fields['options']}:
                return {'applied': False}
            if fields['kind'] == 'plan':
                self.plan(choice != 'n')
                return {'applied': obj.is_set()}
            if fields.get('tool') == 'dashboard':
                self.approval(choice == 'y')
                return {'applied': not any(self.state.get(k) is obj for k in
                                          ('_pending_dashboard_action', '_pending_approval'))}
            broker = self.state['_kit_confirm_broker']
            pattern, kind, _ = broker.persistence_for(obj)
            if choice == 'p' and not command.get('persist'):
                return {'applied': False}
            return {'applied': broker.resolve(obj, choice != 'n',
                        persist_pattern=pattern if choice == 'p' else None,
                        persist_kind=kind if choice == 'p' else None)}
        return {'error': 'Unsupported dashboard command'}
