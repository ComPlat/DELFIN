"""Adapt a live dashboard tab to the bridge, using its existing UI callbacks."""
from __future__ import annotations
import re
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
                options = [{'key': 'y', 'label': once, 'detail': 'Run this one action.'},
                           {'key': 'n', 'label': 'Deny', 'detail': 'Aktion ablehnen.'}]
                if pattern and broker._persist_callback is not None:
                    options.append({'key': 'p', 'label': 'Allow + Permanent',
                                    'detail': f'{kind}: {pattern}; also for later sessions.'})
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
                                              'detail': 'Run the plan; switch to Accept Edits.'},
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
            return {'error': 'This dashboard session is no longer open.'}
        method = command['method']
        if method == 'session.setTitle':
            from delfin.agent.service.protocol import clean_title
            from delfin.agent.session_store import rename_session
            title = clean_title(command.get('title'))
            if not title:
                return {'error': 'Please give a name.'}
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
            text = str(command.get('text', ''))
            if text.lstrip()[:1] in {'/', '!'}:
                # A phone's line is never a dashboard command or shell line.
                return {'accepted': False, 'error': 'A command from the phone is not run here.'}
            accepted = self.send(text)
            return {'accepted': bool(accepted), **({} if accepted else {'error': 'The dashboard has an input open; '
                                                                                 'message not sent.'})}
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


class PairPanel:
    """The dashboard's phone-pairing panel: the invitation (QR code and
    link) and the waiting-for-the-code step shown in a widget of their own,
    outside the chat transcript and the saved session; removed on success,
    expiry and cancel. Everything else /pair says goes to the chat (``say``),
    and none of that carries the link or the security code.

    ``render(html)``: sets the panel widget ("" hides it). ``on_stage``:
    told the step the panel is at -- "offer" (scan the code), "code" (type
    the code the phone got; the dashboard shows its code field) or None
    (closed). ``ttl_s``: when the panel is removed at the latest if no
    answer comes."""

    def __init__(self, *, render, say, ttl_s=None, clock=time.monotonic, on_stage=None):
        from .pairing import INVITATION_TTL_S
        self.render, self.say, self.clock = render, say, clock
        self.on_stage = on_stage or (lambda stage: None)
        self.ttl_s = INVITATION_TTL_S + 30 if ttl_s is None else float(ttl_s)
        self._deadline = None
        self._lock = threading.Lock()

    @staticmethod
    def transient(text):
        """Whether ``text`` belongs in the panel: it shows a link, a QR code
        or the security code (the invitation and the connected step)."""
        text = str(text or '')
        return ('simplex:/' in text or 'simplex.chat/' in text
                or text.startswith('Phone connected') or 'security code' in text
                or text.startswith('Wrong code'))

    _LINK = re.compile(r'(?:simplex:/|https://simplex\.chat/)\S+')

    def show(self, text):
        if not self.transient(text):
            self.clear()
            if text:
                self.say(text)
            return
        import html
        text = str(text)
        with self._lock:
            self._deadline = self.clock() + self.ttl_s
        link = self._LINK.search(text)
        if link:
            from . import qr
            # The steps and the link as text, the code as a picture: the
            # text code's block lines are dropped beside it.
            words = [line for line in text.splitlines() if not set(line) <= set('█▀▄ ')]
            try:
                picture = qr.svg(link.group(0), scale=4)
            except Exception:  # noqa: BLE001 - the link below still works
                picture = ''
            body = ('<div style="display:flex;gap:12px;flex-wrap:wrap;align-items:flex-start">'
                    f'<div>{picture}</div><pre style="white-space:pre-wrap;margin:4px 0;flex:1 1 260px">'
                    + html.escape('\n'.join(words)) + '</pre></div>')
            stage = 'offer'
        else:
            body = ('<pre style="white-space:pre-wrap;margin:4px 0">' + html.escape(text) + '</pre>'
                    '<div><b>Type the 6-digit code from your phone in the field below.</b></div>')
            stage = 'code'
        self.render('<div class="delfin-messenger-pair">' + body + '</div>')
        self.on_stage(stage)
        if self.ttl_s > 0:
            timer = threading.Timer(self.ttl_s, self.expire_if_due)
            timer.daemon = True
            timer.start()

    def clear(self):
        with self._lock:
            self._deadline = None
        self.render('')
        self.on_stage(None)

    def expire_if_due(self):
        with self._lock:
            due = self._deadline is not None and self.clock() >= self._deadline
        if due:
            self.clear()
