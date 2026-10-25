"""A private group per session, plus an explicitly labelled shared forum.

Only the paired human's received messages are commands. Agent notifications
are never fed back into other agents. Stable group/session IDs determine the
recipient; names are presentation only.
"""
from __future__ import annotations

import asyncio
import time
from collections import deque

from . import simplex, store
from .bridge_client import BridgeClient, default_socket_path

HELP = (
    'Einzelchat: Text geht nur an diese Sitzung.\n'
    'Forum: Text geht an alle verbundenen Sitzungen.\n'
    '/status – Stand; /update – kurzen Bericht anfordern; /pending – offene Anfragen\n'
    '/answer CODE Text – Rückfrage beantworten\n'
    '/approve CODE oder /deny CODE – konkrete Anfrage im Einzelchat entscheiden\n'
    '/stop – aktuellen Lauf anhalten (im Forum: verbundene Sitzungen)\n'
    '/sessions – Sitzungen; /chat SITZUNGS-ID – Einzelchat einrichten\n'
    '/forum – Forum einrichten; /sync – aktuelle Sitzungen verbinden\n'
    'Freigaben gelten einmalig; Planfreigaben wechseln in den Standardmodus.'
)
UPDATE = 'Bitte gib mir einen kurzen Zwischenstand: erledigt, aktuell in Arbeit, Hindernisse, nächster Schritt. Höchstens fünf Sätze.'


def label(room: dict) -> str:
    return f"{room['title']} · {room['session_id']}"


class Rooms:
    def __init__(self, daemon, bridge=None):
        self.daemon = daemon
        self.bridge = bridge or BridgeClient(default_socket_path())
        self.owner = store.paired_contact()
        self._rate = deque()
        self._shown = {}
        self._setup = asyncio.Lock()

    def paired(self):
        return bool(self.owner) and store.paired_contact() == self.owner

    async def ask(self, command):
        if not self.paired():
            raise ValueError('Handy wurde entkoppelt.')
        return await asyncio.to_thread(self.bridge.ask, command)

    async def session(self, sid, action, **kwargs):
        return await self.ask({'cmd': 'session', 'session': sid, 'action': action, **kwargs})

    async def ensure_room(self, sid=None):
        async with self._setup:
            if not self.paired():
                raise ValueError('Handy wurde entkoppelt.')
            existing = next(((chat, r) for chat, r in store.rooms().items()
                             if r['session_id'] == sid), None)
            if existing and existing[1].get('ready'):
                return existing[0]
            if sid is not None:
                reply = await self.session(sid, 'status')
                if reply.get('error'):
                    raise ValueError(reply['error'])
                title = str(reply.get('title') or sid)
                name = f'DELFIN {title[:40]} {sid[:12]}'
            else:
                title, name = 'DELFIN Forum – alle verbundenen Sitzungen', 'DELFIN Forum'
            chat = existing[0] if existing else await self.daemon.create_group(name)
            store.put_room(chat, sid, title, owner=self.owner)
            await self.daemon.invite(chat, self.owner)
            store.put_room(chat, sid, title, owner=self.owner, ready=True)
            return chat

    async def sync(self):
        await self.ensure_room()
        reply = await self.ask({'cmd': 'sessions'})
        if reply.get('error'):
            raise ValueError(reply['error'])
        for row in reply.get('sessions', []):
            if row.get('id') and not row.get('elsewhere'):
                await self.ensure_room(row['id'])

    async def send(self, chat, text):
        if not self.paired():
            return
        if chat != '@' + self.owner and not store.rooms().get(chat, {}).get('ready'):
            return
        # Split at message boundaries; never truncate an approval preview silently.
        for start in range(0, len(text), 2000):
            if self.paired():
                await self.daemon.send_chat(chat, text[start:start + 2000])

    async def handle(self, message):
        if not self.paired() or message['sender'] != self.owner:
            return
        chat = message['chat']
        room = store.rooms().get(chat)
        direct = chat == '@' + self.owner
        if not direct and (not room or not room.get('ready')):
            return
        text = message['text'].strip()
        now = time.monotonic()
        while self._rate and now - self._rate[0] > 60:
            self._rate.popleft()
        if len(self._rate) >= 15 or len(text) > 4096 or not text:
            return
        self._rate.append(now)
        if not store.claim_message(chat, message['id']):
            return
        try:
            result = await self.dispatch(text, room)
        except Exception:
            # No raw service exceptions or credentials in phone error messages.
            result = 'Nicht ausgeführt oder Verbindung unterbrochen. Bitte /status prüfen; nicht automatisch erneut gesendet.'
        await self.send(chat, result)

    async def dispatch(self, text, room):
        word, _, rest = text.partition(' ')
        command = word.lower().lstrip('/')
        if command in {'help', 'hilfe'}:
            return HELP
        if command in {'sessions', 'sitzungen'}:
            reply = await self.ask({'cmd': 'sessions'})
            return reply.get('error') or '\n'.join(
                f"{r['id']} · {r.get('name') or r['id']} · {r.get('state')}" for r in reply.get('sessions', [])) or 'Keine Sitzungen.'
        if command == 'forum':
            await self.ensure_room()
            return 'Einladung zum Forum liegt im SimpleX-Chat. Dort geht freier Text an alle verbundenen Sitzungen.'
        if command == 'sync':
            await self.sync()
            return 'Forum und Einzelchats eingerichtet. Bitte die Gruppeneinladungen annehmen.'
        if command == 'chat' and rest.strip():
            await self.ensure_room(rest.strip())
            return 'Einladung zum festen Sitzungs-Chat gesendet. Bitte in SimpleX annehmen.'
        if room is None:
            return 'Dieser Chat dient nur der Einrichtung. /sessions zeigt IDs, /chat ID verbindet eine Sitzung, /sync richtet Forum und Einzelchats ein.'
        sid = room['session_id']
        if command in {'approve', 'deny', 'answer'}:
            if sid is None:
                return 'Bitte die Anfrage im zugehörigen Einzelchat beantworten oder freigeben. Im Forum gibt es keine Sammelfreigabe.'
            code, _, answer = rest.strip().partition(' ')
            if not code or (command == 'answer' and not answer):
                return 'Syntax: /approve CODE, /deny CODE oder /answer CODE Text. /pending zeigt die Codes.'
            reply = await self.session(sid, command, id=code, **({'text': answer} if command == 'answer' else {}))
            return reply.get('error') or ('Entscheidung übernommen.' if reply.get('ok') else 'Nicht übernommen.')
        if text.startswith('@'):
            return 'Dieser Chat hat feste Empfänger. Bitte den Einzelchat der gewünschten Sitzung öffnen.'
        if text.startswith('/') and command not in {'status', 'pending', 'stop', 'update'}:
            return 'Unbekannter Befehl. /help zeigt die verfügbaren Befehle.'
        targets = [room] if sid else [r for r in store.rooms().values() if r['session_id'] and r.get('ready')]
        if not targets:
            return 'Keine Sitzungen verbunden. /sync verbindet die aktuellen Sitzungen.'
        results = []
        for target in targets:
            # Revalidate before each dispatch: a revoke must stop a broadcast too.
            if not self.paired():
                break
            target_sid = target['session_id']
            try:
                if command == 'pending':
                    reply = await self.session(target_sid, 'pending')
                    result = reply.get('error') or self.pending_text(reply.get('items', []))
                elif command == 'status':
                    reply = await self.session(target_sid, 'status')
                    result = reply.get('error') or '\n'.join(filter(None, [reply.get('state'), reply.get('activity'), reply.get('answer')]))
                else:
                    action = 'stop' if command == 'stop' else 'update' if command == 'update' else 'prompt'
                    task = UPDATE if command == 'update' else text
                    reply = await self.session(target_sid, action, **({'text': task} if action in {'prompt', 'update'} else {}))
                    result = reply.get('error') or ('Lauf wird angehalten.' if action == 'stop' else 'Nachricht übergeben.')
            except Exception:
                result = 'Zustellung unklar: Verbindung unterbrochen. Vor erneutem Senden /status prüfen.'
            results.append(f'[{label(target)}]\n{result}')
        return '\n\n'.join(results)

    @staticmethod
    def pending_text(items):
        lines = []
        for item in items:
            code, kind = item['id'], item.get('kind')
            text = f"{kind} · {code}\n{item.get('tool') or ''}\n{item.get('text') or ''}"
            if kind == 'question':
                text += f'\n/answer {code} Antworttext'
            elif kind in {'confirm', 'plan'}:
                if item.get('phoneApproval'):
                    text += f'\n/approve {code} – ' + ('Plan starten im Standardmodus' if kind == 'plan' else 'einmal erlauben')
                else:
                    text += '\nVollständige Aktion am Dashboard prüfen; Handyfreigabe nicht verfügbar.'
                text += f'\n/deny {code} – ablehnen'
            lines.append(text)
        return '\n\n'.join(lines) or 'Nichts offen.'

    async def refresh(self):
        """Bounded summaries and pending requests only; never fetch a transcript."""
        bindings = store.rooms()
        forums = [chat for chat, r in bindings.items() if r['session_id'] is None and r.get('ready')]
        for chat, room in bindings.items():
            sid = room['session_id']
            if not sid or not room.get('ready') or not self.paired():
                continue
            status = await self.session(sid, 'status')
            next_shown = dict(self._shown)
            notices = []
            if status.get('error'):
                fingerprint = status['error']
                if self._shown.get((chat, 'status')) != fingerprint:
                    notices.append(fingerprint)
                    next_shown[(chat, 'status')] = fingerprint
            else:
                state, answer = status.get('state'), status.get('answer') or ''
                fingerprint = (state, answer)
                if state in {'idle', 'done', 'failed', 'error', 'ended'} and self._shown.get((chat, 'status')) != fingerprint:
                    notices.append(str(state) + ('\n' + answer if answer else ''))
                next_shown[(chat, 'status')] = fingerprint
                pending = await self.session(sid, 'pending')
                items = pending.get('items', [])
                if items and self._shown.get((chat, 'pending')) != items:
                    notices.append(self.pending_text(items))
                next_shown[(chat, 'pending')] = items
            for notice in notices:
                text = f'[{label(room)}]\n{notice}'
                for destination in [chat, *forums]:
                    await self.send(destination, text)
            self._shown = next_shown

    async def run(self):
        await self.daemon.subscribe()
        async def notifications():
            while self.paired():
                try:
                    await self.refresh()
                except Exception:
                    pass  # Retry summaries after a transient service/daemon failure.
                await asyncio.sleep(5)
        task = asyncio.create_task(notifications())
        try:
            async for event in self.daemon.events():
                for message in simplex.messages_of(event):
                    await self.handle(message)
        finally:
            task.cancel()
            await asyncio.gather(task, return_exceptions=True)
