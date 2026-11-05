"""What the phone is shown: compact, clean texts in SimpleX markdown.

Input: the bridge's items (questions, approvals) and the gateway's state.
Output: the message text, in the person's language (lang.t), with SimpleX
markdown (``*bold*``, ``_italic_``, `` `code` ``; see docs/design/
messenger-simplex.md for the official syntax). No raw ids, hashes or
internal tags are rendered.
"""
from __future__ import annotations

from delfin.agent.approval_answers import answer_hint, answer_range

from .lang import english, t


def hint(n: int, multi: bool) -> str:
    """The answer line of a question of ``n`` options, in the phone's language."""
    if english():
        return answer_hint(n, multi)
    numbers = answer_range(n, joiner=t('oder'))
    if multi:
        return t('Antworte mit einer oder mehreren von {range}, durch Kommas getrennt '
                 '(z. B. 1,3), oder mit eigenem Text.').format(range=numbers)
    return t('Antworte mit {range} oder mit eigenem Text.').format(range=numbers)


def question_card(item: dict) -> str:
    """A question: the question in bold, the options numbered, one hint."""
    choices = [str(c) for c in item.get('choices') or []]
    question = str(item.get('question') or item.get('text') or '').strip()
    if choices and not item.get('question'):
        # The bridge's preview carries the options after the question.
        lines = question.splitlines()
        while lines and lines[-1].split(':')[0].strip() in choices:
            lines.pop()
        question = '\n'.join(lines).strip()
    out = [f'❓ *{question}*' if question else '❓']
    out += [f'{i}. {c}' for i, c in enumerate(choices, 1)]
    out.append(hint(len(choices), bool(item.get('multiSelect'))) if choices
               else t('Bitte mit deiner Antwort antworten.'))
    return '\n'.join(out)


#: The tool of an approval, in words (the card's first line).
_WHAT = {
    'bash': 'einen Befehl ausführen', 'write_file': 'eine Datei schreiben',
    'edit_file': 'eine Datei ändern', 'multi_edit': 'eine Datei ändern',
    'apply_patch': 'einen Patch anwenden', 'read_file': 'eine Datei lesen',
    'notebook_edit': 'ein Notebook ändern',
}
_WHAT_BY_ACTION = {'shell': 'einen Befehl ausführen', 'write': 'eine Datei ändern',
                   'read': 'eine Datei außerhalb des Arbeitsordners lesen'}

#: A preview of at most this many lines is shown on the card, in a code
#: block; a longer one on ``show``.
SHORT_PREVIEW = 6
#: A final answer longer than this is cut, with a ``more`` hint.
ANSWER_CUT = 900


def _code(text: str) -> str:
    return '```\n' + text.replace('```', "'''") + '\n```'


def approval_card(item: dict, machine: str = '') -> str:
    """An approval (or a plan) as one compact card: what (a command names
    the ``machine`` it runs on), why, the action when it is short, and how
    to answer. No ids, no rule texts."""
    preview = str(item.get('text') or '').strip('\n')
    lines = preview.splitlines()
    tool = str(item.get('tool') or '')
    if item.get('kind') == 'plan':
        what = 'einen Plan freigeben'
    else:
        what = _WHAT.get(tool) or _WHAT_BY_ACTION.get(item.get('action') or '') or ''
    head = '🔐 *' + (t(what) if what else f'`{tool or "?"}`') + '*'
    if machine and item.get('action') == 'shell':
        head += t(' auf {machine}').format(machine=machine)
    if lines and not preview.startswith('$') and len(lines) > SHORT_PREVIEW:
        head += ': ' + lines[0][:80] + t(' ({n} Zeilen)').format(n=len(lines) - 1)
    out = [head]
    if item.get('purpose'):
        out.append(t('Grund laut Agent (ungeprüft): {purpose}').format(purpose=item['purpose']))
    if lines and len(lines) <= SHORT_PREVIEW:
        out.append(_code(preview))
    elif lines:
        out.append(t('Ganze Aktion: mit `show` antworten.'))
    extra = [o for o in item.get('options') or [] if o.get('key') not in {'y', 'n', 'd', 'e', '?', 'c'}]
    if any(o.get('key') == 'c' for o in item.get('options') or []):
        out.append(t('\nBefehl ändern: mit „ändern: NEUER BEFEHL“ antworten (getippt). '
                     'Dein Befehl läuft dann einmal, genau so.').strip())
    if extra:
        out.append(t('Weitere Antworten: ') + ' · '.join(
            f"`{o['key']}` {o.get('label', '')}".strip() for o in extra))
    if not item.get('phoneApproval'):
        out.append(t('Nur am Rechner freizugeben.'))
    else:
        out.append(t('Antworte mit *ja* oder *nein*.'))
    return '\n'.join(out)


def preview_block(item: dict) -> str:
    """The whole action of an approval, on ``show``."""
    return _code(str(item.get('text') or '').strip('\n') or '-')


def cut(text: str, limit: int = ANSWER_CUT) -> tuple[str, str]:
    """``text`` cut at a paragraph, sentence or word end at most ``limit``
    characters in, and the rest ('' when it fits)."""
    text = str(text or '').strip()
    if len(text) <= limit:
        return text, ''
    head = text[:limit]
    for mark in ('\n\n', '\n', '. ', ' '):
        at = head.rfind(mark)
        if at >= limit // 2:
            head = head[:at + (1 if mark == '. ' else 0)]
            break
    return head.rstrip(), text[len(head):].strip()


def more_hint() -> str:
    return '\n\n… ' + t('ganze Antwort: `more`')


#: A session's state, in words (status).
_STATE = {'idle': 'Bereit.', 'running': '⏳ Arbeitet …', 'needs_input': 'Wartet auf eine Antwort.',
          'ended': 'Sitzung beendet.'}


def status_text(reply: dict) -> str:
    """``status``: the state in words, what it is doing, the last answer
    (cut)."""
    state = str(reply.get('state') or '')
    out = [t(_STATE[state]) if state in _STATE else state]
    if reply.get('activity'):
        out.append(str(reply['activity']))
    if reply.get('answer'):
        out.append(cut(str(reply['answer']), 400)[0])
    return '\n'.join(filter(None, out))


def in_forum(title: str, body: str) -> str:
    """A session's post in the forum: its name in bold, then the text."""
    return f'*{title}*\n{body}'
