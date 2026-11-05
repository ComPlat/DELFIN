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
