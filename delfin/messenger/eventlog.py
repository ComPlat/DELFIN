"""One content-free line per message the gateway handles, for the journal.

Input: what the gateway knows about a message once it handled it. Output,
on stderr (the unit's journal), one line::

    [gateway] event in kind=prompt person=anna room=r3f2a session=s7c1e result=delivered ms=42

* direction ``in`` (from a phone) or ``out`` (to a chat);
* ``kind``: prompt, a command's name, approval_answer, question_answer,
  file, voice; out: text, card, question;
* the person (``owner`` for the account owner, ``-`` when unknown);
* the room and the session as short tags (four hex digits of a hash),
  never their ids;
* ``result``: delivered, refused:REASON, no_session, starting, failed:TYPE;
* ``ms``: the time the gateway took.

Never a message text, link, code, file name or title: only the fields
above, each from a fixed alphabet. On by default; DELFIN_MESSENGER_EVENTS=off
switches it off. The same line within REPEAT_S seconds is counted, not
printed, and the next one says how many were left out (``repeated=N``).

Rejected: logging the first characters of a message (a prompt can hold a
secret; field debugging needs the path a message took, not its words).
"""
from __future__ import annotations

import hashlib
import os
import re
import sys
import time

#: Identical lines within this many seconds are collapsed.
REPEAT_S = 30.0
_WORD = re.compile(r'[a-z_]{1,24}')
_RESULT = re.compile(r'[a-z_]{1,24}(?::[A-Za-z_]{1,40})?')
_last: dict[tuple, list] = {}


def enabled() -> bool:
    return os.environ.get('DELFIN_MESSENGER_EVENTS', '').strip().lower() not in {'off', '0', 'no', 'false'}


def tag(prefix: str, value) -> str:
    """A short tag for a room or session id: never the id itself."""
    if value is None or value == '':
        return '-'
    return prefix + hashlib.sha256(str(value).encode()).hexdigest()[:4]


def _person(person) -> str:
    if person is None:
        return '-'
    if person == '':
        return 'owner'
    name = re.sub(r'[^A-Za-z0-9_-]', '', str(person))[:32]
    return name or '-'


def emit(direction: str, kind: str, *, person=None, room=None, session=None,
         result: str = 'delivered', ms: float | None = None) -> None:
    """Print the line for one message (see the module text). Never raises."""
    try:
        if not enabled():
            return
        direction = 'out' if direction == 'out' else 'in'
        kind = kind if _WORD.fullmatch(str(kind or '')) else 'other'
        result = result if _RESULT.fullmatch(str(result or '')) else 'other'
        fields = (direction, kind, _person(person), tag('r', room), tag('s', session), result)
        now = time.monotonic()
        seen = _last.get(fields)
        if seen is not None and now - seen[0] < REPEAT_S:
            seen[1] += 1
            return
        left_out = seen[1] if seen is not None else 0
        _last[fields] = [now, 0]
        while len(_last) > 256:
            _last.pop(next(iter(_last)))
        line = ('[gateway] event {} kind={} person={} room={} session={} result={}'.format(*fields)
                + (f' ms={int(ms)}' if ms is not None else '')
                + (f' repeated={left_out}' if left_out else ''))
        print(line, file=sys.stderr, flush=True)
    except Exception:  # noqa: BLE001 - a log line never costs a message
        pass
