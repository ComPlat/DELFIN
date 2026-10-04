"""Turn-continuation detector: announced-but-not-done intent plus open work.

Wave 12, package R1, finding 1: sessions ended turns with "Let me write X …"
while their task list had open items and then sat idle until poked. This
module is the PURE detector — no IO, no engine state, only the last answer's
text and a count of open tasks. The engine wiring that starts at most one
follow-up turn lives as a patch for the operator (engine.py is protected);
this module only answers "should the session be nudged to continue?".

Design is deliberately conservative: the costly failure is a FALSE POSITIVE —
nudging a session whose answer already closed the turn ("Let me know…",
"I'm done…", "I'll wait…"). A false negative (missing a real announcement)
only reproduces today's behaviour: the session stays idle, nothing breaks.
So hand-backs to the user are excluded first; a genuine first-person intent
to keep acting in this session is then required.
"""

from __future__ import annotations

import re

#: The engine starts at most ONE follow-up turn per wake cycle. Stated here
#: so the loop guard and its test share one constant; enforcement is the
#: engine's job (this module never loops on its own).
CONTINUATION_LIMIT = 1

#: Signals a future first-person ACTION in this session. `let me` is the
#: strongest ("let me write…"); the i-will / i'll / i'm-going-to forms are the
#: explicit future. Anchored with ``\\b`` on the leading ``i`` so "ill",
#: "Iwille" or "the suite will" never match. Case-insensitive.
_RE_SIGNAL = re.compile(
    r"\blet me\b"
    r"|\bi(?:'ll| will|'m going to| am going to)\b",
    re.IGNORECASE,
)

#: Phrases that say the turn is being CLOSED or handed back to the user,
#: even though they contain a signal word ("let me know", "i'll be here").
#: Presence of any of these suppresses the continuation regardless of other
#: words: a closed turn must not be reopened by the detector.
_HANDBACK_RES = tuple(
    re.compile(rf"\b{re.escape(phrase)}\b", re.IGNORECASE)
    for phrase in (
        "let me know",
        "let me ask",
        "let me hear",
        "let me get back",
        "i will wait",
        "i'll wait",
        "i will be here",
        "i'll be here",
        "i will stop",
        "i'll stop",
        "i'm done",
        "i am done",
        "i will keep you posted",
        "i'll keep you posted",
        "i will let you",
        "i'll let you know",
        "i will get back",
        "i'll get back",
        "waiting for your",
    )
)

#: The short system note injected into the at-most-one follow-up turn.
_CONTINUATION_NOTE = (
    "System note: this turn ended announcing further work (\"Let me …\", "
    "\"I will …\") but the task list is not empty. Continue the announced "
    "work now. Do not end another turn purely by announcing what you will "
    "do; either do it, or if it truly cannot proceed mark the tasks blocked "
    "and say what is missing."
)


def announced_action(text: str) -> bool:
    """Whether the final answer announces an agent ACTION in this session.

    True for a first-person future intent ("Let me fix…", "I will run…",
    "Next I'll…"). False for an empty/whitespace text, a hand-back to the
    user ("Let me know…", "What should I do?"), waiting, or completion.
    """
    if not text or not text.strip():
        return False
    for rx in _HANDBACK_RES:
        if rx.search(text):
            return False
    return bool(_RE_SIGNAL.search(text))


def should_continue(answer: str, open_tasks: int) -> tuple[bool, str]:
    """Whether the session should be nudged to continue, and the note.

    Both conditions must hold: the answer announces more agent work AND open
    tasks remain (``open_tasks`` > 0; a negative or zero count never nudges).
    Returns ``(True, note)`` when both, else ``(False, "")``.
    """
    n = max(0, open_tasks or 0)
    if n == 0 or not announced_action(answer):
        return False, ""
    return True, _CONTINUATION_NOTE
