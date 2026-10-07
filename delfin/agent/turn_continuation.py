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
    r"|\bi(?:'ll| will|'m going to| am going to)\b"
    # A bare participle closing the turn: "Pushing the branch now.",
    # "Running the tests now." No first person, no future auxiliary, and
    # the same meaning. Report 20261005-140408 ended one of its four
    # turns exactly so, and that turn was the one the detector missed.
    #
    # Two things keep this narrow. It must END the text (optionally
    # followed by closing punctuation), so a participle in the middle of
    # a report -- "Pushing it would need a grant, so I stopped" -- does
    # not match. And it requires the word "now", which is what makes the
    # sentence a statement of what is happening rather than a heading or
    # a caption ("Pushing the branch" as a task title).
    r"|\b\w+ing\b[^.!?\n]{0,80}\bnow\b[\s.!?\u2026]*$",
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
        "will follow up later",
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

    Review finding (nacht-s12): real transcripts use the U+2019 right
    single quote, so "I'll…" and "I'm…" fell through the ASCII-only class.
    The text is normalised to ASCII apostrophes BEFORE any scan, so the
    signal and the hand-back table both see one shape.
    """
    if not text or not text.strip():
        return False
    text = (str(text).replace("\u2019", "'")
            .replace("\u2018", "'"))
    for rx in _HANDBACK_RES:
        if rx.search(text):
            return False
    return bool(_RE_SIGNAL.search(text))


def should_continue(answer: str, open_tasks: int) -> tuple[bool, str]:
    """Whether the session should be nudged to continue, and the note.

    One condition: the answer announces agent work it did not do.
    ``open_tasks`` is accepted and reported in the note, but no longer
    required.

    It used to be required, and that is what reports 20261005-135825 and
    20261005-140408 are. The user typed "push den branch", then "push",
    then "push den brnach", then "push"; each turn answered "Let me push
    the branch now" and called nothing. The task list read
    ``{"completed": 5}``, so ``open_tasks`` was 0 and no follow-up fired
    -- the detector asked whether work was left ON THE LIST when the
    question is whether THIS TURN announced an action it did not take.
    An announcement is the signal; a task row is at most a hint about
    where to resume.

    The hand-back table above is what keeps this from nagging, and it
    carries the whole weight now that the task count does not: a turn
    that gives the work back ("Let me know...", "I'll keep you posted")
    announces nothing and is not nudged, with or without open tasks.
    The engine's own latch (``_continuation_fired``) still allows at most
    one note per turn, so dropping the count cannot loop.
    """
    if not announced_action(answer):
        return False, ""
    n = max(0, open_tasks or 0)
    if n:
        return True, _CONTINUATION_NOTE
    return True, _CONTINUATION_NOTE_NO_TASKS


#: The note for an announcement with an empty task list. It must not talk
#: about tasks: a model told to "continue the announced work" while the
#: list says everything is done has been given two contradictory facts,
#: and the recorded sessions spent their turns re-reading the list.
_CONTINUATION_NOTE_NO_TASKS = (
    "System note: this turn ended by announcing an action (\"Let me \u2026\", "
    "\"I will \u2026\") and did not take it. Take it now, in this turn. If it "
    "cannot be taken -- a gate refused it, something is missing -- say that "
    "to the user in one sentence, naming what blocked it and what they can "
    "do, instead of announcing it again."
)


def blocked_by_plan_mode(answer: str, *, refusals: int) -> str:
    """One sentence for the USER when plan mode refused the announced act.

    Input: the final answer, and how many calls plan mode refused this
    turn. Output: the sentence to append to what the user reads, or "".

    Report 20261005-140408. The session was in plan mode, `bash` was
    refused read-only, and that refusal is addressed to the model --
    "call exit_plan_mode with the plan". The person read "Let me push the
    branch now" four times and saw no push; the model itself worked out
    the cause one turn later and still did not say it.

    Nudging the model is the wrong answer here and that is why this is
    separate from the note above: in a read-only session the next turn
    cannot act either, so repeating the attempt only spends turns. The
    person is the only one who can lift it.

    Silent unless BOTH hold: plan mode actually refused something, and
    the turn announced an action. A read-only investigation in plan mode
    is the normal case and must not carry a warning.
    """
    if refusals <= 0 or not announced_action(answer):
        return ""
    return (
        "Note: this session is in plan mode, which is read-only, and "
        f"{refusals} call{'s' if refusals != 1 else ''} were refused "
        "because of it \u2014 so the action named above did not run. "
        "Nothing will run until the plan is approved (exit_plan_mode) or "
        "you leave plan mode (`/mode solo`, or Perms in the dashboard)."
    )
