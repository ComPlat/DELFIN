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
        "waiting on",
        "waiting for",
        "ich warte",
        "warte auf",
        "warte ich",
        # German: the turn addresses the user -- a request, an
        # imperative, a question, an offer, or the user's approval as
        # the precondition. Each of these occurred in an archived
        # closing that also carried a first-person verb ("Sag kurz Go,
        # dann pushe ich ...", "Danach prüfe ich gern ...").
        "mich wissen",
        "bescheid",
        "sag",
        "sag mir",
        "sage mir",
        "sagen",
        "bestätige",
        "bestätigen",
        "bitte",
        "gern",
        "gerne",
        "soll ich",
        "sollen wir",
        "möchtest du",
        "möchten sie",
        "willst du",
        "wollen sie",
        "melde mich",
        "mich melden",
        "warten",
        "rückmeldung",
        "freigabe",
        "freigibst",
        "freigeben",
        "akzeptierst",
        "akzeptieren",
        "zustimmung",
        "grünes licht",
        "sofern",
        "vorausgesetzt",
    )
)

#: German first person + an immediate marker, or the future, or the
#: idioms that mean "I'll take this on". Present tense with "jetzt" is
#: how German announces ("Ich prüfe jetzt die Logs."), so the verb is
#: any first-person form ending in -e except the verbs that state a
#: belief, a perception, a capacity or a stop ("Ich sehe jetzt, dass
#: ...", "Ich habe jetzt ..."). The verb must be lower case, so a
#: capitalised noun after "Ich" is never taken for one.
#:
#: Rejected: a closed list of action verbs -- the archive's own
#: announcements used "verschaffe", "extrahiere", "beginne", none of
#: which a list written beforehand held; the excluded verbs are the
#: smaller and more stable set. Rejected: the bare present without a
#: marker ("Ich prüfe die Logs.") -- archived closings use it for what
#: the session already does or keeps doing ("Ich ignoriere das Ereignis
#: weiter", "ich markiere sie als unverifiziert").
_DE_NOT_ACTION = (
    r"habe|hatte|bin|sehe|glaube|denke|meine|finde|vermute|hoffe|"
    r"weiß|kenne|verstehe|erkenne|merke|brauche|warte|bleibe|lasse|"
    r"halte|gebe|übergebe|sage|frage|bitte|danke|empfehle|schlage|"
    r"rate|höre|beende|stoppe|pausiere|ignoriere|behalte|gehe|"
    r"fasse|nehme|stelle|melde|möchte|würde|könnte|sollte|müsste|"
    r"freue|stimme"
)
_DE_VERB = rf"(?!(?:{_DE_NOT_ACTION})\b)[a-zäöüß]+e"
_DE_NOW = (
    r"(?i:jetzt|nun|gleich|sofort|direkt|zuerst|zunächst|erst einmal"
    r"|erstmal|als nächstes|als erstes|kurz|noch|mal)"
)
_RE_SIGNAL_DE = re.compile(
    # "Ich prüfe jetzt die Logs.", "Ich analysiere zuerst, was ..."
    rf"\b(?i:ich)\s+{_DE_VERB}\b[^.!?,;:\n]{{0,60}}?\b{_DE_NOW}\b"
    # "Als Nächstes starte ich ...", "Zuerst verschaffe ich mir ..."
    rf"|\b(?i:jetzt|nun|gleich|sofort|zuerst|zunächst|erst einmal"
    rf"|erstmal|als nächstes|als erstes)\s+{_DE_VERB}\s+(?i:ich)\b"
    # "Ich beginne mit der Implementierung:" -- the colon promises what
    # follows, and nothing followed.
    rf"|\b(?i:ich)\s+{_DE_VERB}\b[^.!?\n]{{0,80}}:\s*$"
    r"|\b(?i:ich werde|werde ich)\b"
    r"|\b(?i:lass|lasst|lassen sie) mich\b"
    r"|\b(?i:ich (?:schaue|sehe|gucke) mir)\b[^.!?,;\n]{0,60}\ban\b"
    r"|\b(?i:ich kümmere mich)\b"
    r"|\b(?i:ich mache mich)\s+(?:\w+\s+)?(?:an|daran|dran)\b"
    r"|\b(?i:ich lege)\s+(?:\w+\s+)?los\b",
)

#: A German signal is void when its clause is negated ("Ich werde nichts
#: pushen") or the sentence is a question or a rule for later ("Ich werde
#: das künftig beachten").
_RE_DE_NEGATION = re.compile(
    r"\b(?:nicht|nichts|kein\w*|nie|niemals)\b", re.IGNORECASE)
_RE_DE_LATER = re.compile(
    r"\b(?:künftig|zukünftig|in zukunft|von nun an|ab sofort|"
    r"beim nächsten mal|nächstes mal|merken|beachten|daran denken|"
    r"erinnern)\b", re.IGNORECASE)
_RE_CLAUSE_END = re.compile(r"[,;:—–]")

#: A future made conditional on someone else's step ("once C names
#: candidates", "as soon as the user approves", "sobald die Freigabe da
#: ist") hands the turn back: the action waits for an input this session
#: cannot produce. Judged per sentence -- only the sentence that carries
#: the signal is cancelled. German markers are matched because users and
#: models write German.
#:
#: Kept narrow where the bare word is ordinary action vocabulary: "once"
#: excludes "once more" / "once again" / "at once"; "when", "after",
#: "until", "if" and "unless" count only before a subject other than the
#: agent ("when you", "after the user", "if you"), so "Let me check if the
#: file exists" and "Let me check when the job ended" still announce.
_RE_CONDITIONAL = re.compile(
    r"(?<!at )\bonce\b(?!\s+(?:more|again)\b)"
    r"|\bas soon as\b|\bwhenever\b|\bprovided (?:that )?\b"
    r"|\b(?:when|after|until|if|unless) (?:you|your|the user|they|he|she"
    r"|session \w+|approval|approved|confirmed)\b"
    r"|\bon your (?:go|signal|approval|confirmation)\b"
    r"|\bsobald\b|\bwenn\b|\bfalls\b|\bnachdem\b|\bbis (?:du|sie|ihr)\b",
    re.IGNORECASE,
)

#: "when C names ...", "after B pushes": a capitalised party other than
#: "I". Case-sensitive, so a separate pattern.
_RE_CONDITIONAL_NAMED = re.compile(
    r"\b(?:[Ww]hen|[Aa]fter|[Uu]ntil) (?!I\b)[A-Z]\w*\b")

#: Sentence boundary: closing punctuation followed by whitespace, or by an
#: upper-case letter with no space between. The second form is the seam
#: left where the text of consecutive rounds is joined ("...first.The
#: file exists").
_RE_SENTENCE_END = re.compile(
    r"(?<=[.!?…])(?:\s+|(?=[A-ZÄÖÜ]))")

#: How many closing sentences of the last paragraph are judged.
CLOSING_SENTENCES = 2

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
    "Next I'll…", "Ich prüfe jetzt…", "Lass mich…"). False for an empty/whitespace text, a hand-back to the
    user ("Let me know…", "What should I do?"), waiting, or completion.

    Review finding (nacht-s12): real transcripts use the U+2019 right
    single quote, so "I'll…" and "I'm…" fell through the ASCII-only class.
    The text is normalised to ASCII apostrophes BEFORE any scan, so the
    signal and the hand-back table both see one shape.

    Only the CLOSING part is judged: the last paragraph, and within it the
    last ``CLOSING_SENTENCES`` sentences. A hand-back phrase or a
    conditional / waiting future ("once C names candidates", "sobald ...")
    anywhere in that part closes the turn -- across the sentence boundary,
    because the action sentence often refers back to the condition ("Ich
    prüfe das, sobald C nennt. I'll check it then.").

    Field report 20261012-082746: judged on the whole answer, 9 of the 10
    answers that drew a note matched, 8 of them falsely -- the signal came
    from narration written before a tool call ("Let me first verify ...")
    or from an announcement already followed by its result, both of which
    end up in the turn's answer. Judged on the closing part, 1 of 10
    matches, the real one. Rejected: the very last sentence only -- "Let me
    run the suite. This takes a minute." announces, and its last sentence
    carries no signal.
    """
    closing = _judged_sentences(text, True)
    if _is_conditional(" ".join(closing)):
        return False
    return any(_announces(s) for s in closing)


def _judged_sentences(text: str, closing_only: bool) -> list[str]:
    """Sentences to judge, apostrophes normalised.

    *closing_only*: the last ``CLOSING_SENTENCES`` sentences of the last
    non-empty paragraph; otherwise every sentence of the text. Returns []
    for an empty text, or when a hand-back phrase sits in the judged part.
    """
    text = str(text or "").replace("\u2019", "'").replace("\u2018", "'")
    if not text.strip():
        return []
    if closing_only:
        text = [p for p in re.split(r"\n\s*\n", text) if p.strip()][-1]
    sentences = [s for s in _RE_SENTENCE_END.split(text.strip()) if s.strip()]
    if closing_only:
        sentences = sentences[-CLOSING_SENTENCES:]
    judged = " ".join(sentences)
    if any(rx.search(judged) for rx in _HANDBACK_RES):
        return []
    return sentences


def _announces(sentence: str) -> bool:
    """A first-person action signal in *sentence* that no conditional or
    waiting future qualifies."""
    sentence = sentence.strip()
    if not (_RE_SIGNAL.search(sentence) or _announces_de(sentence)):
        return False
    return not _is_conditional(sentence)


def _announces_de(sentence: str) -> bool:
    """A German first-person signal in *sentence*, its clause not negated,
    the sentence neither a question nor a rule for later turns.

    The negation is read in the clause the signal opens (up to the next
    comma, colon, semicolon or dash), so "Ich prüfe jetzt, warum das nicht
    läuft" still announces while "Ich werde das nicht pushen" does not.
    """
    if sentence.endswith("?") or _RE_DE_LATER.search(sentence):
        return False
    for m in _RE_SIGNAL_DE.finditer(sentence):
        rest = sentence[m.start():]
        end = _RE_CLAUSE_END.search(rest, m.end() - m.start())
        clause = rest[:end.start()] if end else rest
        if not _RE_DE_NEGATION.search(clause):
            return True
    return False


def _is_conditional(text: str) -> bool:
    return bool(_RE_CONDITIONAL.search(text)
                or _RE_CONDITIONAL_NAMED.search(text))


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
    # A turn that registered a WAIT is not nudged: the wait is the
    # mechanism that brings the session back when its event arrives
    # (wait_registry.due fires it once, with a note naming the event).
    # Nudging on top would have the model re-announce work it is
    # deliberately parked on (package KG1, phase 2). An empty or
    # unreadable registry nudges as before.
    try:
        from .wait_registry import SessionWaits
        if SessionWaits().list():
            return False, _WAIT_NOTE
    except Exception:
        pass
    n = max(0, open_tasks or 0)
    if n:
        return True, _CONTINUATION_NOTE
    return True, _CONTINUATION_NOTE_NO_TASKS


#: The note for a suppressed nudge. It names the mechanism, so the person
#: reading the status line sees WHY the session is not being poked: it is
#: waiting, on purpose, and the wait will wake it.
_WAIT_NOTE = (
    "System note: the nudge is suppressed because this turn registered a "
    "wait (wait_registry); the session wakes when the waited-for event "
    "arrives. The person can clear the wait."
)


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

    Unlike ``announced_action`` this reads every sentence of the turn's
    whole answer: the refused call comes AFTER its announcement, and the
    refusal is itself the evidence that the announced act did not run.
    """
    if refusals <= 0 or not any(
            _announces(s) for s in _judged_sentences(answer, False)):
        return ""
    return (
        "Note: this session is in plan mode, which is read-only, and "
        f"{refusals} call{'s' if refusals != 1 else ''} were refused "
        "because of it \u2014 so the action named above did not run. "
        "Nothing will run until the plan is approved (exit_plan_mode) or "
        "you leave plan mode (`/mode solo`, or Perms in the dashboard)."
    )
