"""A go-ahead the turn did not act on, and a choice asked in prose.

Input: the user's message that opened the turn, the model's final answer
and the names of the tools the turn called. Output:

* :func:`is_go_ahead` -- the message tells the agent to proceed ("ja",
  "ja mach", "setz es um", "go ahead", "warum machst du nicht einfach?");
* :func:`ignored_go_ahead` -- the turn then ended with no change tool
  called at all: it restated a plan, or asked again;
* :func:`is_choice_question` -- the answer ends by asking the user to pick
  between options written out as text.

Measured 2026-10-19 (session 78c36f8c, gpt-5-nano building a browser game):
after 4 of 7 go-aheads the turn called no change tool and answered with
the plan again ("Was ich jetzt mache ..."); the user wrote "warum machst
du nicht einfach??" and "du sollst alles angehen und umsetzen". The
choices -- width 10, 12 or 16 -- came as numbered prose ("bitte eine Zahl
wählen") although ask_user_question was bound; the user answered "sowas
kannst du doch als rückfrage machen". Neither was caught: an answer that
ends in a question counts as handed back, so the announced-action check
stands down on exactly these turns.

How a good agent handles both: a clear instruction is carried out, not
re-confirmed; only a decision the user actually has to make is asked, as
a structured question with the options and a recommendation, and a choice
with a sensible default is made and named instead of asked.

Pure; never raises. German and English on the recognising side.
"""

from __future__ import annotations

import re

#: Tools that change something or run work: calling any of them in the
#: turn is acting on the go-ahead, whatever came of it (a refused write is
#: answered by the change-claim caveat, not here).
ACTING_TOOLS = frozenset({
    "write_file", "edit_file", "multi_edit", "apply_patch", "notebook_edit",
    "bash", "bash_background", "run_tests", "edit_sheet", "task_create",
    "task_update", "subagent", "undo_changes",
})

# A short message that is a yes or an order to proceed. Anchored to the
# start: "ja, aber erst ..." is still a yes to the start of the work; a
# question about state ("was ist der stand?") is not.
_GO_START_RE = re.compile(
    r"(?i)^\W*(?:ja+|jo|jup|jap|jep|yes|yep|yeah|ok(?:ay)?|okey|passt|gut|"
    r"genau|klar|sicher|einverstanden|go|los|leg\s+los|weiter|mach|"
    r"machs?|mache|tu\s+es|tus?|bitte|do\s+it|go\s+ahead|proceed|"
    r"continue|sure|setz|bau|fix|patch|umsetzen|starte?|implementier)\b")
# An order anywhere in a short message ("du sollst alles umsetzen", "warum
# machst du nicht einfach??", "Patch a + B").
_GO_ANY_RE = re.compile(
    r"(?i)\b(?:mach(?:\s+(?:es|das|mal|einfach|weiter))?\s*(?:$|[.!])|"
    r"umsetz\w*|setz\w*\s+(?:es|das|alles)?\s*um\b|"
    # "nicht" as typed in a hurry: nciht, nich, nicth.
    r"warum\s+machst\s+du\s+(?:es\s+|das\s+)?n(?:icht|ciht|ich|icth)\b|"
    r"\bdu\s+(?:musst|sollst)\b[^?]{0,60}\b(?:fix\w*|mach\w*|bau\w*|"
    r"umsetz\w*|änder\w*|aender\w*|test\w*|start\w*|patch\w*)|"
    r"\b(?:baust|machst|startest|fixt|patchst)\s+du\s+(?:das|es)\s+"
    r"(?:jetzt|jetz|nun|gleich)\b|"
    r"mach\s+(?:es|das)\s+einfach|einfach\s+machen|"
    r"just\s+do\s+it|do\s+it\s+now|(?:yes|ok(?:ay)?|sure),?\s+do\s+it|"
    r"go\s+ahead|"
    r"patch\s+\w+(?:\s*\+\s*\w+)?)")
# What turns a short "ja" into something else: a refusal or a question
# about the plan itself.
_NOT_GO_RE = re.compile(
    r"(?i)\b(?:nein|nicht\s+(?:jetzt|machen)|stop+|halt|warte|wait|"
    r"don'?t|cancel|abbrechen|lass\s+(?:es|das))\b")

_MAX_GO_CHARS = 160


# The agent's previous answer offered or asked something: then a bare "ja"
# or "ok" answers it. After a finished report ("Fertig, Tests grün.") the
# same "ok" is an acknowledgement, and pushing the agent to act on it would
# invent work.
_OFFERED_RE = re.compile(
    r"(?i)\?|\b(?:soll\s+ich|sollen\s+wir|shall\s+i|should\s+i|want\s+me\s+to|"
    r"do\s+you\s+want|möchtest\s+du|willst\s+du|bestätig\w*|confirm|"
    r"vorschlag|proposal|plan\b|nächste[rn]?\s+schritt|next\s+step)")


def is_go_ahead(user_text: str, offered: str = "") -> bool:
    """Whether *user_text* tells the agent to go ahead (module docstring).

    An explicit order ("mach", "setz es um", "warum machst du nicht
    einfach") always counts. A bare yes ("ja", "ok", "go") counts only as
    the answer to *offered* -- the agent's previous answer -- when that
    answer asked or proposed something (``_OFFERED_RE`` on its last
    paragraph)."""
    t = " ".join(str(user_text or "").split())
    if not t or len(t) > _MAX_GO_CHARS or _NOT_GO_RE.search(t):
        return False
    if _GO_ANY_RE.search(t):
        return True
    # "ja?" / "ok?" is a question back, not a yes.
    if not _GO_START_RE.search(t) or t.rstrip().endswith("?"):
        return False
    paras = [x for x in re.split(r"\n\s*\n", str(offered or "")) if x.strip()]
    return bool(paras) and bool(_OFFERED_RE.search(paras[-1]))


def ignored_go_ahead(user_text: str, tools_called, offered: str = "") -> bool:
    """A go-ahead after which the turn called no acting tool."""
    if not is_go_ahead(user_text, offered):
        return False
    names = {str(n or "").split("__")[-1] for n in (tools_called or ())}
    return not (names & ACTING_TOOLS)


# Options written out: two or more list items (-, *, •, 1., 1), a), A))
# or an "A / B / C" or "oder" choice, followed by a question or a request
# to choose.
_OPTION_LINE_RE = re.compile(
    r"(?m)^\s*(?:[-*•]|\d{1,2}[.)]|[A-Fa-f][.)])\s+\S")
_INLINE_CHOICE_RE = re.compile(
    r"(?i)\b(?:option|variante)\s+[A-F1-6]\b.*\b(?:option|variante)\s+[A-F1-6]\b"
    r"|\S+\s+(?:oder|or)\s+\S+\s*\?")
_ASKS_TO_PICK_RE = re.compile(
    r"(?i)\?|\b(?:wähle|waehle|wählen|waehlen|aussuchen|entscheide|"
    r"sag\s+mir|gib\s+mir\s+bescheid|choose|pick|which\s+(?:one|option)|"
    r"let\s+me\s+know\s+which)\b")


def is_choice_question(answer: str) -> bool:
    """Whether *answer* ends by asking the user to choose between options
    written out as text.

    Only the closing part is judged -- the last two paragraphs -- so a
    report with a list in its middle that ends "Fertig." is not a choice.
    """
    body = re.sub(r"```.*?```", " ", str(answer or ""), flags=re.S).strip()
    if not body:
        return False
    paras = [p for p in re.split(r"\n\s*\n", body) if p.strip()]
    tail = "\n\n".join(paras[-2:])
    last = paras[-1]
    if not _ASKS_TO_PICK_RE.search(last):
        return False
    return (len(_OPTION_LINE_RE.findall(tail)) >= 2
            or bool(_INLINE_CHOICE_RE.search(last)))


GO_AHEAD_NOTE = (
    "The user already told you to go ahead. Do not restate the plan and do "
    "not ask for confirmation again: make the change now with the tools, "
    "then check that it works (run it, test it, read it back). Keep going "
    "until the request is done. If a decision the user has not made truly "
    "blocks you, ask it with ask_user_question; if something else blocks "
    "you, say in one sentence what it is and what the user can do."
)

CHOICE_NOTE = (
    "You ended your turn asking the user to choose between options written "
    "as text. If one option is the sensible default the user would accept, "
    "do not ask: pick it, say in one line which you picked and why, and "
    "act. If the choice is genuinely theirs, ask it with ask_user_question "
    "instead: one question, 2-4 options, each with a one-line trade-off, "
    "your recommendation first and marked (Recommended). Do not repeat the "
    "options as prose."
)

CHOICE_NOTE_NO_UI = (
    "You ended your turn asking the user to choose between options written "
    "as text. If one option is the sensible default the user would accept, "
    "do not ask: pick it, say in one line which you picked and why, and "
    "act. If the choice is genuinely the user's, keep the question as it "
    "is -- short, the options and your recommendation."
)


__all__ = ["ACTING_TOOLS", "is_go_ahead", "ignored_go_ahead",
           "is_choice_question", "GO_AHEAD_NOTE", "CHOICE_NOTE",
           "CHOICE_NOTE_NO_UI"]
