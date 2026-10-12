"""skill_feedback: does the user's next message correct a skill used this turn?

Packet L2, the learning loop: skill_patch can propose a revision, but
nothing connects "a skill was used AND the user corrected the result" to
a patch proposal. This module is the pure decision half of that link. It
turns

    (the skill used this turn, the user's next message)

into a verdict: does that message correct the skill, and if so what are
the user's own words (which skill_patch will turn into a proposed patch —
always a proposal, never applied by itself).

Security (waves 14/15; this is a memory-poisoning surface):

* What may be learned comes from the USER's own messages only — never a
  tool result, a job log, a file, a web page or another session's
  message. ``feedback`` therefore only recognises a correction when the
  caller marks the message as ``source="user"`` (the user's own next
  message on THIS session). Any other source is data, not instruction,
  and proposes nothing.
* The correction text is the user's own words, verbatim — this module
  passes the message through with at most surrounding-whitespace
  trimming. It never fabricates, paraphrases or rewrites the user.
* A proposed patch must pass the live safety gate: ``patched_skill_safe``
  delegates to ``skill_safety.check`` (the same check skill_proposals
  runs inside ``propose``), so a patched skill that adds a deny-pattern
  command or a command that is not auto-allowable is refused.
* Nothing here activates a skill; only a human accepting a proposal does.

The markers are ``skill_learning.CORRECTION_MARKERS``, the list the
codebase already uses for "this user message reads as a course
correction". Two things are narrower here than there, because here a
match puts the user's words on record against a skill:

* only the user's own prose is read: fenced and inline code, quoted
  spans ("...", '...', „...“, «...»), ``>`` quote lines and pasted
  output (indented lines, tracebacks, ``...Error:`` lines) are dropped
  first. A pasted ``ValueError: wrong number of atoms`` or a quoted
  reviewer sentence is not the user saying the skill was wrong.
* a marker is a whole word: "anders" does not match "anderson"; a German
  inflection ending (falsch -> falsche, falschen) and "-ly" still match.

Measured on 17 realistic lines (German and English): 6 false positives
and 2 false negatives with the substring reading, 2 and 2 with this one.
The two left are semantic ("Instead of a summary, could you also...",
"Kannst du die Ergebnisse anders darstellen?") and need more than a
word list; the steering note they cause tells the model to leave the
skill alone when the correction is about something else.
"""

from __future__ import annotations

import re
from typing import Optional

from .skill_learning import CORRECTION_MARKERS

__all__ = ["feedback", "patched_skill_safe"]

_CODE_RE = re.compile(r"```.*?(?:```|\Z)|`[^`\n]*`", re.DOTALL)
_QUOTED_RE = re.compile(
    r'"[^"\n]{0,400}"'
    r"|(?<![^\W_])'[^'\n]{0,400}'(?![^\W_])"
    r"|„[^“”\n]{0,400}[“”]"
    r"|“[^”\n]{0,400}”"
    r"|‘[^’\n]{0,400}’"
    r"|«[^»\n]{0,400}»"
    r"|»[^«\n]{0,400}«")
_PASTED_LINE_RE = re.compile(
    r"^(?:[ \t]{2,}\S|\t|\s*>"
    r"|\s*Traceback \(most recent call last\)"
    r"|\s*File \"[^\"]*\", line \d"
    r"|\s*[\w.]*(?:Error|Exception|Warning|FATAL|ERROR|WARNING)\b\s*:)")
_MARKER_RE = re.compile("|".join(
    r"(?<![^\W\d_])" + re.escape(m)
    + (r"(?:e[nmrs]?|ly)?(?![^\W\d_])" if m[-1].isalpha() else "")
    for m in sorted(CORRECTION_MARKERS, key=len, reverse=True)))


def _own_words(text: str) -> str:
    """The text with code, quotes and pasted output taken out."""
    text = _QUOTED_RE.sub(" ", _CODE_RE.sub(" ", text))
    return "\n".join(line for line in text.splitlines()
                     if not _PASTED_LINE_RE.match(line))


def _reads_as_correction(text: str) -> bool:
    return bool(_MARKER_RE.search(_own_words(text).lower()))

# A correction is recognised only for a real user message. Tool results,
# files, web pages and other sessions' messages are data, never
# instruction — the caller states which of these the message is.
_USER_SOURCE = "user"


def feedback(
    skill_used: Optional[str],
    user_message: Optional[str],
    *,
    source: str = "unknown",
) -> dict:
    """Verdict: does the user's next message correct the skill used this turn?

    Pure — no I/O, no mutation. Returns a dict the wiring may hand to
    skill_patch:

        {"corrects": bool, "skill": str|None, "user_words": str|None}

    * ``corrects`` is True only when a skill was used, the message is a
      real user message (``source == "user"``) and it reads as a course
      correction.
    * ``skill`` is the used skill's name (trimmed) when ``corrects``.
    * ``user_words`` is the user's own message, trimmed of surrounding
      whitespace only, shared verbatim so skill_patch can build the
      patch from the user's actual words — never a paraphrase.

    A correction arriving inside a tool result or another session's
    message (``source != "user"``) proposes nothing, even when the text
    reads like one. An empty or absent message, or an absent skill,
    proposes nothing.
    """
    skill = skill_used if isinstance(skill_used, str) else ""
    skill = skill.strip()
    text = user_message if isinstance(user_message, str) else ""
    text = text.strip()
    corrects = bool(skill) and source == _USER_SOURCE and bool(text) and (
        _reads_as_correction(text)
    )
    if not corrects:
        return {"corrects": False, "skill": None, "user_words": None}
    return {"corrects": True, "skill": skill, "user_words": text}


def patched_skill_safe(new_text: str) -> tuple[bool, list]:
    """Safety gate for a patched skill: (bool, findinds).

    Delegates to the live ``skill_safety.check`` — the same classifiers
    ``skill_proposals.propose`` runs when a proposal is stored — so a
    patched skill that adds a deny-pattern command or a command that is
    not auto-allowable is refused here exactly as it is at storage time.
    Returns ``(True, [])`` for a clean text, ``(False, findings)``
    otherwise. A checker that cannot run is never a silent pass; it
    reports a finding (skill_safety's own invariant).
    """
    from .skill_safety import check as _skill_safety_check  # noqa: PLC0415

    findings = [str(f) for f in _skill_safety_check(new_text or "")]
    return (not findings), findings
