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

The correction heuristic is ``skill_learning._looks_like_correction`` —
the canonical marker list the codebase already uses for
"this user message reads as a course correction". Reused here verbatim so
the two recognition surfaces stay consistent; the only thing added is the
scope (a skill was used), the user's own words, and the source gate.
"""

from __future__ import annotations

from typing import Optional

from .skill_learning import _looks_like_correction

__all__ = ["feedback", "patched_skill_safe"]

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
        _looks_like_correction(text)
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
