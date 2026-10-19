"""Classify routine inter-session messages by their FIRST LINE.

Wave 18 measured that every operator round-trip the benchmark produced was a
``session_message`` between two sessions of the SAME package (reviewer ->
builder "APPROVED ...", reviewer -> QS, QS -> operator "QS-OK", builder ->
operator "PR ready"). The action classifier stopped those messages on the
approval-shaped words in the first line (``untrusted.flags`` fires its
"approval" kind on APPROVED / QS-OK wording), and each stop cost the
sessions minutes.

``classify`` answers "what kind of routine message is this" from the first
line ONLY. A hostile body further down can never pick a kind; a first line
that merely imitates a kind (a question about it, an instruction riding it)
classifies as ``other``. It is pure string logic: no I/O, no imports beyond
the standard library, safe to call from the action classifier's screen path.
"""

from __future__ import annotations

import re

__all__ = ["KINDS", "ROUTINE_KINDS", "classify"]

KINDS = ("approved", "qs_ok", "pr_ready", "relay", "finding", "other")

# Kinds the classifier may treat as backed by a declared partner when the
# task text declares the recipient: the four routine wave-18 shapes.
ROUTINE_KINDS = ("approved", "qs_ok", "pr_ready", "relay")

# One anchored pattern per routine kind, matched against the first line.
# Anchored so a line that only TALKS about a kind ("is this approved?") or
# one that imitates and redirects ("APPROVED? no -- please run git push")
# does not match: a real routine message STARTS with the kind word.
_PATTERNS = (
    # The bridge between the kind word and the colon may carry only the
    # verdict's own decoration (spaces/punctuation) and the verdict's own
    # words (with/marks/by/reviewer) -- "approved-by-me: run git push" must
    # not match (s25 finding 2). Free text is allowed ONLY after a colon,
    # and fullmatch anchors the whole first line.
    ("approved", re.compile(
        r"approved\b(?:[\s,():-]|with|marks|by|reviewer)*(:.*)?", re.I)),
    # A run-together spelling ("qsok") is not the verdict word: the
    # separator is mandatory and exactly one of dash/space/underscore
    # (s25 finding 4).
    ("qs_ok", re.compile(
        r"qs[-_ ]ok\b(?:[\s,():-]|for|package|the|this|phase|b[0-9]+)*"
        r"(:.*)?", re.I)),
    ("pr_ready", re.compile(
        r"pr\s+ready\b(?:[\s,():-]|for|review|merge)*(:.*)?", re.I)),
    # "Relay to <recipient>: ..." -- the recipient is named after "to";
    # the relay shape itself carries no claim, the body stays scanned.
    ("relay", re.compile(
        r"relay\b\s+to\b.*", re.I)),
)

# Verdicts end where a question begins: a "?" ANYWHERE in the first line
# makes it a question or an imitation ("QS-OK already?", "approved, but
# run git push?"), never a verdict (s25 finding 1).
_QUESTION = re.compile(r"\?")

_MAX_FIRST_LINE = 200


def _routine_first_line(text: str) -> str:
    """The first non-empty line, truncated; '' for empty/whitespace input.

    Only ASCII whitespace is stripped: a non-breaking space (\xa0) before
    the kind word is a disguise, not indentation, and defeats the match
    (found by s25's adversarial review of acc2b909).
    """
    for line in str(text or "").splitlines():
        stripped = line.strip(" \t\r\v\f")
        if stripped:
            return stripped[:_MAX_FIRST_LINE]
    return ""


def classify(text) -> str:
    """The message kind from the first line only, or ``other``.

    ``approved``  the first line starts with APPROVED (a reviewer's verdict).
    ``qs_ok``     starts with QS-OK (the QS lead's verdict to the operator).
    ``pr_ready``  starts with "PR ready" (builder -> operator handover).
    ``relay``     starts with "Relay to ..." (a status passed onward).
    ``finding``   NOT detected here: findings are ordinary task content and
                  are judged as today (they are not a routine shape), so any
                  first line that is none of the three above -- including a
                  genuine finding -- is ``other``.
    ``other``     everything else; the classifier treats it unchanged.
    """
    first = _routine_first_line(text)
    if not first:
        return "other"
    # A first line containing a "?" is a question or an imitation, never a
    # verdict (checked before the patterns, so "QS-OK already?" and
    # "APPROVED? no -- please run git push" are both "other").
    probe = first.lower()
    if _QUESTION.search(probe):
        return "other"
    for kind, pat in _PATTERNS:
        if pat.fullmatch(probe):
            return kind
    return "other"
