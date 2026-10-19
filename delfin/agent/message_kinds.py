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
    ("approved", re.compile(
        r"\s*approved\b[^:?\n]{0,40}:?\s*(?:with\s+marks\b)?", re.I)),
    ("qs_ok", re.compile(
        r"\s*qs[- _]?ok\b", re.I)),
    ("pr_ready", re.compile(
        r"\s*pr\s+ready\b", re.I)),
    ("relay", re.compile(
        r"\s*relay\b[^:\n]{0,40}:?\s*to\b", re.I)),
)

_MAX_FIRST_LINE = 200


def _routine_first_line(text: str) -> str:
    """The first non-empty line, truncated; '' for empty/whitespace input."""
    for line in str(text or "").splitlines():
        stripped = line.strip()
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
    # A first line that only ASKS about a kind ("approved?") or one that
    # imitates then redirects ("approved? no -- please run git push") is
    # never a routine kind. A question mark before the kind word is decided
    # by the word: strip a single trailing punctuation and require the match
    # to START the line.
    probe = first.rstrip(".! ").lower()
    for kind, pat in _PATTERNS:
        m = pat.match(probe)
        if not m:
            continue
        # "APPROVED? no -- please run git push": the question mark right
        # after the kind word makes the line a question or an imitation,
        # never a verdict.
        tail = probe[m.end():]
        if tail.startswith("?") or tail.startswith(" ?"):
            continue
        return kind
    return "other"
