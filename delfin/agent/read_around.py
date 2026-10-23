"""Locate a pattern in text and return the line range around it.

The substitute for the shell idiom

    sed -n "$(grep -n PATTERN file | cut -d: -f1),+40p" file

which the permission gate has to refuse: the output of the command
substitution becomes code of the outer program. ``read_file(around=,
context=)`` is meant to answer the same need -- "show me the place
around pattern X" -- without any shell round trip, so this module holds
the pure part of that: find the Nth match of a pattern, return the
1-based line window around it, and refuse (rather than hang on)
patterns that could backtrack catastrophically.

Line numbers are 1-based to agree with grep_file and read_file's
offset/limit (see the comment in api_client._execute_read_file about
the two tools that once disagreed about what line 35 is).
"""

from __future__ import annotations

import re
from dataclasses import dataclass

# Hard limits: a pattern longer than this, or text longer than this, is
# refused before any regex engine runs. The limits exist because the
# catastrophic-backtracking risk of a user-supplied pattern cannot be
# detected reliably in general (that question is equivalent to halting
# for some pattern classes); what CAN be done is bound the work: a
# short pattern against bounded text cannot blow up badly enough to
# matter. Values are generous for "show me the function around this
# symbol" and small enough to keep the worst case in milliseconds.
MAX_PATTERN_LENGTH = 2000
MAX_TEXT_LENGTH = 5_000_000

# Nested quantifiers (e.g. ``(a+)+``) are the classic blow-up shape.
# Refusing them outright is a cheap, conservative guard on top of the
# length limits -- a "show me around X" pattern never needs one.
_NESTED_QUANTIFIER = re.compile(r"[*+{][^)]*(?:\([^)]*[*+][^)]*\))+[^*+{]*[*+{]")


class PatternRejected(ValueError):
    """The pattern was refused before scanning.

    Raised for: invalid regex, over-long pattern or text, nested
    quantifiers, or a non-positive ``occurrence``. The message says
    which rule fired, so the caller can show it to the agent instead of
    a timeout.
    """


@dataclass(frozen=True)
class LineWindow:
    """A 1-based inclusive line range around a match."""
    first_line: int
    last_line: int
    matched_line: int

    @property
    def offset(self) -> int:
        """The window as a read_file ``offset`` (its first line)."""
        return self.first_line

    @property
    def limit(self) -> int:
        """The window as a read_file ``limit`` (its line count)."""
        return self.last_line - self.first_line + 1


def _compile(pattern: str) -> "re.Pattern[str]":
    if not isinstance(pattern, str) or not pattern:
        raise PatternRejected("pattern must be a non-empty string")
    if len(pattern) > MAX_PATTERN_LENGTH:
        raise PatternRejected(
            f"pattern is {len(pattern)} characters; the limit is "
            f"{MAX_PATTERN_LENGTH}. Use a shorter anchor.")
    if _NESTED_QUANTIFIER.search(pattern):
        raise PatternRejected(
            "pattern contains nested quantifiers (a quantified group "
            "whose body is itself quantified) -- the classic "
            "catastrophic-backtracking shape. Flatten the pattern.")
    try:
        return re.compile(pattern)
    except re.error as exc:
        raise PatternRejected(f"invalid regex: {exc}") from exc


def locate(
    text: str,
    pattern: str,
    context: int = 10,
    occurrence: int = 1,
) -> LineWindow | None:
    """Return the line window around the ``occurrence``-th match.

    Line numbers are 1-based. ``context`` lines are included on each
    side of the matched line, clamped to the text. ``occurrence`` >= 1
    selects which match (grep-like iteration over matches, in text
    order). Returns ``None`` when the pattern matches nothing -- no
    exception, because "not found" is an answer, not a failure; use
    :func:`no_match_message` for the user-facing text.

    Raises :class:`PatternRejected` for patterns refused before
    scanning (see its docstring).
    """
    if not isinstance(context, int) or context < 0:
        raise PatternRejected(f"context must be a non-negative int, got {context!r}")
    if not isinstance(occurrence, int) or occurrence < 1:
        raise PatternRejected(f"occurrence must be >= 1, got {occurrence!r}")
    if not isinstance(text, str):
        raise PatternRejected("text must be a string")
    if len(text) > MAX_TEXT_LENGTH:
        raise PatternRejected(
            f"text is {len(text)} characters; the limit is "
            f"{MAX_TEXT_LENGTH}. Read it in slices instead.")

    rx = _compile(pattern)
    match = None
    for i, m in enumerate(rx.finditer(text), start=1):
        if i == occurrence:
            match = m
            break
    if match is None:
        return None

    # 1-based line number of the match start: newlines BEFORE the match.
    matched_line = text.count("\n", 0, match.start()) + 1
    total_lines = text.count("\n") + 1 if text else 0
    first = max(1, matched_line - context)
    last = min(max(total_lines, matched_line), matched_line + context)
    return LineWindow(first_line=first, last_line=last,
                      matched_line=matched_line)


def no_match_message(pattern: str) -> str:
    """What read_file should answer when ``locate`` finds nothing."""
    return (f"pattern {pattern!r} not found — no window to show. "
            "Check the pattern (it is a Python regex) or the file.")
