"""Fuzzy edit-block matching, inspired by Aider's editblock coder.

When an LLM proposes an `edit_file(old_string=..., new_string=...)` call,
exact-string match is the gold standard. But LLMs frequently emit
old_string with subtle whitespace drift:

    * leading indent off by 4 spaces (got the function body, forgot the
      extra indent inside an `if` block)
    * tabs converted to spaces (or vice versa) in the LLM's view of the file
    * trailing whitespace dropped on some lines
    * extra blank line between blocks

This module provides a single entry point ``fuzzy_replace(haystack, old,
new)`` that attempts progressively more lenient matching strategies and
returns the patched text plus a human-readable explanation of *what kind*
of fuzzy match succeeded — so the agent can correct its mental model on
the next call.

Strategies, tried in order:

1. **Exact** — included for completeness; same as ``str.replace``.
2. **Whitespace-normalized per-line** — strip trailing whitespace and
   collapse runs of `[ \\t]+` to a single space on each line, then look
   for a unique match window. If found, splice in ``new`` while
   re-indenting it to match the file's actual indent.
3. **Dedent + reindent** — strip the common leading whitespace from
   ``old`` and from each candidate window, compare the dedented forms.
   This catches the "indent drift" case where the LLM produced the right
   block but at the wrong indentation level.

The matcher REJECTS the edit when the fuzzy form matches more than one
window (ambiguous). It is up to the caller to either error out or ask the
user — silently picking one would be unsafe.

The caller is also responsible for honoring ``replace_all``; this module
deals with the *first* fuzzy match only. (For multi-occurrence cases
exact match is reliable enough; fuzzy + replace_all is asking for
trouble.)
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Optional


_WS_RUN = re.compile(r"[ \t]+")


def _leading_ws(s: str) -> str:
    return s[: len(s) - len(s.lstrip(" \t"))]


def _common_leading_ws(lines: list[str]) -> str:
    """Longest leading whitespace shared by all non-blank lines."""
    non_blank = [ln for ln in lines if ln.strip()]
    if not non_blank:
        return ""
    prefix = _leading_ws(non_blank[0])
    for ln in non_blank[1:]:
        lw = _leading_ws(ln)
        # Shrink prefix to the common head of (prefix, lw).
        i = 0
        while i < len(prefix) and i < len(lw) and prefix[i] == lw[i]:
            i += 1
        prefix = prefix[:i]
        if not prefix:
            break
    return prefix


def _normalize_line(ln: str) -> str:
    """Collapse internal whitespace runs to one space, strip trailing."""
    return _WS_RUN.sub(" ", ln).rstrip()


def _dedent_lines(lines: list[str], indent: str) -> list[str]:
    if not indent:
        return list(lines)
    out = []
    for ln in lines:
        if ln.startswith(indent):
            out.append(ln[len(indent):])
        elif not ln.strip():
            out.append(ln.lstrip(" \t"))  # blank lines: drop indent
        else:
            out.append(ln)
    return out


def _reindent(text: str, from_indent: str, to_indent: str) -> str:
    """Replace each line's `from_indent` prefix with `to_indent`.

    Lines that don't start with `from_indent` are still given the
    `to_indent` prefix if they have non-whitespace content; blank lines
    are left empty.
    """
    if from_indent == to_indent:
        return text
    keep_trailing_nl = text.endswith("\n")
    lines = text.split("\n")
    if keep_trailing_nl:
        # text.split("\n") with trailing "\n" gives empty last element;
        # drop it so we can rejoin cleanly later.
        lines = lines[:-1]
    out = []
    for ln in lines:
        if not ln.strip():
            out.append("")
            continue
        if from_indent and ln.startswith(from_indent):
            # Strip exactly the common-old-indent, prepend new common indent.
            out.append(to_indent + ln[len(from_indent):])
        elif not from_indent:
            # No old common indent to strip — line's existing leading
            # whitespace IS the relative indent within the block. Just
            # prepend the new common indent.
            out.append(to_indent + ln)
        else:
            # Line lacks the expected from_indent (rare; mismatched
            # indentation within `new`). Conservative: strip all leading
            # whitespace and prepend to_indent.
            out.append(to_indent + ln.lstrip(" \t"))
    rejoined = "\n".join(out)
    if keep_trailing_nl:
        rejoined += "\n"
    return rejoined


@dataclass
class FuzzyMatch:
    """Result of a successful fuzzy match.

    ``new_text`` is the full file with the splice applied. ``strategy``
    names the lenient strategy that produced it ("ws_normalized" or
    "dedent_reindent") so the tool result can mention it. ``indent_shift``
    describes the indent change applied to ``new_string`` (empty if
    none).
    """

    new_text: str
    strategy: str
    indent_shift: str
    matched_text: str  # the actual text in the file that matched


def fuzzy_replace(
    haystack: str, old: str, new: str
) -> Optional[FuzzyMatch]:
    """Try fuzzy strategies; return None if none yields a unique match.

    Caller should only invoke this AFTER confirming ``old`` does not
    appear exactly in ``haystack``.
    """
    if not old:
        return None

    old_lines = old.splitlines(keepends=False)
    if not old_lines:
        return None
    n = len(old_lines)

    haystack_lines_kept = haystack.splitlines(keepends=True)
    haystack_lines = haystack.splitlines(keepends=False)

    if len(haystack_lines) < n:
        return None

    # Strategy 2: whitespace-normalized per-line match.
    norm_old = [_normalize_line(ln) for ln in old_lines]

    ws_hits: list[int] = []
    for i in range(len(haystack_lines) - n + 1):
        window = haystack_lines[i : i + n]
        if [_normalize_line(ln) for ln in window] == norm_old:
            ws_hits.append(i)

    if len(ws_hits) == 1:
        i = ws_hits[0]
        # Compute char offsets in original haystack using the keepends list.
        start = sum(len(line) for line in haystack_lines_kept[:i])
        end = start + sum(len(line) for line in haystack_lines_kept[i : i + n])
        matched_text = haystack[start:end]

        # Detect indent shift: compare common indent of file window vs old.
        file_indent = _common_leading_ws(haystack_lines[i : i + n])
        old_indent = _common_leading_ws(old_lines)
        new_text_block = _reindent(new, old_indent, file_indent) if new else ""

        # Preserve a trailing newline if the matched block ended with one
        # but `new` does not (so we don't accidentally collapse two lines
        # into one).
        if matched_text.endswith("\n") and not new_text_block.endswith("\n"):
            new_text_block = new_text_block + "\n"

        patched = haystack[:start] + new_text_block + haystack[end:]
        return FuzzyMatch(
            new_text=patched,
            strategy="ws_normalized",
            indent_shift=(
                f"{old_indent!r} -> {file_indent!r}"
                if old_indent != file_indent else ""
            ),
            matched_text=matched_text,
        )

    if len(ws_hits) > 1:
        # Ambiguous — refuse rather than guess.
        return None

    # Strategy 3: dedent + reindent. Strip common indent from `old` and
    # compare against dedented windows.
    old_indent = _common_leading_ws(old_lines)
    old_dedented = _dedent_lines(old_lines, old_indent)
    norm_old_ded = [_normalize_line(ln) for ln in old_dedented]

    dr_hits: list[tuple[int, str]] = []
    for i in range(len(haystack_lines) - n + 1):
        window = haystack_lines[i : i + n]
        win_indent = _common_leading_ws(window)
        win_ded = _dedent_lines(window, win_indent)
        if [_normalize_line(ln) for ln in win_ded] == norm_old_ded:
            dr_hits.append((i, win_indent))

    if len(dr_hits) == 1:
        i, win_indent = dr_hits[0]
        start = sum(len(line) for line in haystack_lines_kept[:i])
        end = start + sum(len(line) for line in haystack_lines_kept[i : i + n])
        matched_text = haystack[start:end]
        new_block = _reindent(new, old_indent, win_indent) if new else ""
        if matched_text.endswith("\n") and not new_block.endswith("\n"):
            new_block = new_block + "\n"
        patched = haystack[:start] + new_block + haystack[end:]
        return FuzzyMatch(
            new_text=patched,
            strategy="dedent_reindent",
            indent_shift=(
                f"{old_indent!r} -> {win_indent!r}"
                if old_indent != win_indent else ""
            ),
            matched_text=matched_text,
        )

    return None
