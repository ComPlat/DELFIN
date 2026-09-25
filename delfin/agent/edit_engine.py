"""The decision core of an edit: text in, result out. No file I/O.

Everything ``edit_file`` / ``multi_edit`` must decide before any
permission, path or write is involved lives here, pure:

* ``apply_edit`` - exact ``old_string`` -> ``new_string``, once or
  ``replace_all``. An ambiguous match without ``replace_all`` is an
  error that names the line numbers of ALL matches (today's message
  names only the count - the s7 report).
* near-miss diagnosis - when ``old_string`` does not match exactly, the
  engine says WHERE it almost stands (indentation shift, trailing
  whitespace, tabs vs spaces) with the line number and the actual text
  there. It NEVER replaces approximately: the drifted text survives
  untouched. This replaces the silent whitespace-tolerant fallback in
  api_client.py (``editblock.fuzzy_replace`` at the ``edit_file`` call
  site), which applied drifted matches without asking.
* ``apply_multi_edit`` - all edits applied against the intermediate
  state; if any fails, NONE is applied and the result names the failing
  index and why. (api_client's in-memory loop already behaves this way;
  this module keeps the guarantee for the future.)
* python syntax regression - if the text parsed with ``ast.parse``
  before the edit and does not after, the syntax error's line is
  reported in ``syntax_regression``. Reported only: whether that blocks
  the edit is the caller's decision (see the docstring of EditResult).
* shape preservation - nothing outside a matched occurrence changes;
  CRLF pairs and the final-newline state survive as they are.

The engine operates on the NORMALISED text that
``text_files.read_text_file`` hands out (LF endings, no BOM). Restoring
the file's own convention is the writer's job, not ours.
"""

from __future__ import annotations

import ast
import re
from dataclasses import dataclass, field
from typing import Optional, Sequence


@dataclass(frozen=True)
class NearMiss:
    """Where an old_string almost matched, and how it differs.

    ``line`` is 1-based, the first line of the near-matching region.
    ``actual_text`` is the verbatim text at that spot (including its
    line ending when the region spans lines) - the agent reads it,
    copies it verbatim, and its next edit matches.
    """

    line: int
    reason: str
    actual_text: str


@dataclass(frozen=True)
class SyntaxRegression:
    """A file that parsed before the edit and does not after.

    ``line`` is 1-based, straight from the SyntaxError. The engine
    REPORTS this and still returns the edited text; blocking the edit
    is a policy decision that belongs to the caller.
    """

    line: int
    message: str


@dataclass(frozen=True)
class EditResult:
    """One edit's outcome. Either the edit landed, or the fields say
    exactly why not.

    On success: ``applied=True``, ``new_text`` the whole edited text,
    ``error=None``, ``near_misses`` empty, ``match_lines`` the lines
    that were replaced.

    On failure: ``applied=False``, ``new_text=None`` (NEVER a partial or
    approximate result), ``error`` the human-readable reason,
    ``near_misses`` what the diagnosis found (may be empty when nothing
    is close), ``match_lines`` the lines of every ambiguous match.

    ``syntax_regression`` is independent of applied: an edit that lands
    can still leave Python unparseable, and the caller decides what to
    do about it.
    """

    applied: bool
    new_text: Optional[str]
    error: Optional[str] = None
    near_misses: Sequence[NearMiss] = ()
    match_lines: Sequence[int] = ()
    failed_index: Optional[int] = None
    syntax_regression: Optional[SyntaxRegression] = None


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _match_lines(text: str, old: str) -> list[int]:
    """1-based line numbers of every start of an exact ``old`` match."""
    lines: list[int] = []
    # Line starts let us translate a char offset to a line number
    # without scanning for every match.
    starts = [0]
    for i, ch in enumerate(text):
        if ch == "\n":
            starts.append(i + 1)
    pos = 0
    while True:
        idx = text.find(old, pos)
        if idx < 0:
            break
        # largest start <= idx
        lo, hi = 0, len(starts) - 1
        while lo < hi:
            mid = (lo + hi + 1) // 2
            if starts[mid] <= idx:
                lo = mid
            else:
                hi = mid - 1
        lines.append(lo + 1)
        pos = idx + 1
    return lines


def _line_start_offsets(text: str) -> list[int]:
    starts = [0]
    for i, ch in enumerate(text):
        if ch == "\n":
            starts.append(i + 1)
    return starts


def _offset_to_line(starts: list[int], offset: int) -> int:
    lo, hi = 0, len(starts) - 1
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if starts[mid] <= offset:
            lo = mid
        else:
            hi = mid - 1
    return lo + 1


_WS_RUN = re.compile(r"[ \t]+")


def _norm(line: str) -> str:
    """Collapse internal whitespace runs to one space, strip the ends."""
    return _WS_RUN.sub(" ", line).strip()


def _leading_ws(line: str) -> str:
    return line[: len(line) - len(line.lstrip(" \t"))]


def _diagnose_near_misses(text: str, old: str) -> list[NearMiss]:
    """Find regions that match ``old`` modulo whitespace drift.

    Reuses the two lenient comparison forms from editblock
    (whitespace-normalised per line; dedented) but ONLY to diagnose:
    every hit becomes a NearMiss naming line, reason and the verbatim
    text, and nothing is ever replaced. Capped at a few hits so a
    pathological file cannot flood the result.
    """
    if not old:
        return []
    old_lines = old.splitlines()
    if not old_lines or not any(ln.strip() for ln in old_lines):
        return []
    n = len(old_lines)
    text_lines = text.splitlines()
    if len(text_lines) < n:
        return []

    norm_old = [_norm(ln) for ln in old_lines]
    starts = _line_start_offsets(text)

    misses: list[NearMiss] = []
    seen_spans: set[tuple[int, int]] = set()

    def _record(i: int, line: int, reason: str) -> None:
        # span of lines [i, i+n) as char offsets in text
        start_off = starts[i]
        end_off = (
            starts[i + n] if i + n < len(starts) else len(text)
        )
        span = (start_off, end_off)
        if span in seen_spans:
            return
        seen_spans.add(span)
        misses.append(NearMiss(
            line=line,
            reason=reason,
            actual_text=text[start_off:end_off],
        ))

    def _first_diff(window: list[str]) -> int:
        for j, (a, b) in enumerate(zip(window, old_lines)):
            if a != b:
                return j
        return 0

    for i in range(len(text_lines) - n + 1):
        window = text_lines[i:i + n]
        if [_norm(ln) for ln in window] != norm_old:
            continue
        # Same modulo whitespace: name the differences.
        reasons = []
        old_indent = _leading_ws(old_lines[0])
        win_indent = _leading_ws(window[0])
        if old_indent != win_indent:
            if "\t" in old_indent + win_indent and (
                    old_indent.strip("\t") == ""
                    or win_indent.strip("\t") == ""):
                reasons.append(
                    f"tabs vs spaces: file has {win_indent!r}, "
                    f"old_string has {old_indent!r}")
            else:
                reasons.append(
                    f"indentation: file has {win_indent!r}, old_string "
                    f"has {old_indent!r}")
        for a, b in zip(window, old_lines):
            if a != b and _norm(a) == _norm(b):
                la, lb = _leading_ws(a), _leading_ws(b)
                if la != lb:
                    if "\t" in la + lb:
                        reasons.append(
                            f"tabs vs spaces: file has {la!r}, "
                            f"old_string has {lb!r}")
                    else:
                        reasons.append(
                            f"indentation: file has {la!r}, "
                            f"old_string has {lb!r}")
                if a.rstrip() == b.rstrip():
                    reasons.append("trailing whitespace differs")
                break
        if not reasons:
            reasons.append("whitespace differs (normalized form matches)")
        _record(i, i + 1 + _first_diff(window),
                "; ".join(dict.fromkeys(reasons)))

    if misses:
        return misses[:4]

    # Dedented form: old at a different indent level entirely.
    def _dedent(lines: list[str]) -> list[str]:
        non_blank = [ln for ln in lines if ln.strip()]
        if not non_blank:
            return list(lines)
        prefix = _leading_ws(non_blank[0])
        for ln in non_blank[1:]:
            lw = _leading_ws(ln)
            j = 0
            while j < len(prefix) and j < len(lw) and prefix[j] == lw[j]:
                j += 1
            prefix = prefix[:j]
            if not prefix:
                break
        out = []
        for ln in lines:
            if ln.startswith(prefix):
                out.append(ln[len(prefix):])
            elif not ln.strip():
                out.append("")
            else:
                out.append(ln)
        return out

    norm_old_ded = [_norm(ln) for ln in _dedent(old_lines)]
    for i in range(len(text_lines) - n + 1):
        window = text_lines[i:i + n]
        win_ded = _dedent(window)
        if [_norm(ln) for ln in win_ded] != norm_old_ded:
            continue
        wi = _leading_ws(window[0]) if window[0].strip() else ""
        oi = _leading_ws(old_lines[0]) if old_lines[0].strip() else ""
        _record(i, i + 1 + _first_diff(window),
                f"indentation: file has {wi!r}, old_string has {oi!r}")
        if len(misses) >= 4:
            break
    return misses


def _syntax_regression(
    before: str, after: str, is_python: bool,
) -> Optional[SyntaxRegression]:
    """ast.parse(before) ok and ast.parse(after) not -> the error, else None.

    A file that did not parse before cannot "regress" - pre-existing
    damage is not this edit's doing, and reporting it as such would send
    the agent hunting a bug it did not introduce.
    """
    if not is_python:
        return None
    try:
        ast.parse(before)
    except SyntaxError:
        return None
    try:
        ast.parse(after)
    except SyntaxError as exc:
        return SyntaxRegression(
            line=(exc.lineno or 0) or 1, message=exc.msg or str(exc))
    return None


# ---------------------------------------------------------------------------
# public API
# ---------------------------------------------------------------------------


def apply_edit(
    text: str,
    old_string: str,
    new_string: str,
    *,
    replace_all: bool = False,
    is_python: bool = False,
) -> EditResult:
    """Apply one exact replacement to ``text``. See the module docstring
    for the full contract; on failure ``new_text`` is None and the
    fields say why."""
    if not old_string:
        return EditResult(
            applied=False, new_text=None,
            error="old_string is required")
    if old_string == new_string:
        return EditResult(
            applied=False, new_text=None,
            error="new_string must differ from old_string")

    count = text.count(old_string)
    if count == 0:
        snippet = old_string if len(old_string) <= 60 else (
            old_string[:57] + "...")
        misses = tuple(_diagnose_near_misses(text, old_string))
        parts = [f"old_string {snippet!r} not found (no exact match)."]
        if misses:
            parts.append(" It nearly stands at:")
            for m in misses:
                parts.append(
                    f"\n  line {m.line}: {m.reason}; actual text there:"
                    f"\n{m.actual_text}")
        parts.append(
            " Re-read the file and copy the target block verbatim."
            " Nothing was replaced.")
        return EditResult(
            applied=False, new_text=None, error="".join(parts),
            near_misses=misses)
    if count > 1 and not replace_all:
        lines = _match_lines(text, old_string)
        return EditResult(
            applied=False, new_text=None,
            error=(
                f"old_string matches {count} times, on lines "
                f"{', '.join(str(n) for n in lines)}. Provide more "
                "surrounding context to make it unique, or set "
                "replace_all=True."),
            match_lines=lines)

    new_text = (
        text.replace(old_string, new_string)
        if replace_all else text.replace(old_string, new_string, 1)
    )
    return EditResult(
        applied=True, new_text=new_text,
        match_lines=_match_lines(text, old_string),
        syntax_regression=_syntax_regression(
            text, new_text, is_python),
    )


def apply_multi_edit(
    text: str,
    edits: Sequence[dict],
    *,
    is_python: bool = False,
) -> EditResult:
    """Apply a list of edits atomically against the intermediate state.

    Every edit is validated and applied to the in-memory text in order;
    the FIRST failure aborts the batch with applied=False, new_text=None
    and failed_index set (0-based). Only when all edits succeed is the
    final text returned - the caller writes once, or not at all.
    """
    if not isinstance(edits, Sequence) or not edits:
        return EditResult(
            applied=False, new_text=None,
            error="edits must be a non-empty list")
    current = text
    for i, ed in enumerate(edits):
        if not isinstance(ed, dict):
            return EditResult(
                applied=False, new_text=None,
                error=f"edit #{i + 1} is not an object",
                failed_index=i)
        r = apply_edit(
            current,
            ed.get("old_string", "") or "",
            ed.get("new_string", "") or "",
            replace_all=bool(ed.get("replace_all", False)),
            is_python=False,  # checked once below, on the final text
        )
        if not r.applied:
            return EditResult(
                applied=False, new_text=None,
                error=f"edit #{i + 1}: {r.error}",
                near_misses=r.near_misses,
                match_lines=r.match_lines,
                failed_index=i)
        current = r.new_text
    return EditResult(
        applied=True, new_text=current,
        syntax_regression=_syntax_regression(text, current, is_python),
    )
