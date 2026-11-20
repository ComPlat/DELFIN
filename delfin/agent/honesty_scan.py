"""The change records the honesty scan reads at turn end.

Input: the engine (its session key, its workspace root, and the mark
taken when the turn began) and the undo journal that the write path
fills for every journalled edit.

Output: a list of ``(path, new_start, new_end, source)`` records, the
shape :func:`delfin.agent.test_honesty.check` takes and
:func:`delfin.agent.verify_guard.scan_for_dishonest_tests` forwards.
``path`` is workspace-relative (the scan's ``test_zz`` rule is anchored
at ``^tests/``), ``source`` is the file's NEW text, and
``new_start..new_end`` are the lines this turn ADDED, so a pre-existing
``except:`` elsewhere in the same file stays quiet.

Mechanism: the journal keeps a byte-exact pre-image per change
(``change_journal``), so the change itself is ``pre-image -> the file on
disk now``; a ``difflib`` opcode pass gives the added line range.
:func:`mark` records the journal's last sequence number at turn start,
and only records above it are read, so a finding is never about another
turn's work. That is what makes the report arrive while the work is
still cheap to fix.

Why read the journal rather than reimplement it: the pre-image is
already stored and already diffed by ``edit_signals`` and
``pending_changes``; a second store of the same history would be a
second mechanism for one job.

Rejected: scanning the whole session's journal (it re-reports findings
for work the user has since accepted, and grows without bound); reading
a vcs diff (the workspace need not be a git repository, and work the
turn has not committed would be missing from it); reporting only the
paths with no line range (kind 1 would then fire on any old ``except:``
in a file the turn touched — the false alarm the diff-scoping in
``test_honesty`` exists to prevent).

Cost: one ``difflib`` pass per journalled text change of the turn.
Never raises: a broken journal degrades to no records, and therefore to
no report.
"""

from __future__ import annotations

import difflib
import os
from pathlib import Path
from typing import Any, Optional

#: Text above this is not diffed or parsed: a change larger than one turn
#: can be reviewed as a test change is not the case this scan is for.
MAX_TEXT = 400_000


def mark(engine: Any) -> int:
    """The journal's last sequence number now — take it at TURN START.

    Records whose ``seq`` is above the returned mark are this turn's.
    Returns 0 when the journal is unreadable, which makes the scan read
    the whole session rather than nothing; the caller never has to
    distinguish the two to stay safe.
    """
    try:
        from . import change_journal as _cj

        session = str(getattr(engine, "_edit_signal_key", lambda: "")() or "")
        if not session:
            return 0
        records = _cj.list_changes(session) or []
        return int(records[-1].get("seq", 0) or 0) if records else 0
    except Exception:
        return 0


def _added_range(pre: Optional[str], post: str, created: bool) -> tuple:
    """The ``(first, last)`` line numbers this change added to *post*.

    A created file is new from its first line to its last. A change whose
    pre-image is unusable (truncated, binary, dropped at the cap) is
    reported as the whole file, the reading that cannot miss a new
    violation; a change that only removed lines adds nothing and is
    reported as line 1..1, so nothing outside it is judged.
    """
    total = post.count("\n") + (0 if post.endswith("\n") or not post else 1)
    if created or pre is None:
        return (1, total)
    a = pre[:MAX_TEXT].splitlines()
    b = post[:MAX_TEXT].splitlines()
    first = last = 0
    for tag, _i1, _i2, j1, j2 in difflib.SequenceMatcher(
            None, a, b, autojunk=False).get_opcodes():
        if tag in ("replace", "insert") and j2 > j1:
            first = j1 + 1 if not first else first
            last = max(last, j2)
    return (first, last) if first else (1, 1)


def _rel(path: str, root: str) -> str:
    """*path* relative to *root* when it is inside it, else unchanged."""
    try:
        return os.path.relpath(path, root) if root else path
    except ValueError:
        return path


def _pre_text(sdir: Path, rec: dict) -> Optional[str]:
    """A journalled record's pre-image as text, or None when unusable."""
    try:
        if rec.get("truncated") or rec.get("binary") or rec.get("dropped"):
            return None
        data = (sdir / str(rec.get("pre_file", ""))).read_bytes()
    except OSError:
        return None
    try:
        return data.decode("utf-8")
    except UnicodeDecodeError:
        return None


def turn_records(engine: Any, since: int = 0) -> list:
    """Records for the changes journalled after *since*; ``[]`` on any
    failure. See the module docstring for the shape and the reason."""
    try:
        from . import change_journal as _cj

        session = str(getattr(engine, "_edit_signal_key", lambda: "")() or "")
        if not session:
            return []
        root = str(getattr(engine, "repo_dir", "") or "")
        sdir = _cj._session_dir(session)
        out: list = []
        for rec in _cj.list_changes(session) or []:
            if int(rec.get("seq", 0) or 0) <= int(since or 0):
                continue
            if rec.get("undone") or rec.get("deleted"):
                continue
            path = str(rec.get("path", "") or "")
            if not path:
                continue
            try:
                post = Path(path).read_text(encoding="utf-8")
            except (OSError, UnicodeDecodeError):
                continue          # gone or unreadable now: nothing to judge
            if len(post) > MAX_TEXT:
                continue
            created = bool(rec.get("created"))
            pre = None if created else _pre_text(sdir, rec)
            first, last = _added_range(pre, post, created)
            out.append((_rel(path, root), first, last, post))
        return out
    except Exception:
        return []


__all__ = ["turn_records", "mark", "MAX_TEXT"]
