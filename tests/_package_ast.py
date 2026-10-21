"""One parse per package file and test session, shared by the whole-tree guards.

Input: a path to a ``.py`` file and the ``errors`` mode the caller would
have passed to ``read_text(encoding="utf-8", errors=...)``. Output: the
``ast.Module`` that ``ast.parse(path.read_text(encoding="utf-8",
errors=errors))`` returns, with the same exceptions (``UnicodeDecodeError``
under ``"strict"``, ``SyntaxError``); :func:`source` returns the decoded
text. An entry is reused only while the file's ``st_mtime_ns`` and
``st_size`` are unchanged, so a file written during the session is parsed
again. Failures are not cached: the next call raises afresh.

Why: eight guard tests in six files (solo safety invariants, argv/shell,
real-home constants, scheduler calls, containment labels, syntax
integrity) each parsed all ~805 files of ``delfin/``; one pass is 2.1 s
on an idle machine and the gate made about nine of them. The guards stay
in the fast run and keep their own predicates; only the parse is shared.

The first miss for a file under ``delfin/`` parses the whole package with
the cyclic garbage collector paused. Parsed one file at a time into a
growing cache, the 2.6 M long-lived nodes set off 14 full collections on
the way (7.8 s instead of 2.8 s for the same parse, measured with a
1.5 M-object heap standing in for the suite's).

Cost: holding every tree of the package is ~740 MB RSS (2.6 M nodes,
measured 2026-10-15). ``tests/conftest.py`` calls :func:`release` after
the last collected test whose module uses this one, so the trees live
from the first consumer to the last, not to the end of the session.

Callers must not modify a returned tree; it is shared.
"""

from __future__ import annotations

import ast
import gc
from pathlib import Path

_PACKAGE = Path(__file__).resolve().parents[1] / "delfin"

#: str(path) -> (mtime_ns, size, text, tree); only files that decode as
#: strict UTF-8 are cached, so every errors mode reads the same text.
_CACHE: dict[str, tuple[int, int, str, ast.Module]] = {}


_PRIMED: list[bool] = []


def _prime() -> None:
    """Parse every file of the package once, the collector paused."""
    _PRIMED.append(True)
    was_enabled = gc.isenabled()
    gc.disable()
    try:
        for one in sorted(_PACKAGE.rglob("*.py")):
            try:
                _entry(one, "strict")
            except (OSError, SyntaxError, UnicodeDecodeError, ValueError):
                continue
    finally:
        if was_enabled:
            gc.enable()


def _entry(path, errors: str):
    p = Path(path)
    st = p.stat()
    resolved = p.resolve()
    key = str(resolved)
    hit = _CACHE.get(key)
    if hit is not None and hit[0] == st.st_mtime_ns and hit[1] == st.st_size:
        return hit
    if not _PRIMED and resolved.is_relative_to(_PACKAGE):
        _prime()
        hit = _CACHE.get(key)
        if hit is not None and hit[0] == st.st_mtime_ns and hit[1] == st.st_size:
            return hit
    raw = p.read_bytes()
    try:
        text, cacheable = raw.decode("utf-8"), True
    except UnicodeDecodeError:
        if errors == "strict":
            raise
        text, cacheable = raw.decode("utf-8", errors=errors), False
    # read_text() decodes in text mode: universal newlines.
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    entry = (st.st_mtime_ns, st.st_size, text, ast.parse(text))
    if cacheable:
        _CACHE[key] = entry
    return entry


def parse(path, *, errors: str = "strict") -> ast.Module:
    """The tree of *path*, parsed at most once while the file is unchanged."""
    return _entry(path, errors)[3]


def source(path, *, errors: str = "strict") -> str:
    """The text the tree of :func:`parse` was made from."""
    return _entry(path, errors)[2]


def release() -> None:
    """Drop every cached tree."""
    _CACHE.clear()
    _PRIMED.clear()
