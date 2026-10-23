"""scope_route -- the one correct route for a write outside the session's
write scope.

When a write is refused because the path is not covered by
``write_allow_globs``, the sanctioned route is the operator-built
protected-file mechanism: two copies of the file under ``.gate/px/`` in the
workspace (``orig_`` = the file as it is, ``new_`` = the file with only the
wanted change) plus a one-line message to the operator, who applies the
change through the write gate. Every other route -- a script that writes,
``git checkout <branch> -- <file>`` -- reaches the same file around the
gate and is refused or a violation.

:func:`plan` returns ``None`` for a path inside the scope (nothing to
route) and a dict with ``orig_copy``, ``new_copy`` and ``message`` for a
path outside it.

Safety: the message is built from STATIC strings. The only dynamic part is
the file name, sanitised to a plain basename in ``[A-Za-z0-9_.+-]`` (dot
runs collapsed, no ``..``), so hostile text in a path -- newline injection,
``$()``, quotes, a traversal -- can never put a line, a command or a
directive into what the agent would send. An unsanitisable or
out-of-workspace path gets a generic message that names no path at all.
"""

from __future__ import annotations

import posixpath
import re

_SAFE_NAME = re.compile(r"[A-Za-z0-9_.+-]{1,64}")

#: A copy path inside ``.gate/px/``: relative to the workspace, plain.
_COPY_PREFIX = ".gate/px/"


def _glob_regex(pattern: str) -> re.Pattern:
    """A scope glob as a regex where ``**`` spans path separators and a
    single ``*`` and ``?`` do not (fnmatch's ``*`` crosses ``/``, which
    would let ``delfin/*`` swallow ``delfin/agent/x.py``)."""
    out = []
    i = 0
    while i < len(pattern):
        if pattern.startswith("**", i):
            out.append(".*")
            i += 2
        elif pattern[i] == "*":
            out.append("[^/]*")
            i += 1
        elif pattern[i] == "?":
            out.append("[^/]")
            i += 1
        else:
            out.append(re.escape(pattern[i]))
            i += 1
    return re.compile("".join(out) + r"\Z")


def _in_scope(rel_path: str, scope: tuple[str, ...]) -> bool:
    """True when ``rel_path`` (workspace-relative, '/'-separated) matches
    one scope glob. ``**`` spans any depth, ``*`` one level, like the
    write gate reads them. An EMPTY scope matches everything: that is the
    gate's own reading (_write_in_scope, api_client.py:13883-13884 -- no
    write scope set means no restriction), mirrored here so the two never
    disagree on which writes are allowed."""
    norm = posixpath.normpath(rel_path.replace("\\", "/"))
    if not (scope or ()):  # empty scope == no restriction (gate: :13883)
        return True
    for pattern in scope or ():
        pattern = (pattern or "").strip()
        if not pattern:
            continue
        if _glob_regex(pattern).fullmatch(norm):
            return True
    return False


def _safe_basename(path: str) -> str:
    """The sanitised basename of ``path``, or "" when there is none.

    Takes the LAST path component only (a directory name never reaches the
    message), collapses repeated dots, and keeps the match bounded.
    """
    name = posixpath.basename(path.replace("\\", "/").rstrip("/"))
    match = _SAFE_NAME.search(name)
    if not match:
        return ""
    cleaned = match.group(0)
    # "a..b" survives the charset; collapse dot runs so no ".." remains.
    while ".." in cleaned:
        cleaned = cleaned.replace("..", ".")
    return cleaned.lstrip(".").lstrip("-") or ""


def _copies_for(rel_path: str, workspace: str) -> tuple[str, str]:
    """The orig_/new_ copy paths for a workspace-relative path, kept
    distinct between files with the same basename by including the parent
    directory (sanitised the same way)."""
    parent = posixpath.dirname(posixpath.normpath(rel_path))
    stem = _safe_basename(rel_path)
    if not stem:
        stem = "unnamed"
    dirs = [d for d in parent.split("/") if d]
    tag = "_".join(_safe_basename(d) or "" for d in dirs[-2:])
    name = f"{tag}_{stem}" if tag else stem
    base = posixpath.join(workspace, _COPY_PREFIX)
    return (posixpath.join(base, f"orig_{name}"),
            posixpath.join(base, f"new_{name}"))


def plan(path: str, scope: tuple[str, ...], workspace: str) -> dict | None:
    """The route for writing ``path``, or ``None`` when it is in scope.

    ``path`` may be workspace-relative or absolute; ``scope`` is the
    ``write_allow_globs`` tuple; ``workspace`` is the session's workspace
    root (absolute or workspace-relative as the caller knows it).

    Returns a dict ``{"path", "orig_copy", "new_copy", "message"}``.
    ``orig_copy``/``new_copy`` are absolute paths under the workspace's
    ``.gate/px/`` when the path is a file inside the workspace, else
    ``None`` (the message then asks the operator directly, without
    pretending a host file is project source).
    """
    path = (path or "").strip()
    workspace = (workspace or "").strip()
    if not path:
        return None

    is_absolute = path.startswith("/") or path.startswith("~")
    if is_absolute:
        rel = posixpath.relpath(path, workspace) \
            if workspace and path.startswith(workspace + "/") else None
    else:
        rel = posixpath.normpath(path.replace("\\", "/"))

    if rel is not None and not rel.startswith(".."):
        if rel in (".", ""):
            # The path IS the workspace itself (e.g. "delfin/..") -- a
            # directory, not a file; the gate resolves it to the workspace
            # root. Nothing to plan.
            return None
        if _in_scope(rel, scope):
            return None
        copies = _copies_for(rel, workspace)
        # The message names the route, never the raw path text: every
        # fragment it echoes comes from _safe_basename, so hostile text
        # (newlines, $(), quotes) cannot land in what the agent sends.
        stem = _safe_basename(rel) or "unnamed"
        tag = posixpath.basename(copies[0]).removeprefix("orig_")
        return {
            "path": rel,
            "orig_copy": copies[0],
            "new_copy": copies[1],
            "message": (f"commit please: file {stem} -> two copies "
                        f"{tag} and new_{tag} under .gate/px/"),
        }

    # Absolute outside the workspace, or an unusable relative path: no
    # copies, a one-line generic message, no path text echoed.
    return {
        "path": path,
        "orig_copy": None,
        "new_copy": None,
        "message": ("commit please: a protected change outside this "
                    "workspace's write scope -- the operator applies it "
                    "from .gate/px/ copies; the path was not echoed for "
                    "safety"),
    }
