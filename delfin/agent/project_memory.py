"""Project-instruction file loader: the one loader every entry path uses.

Walks up from the agent's working directory and concatenates every
project-memory file it finds. Recognised filenames: ``DELFIN.MD`` /
``DELFIN.md`` (DELFIN's own project-instructions file) and ``AGENTS.md``
(the tool-neutral cross-agent standard), from the workspace and its
parents; ``CLAUDE.md`` and ``CONTRIBUTING.md`` from a workspace root
only. CLAUDE.md is read because a repository written for another agent
states its conventions there (build/test commands, layout, rules), and
a project's conventions are the project's, whichever agent they were
written for; like every file here it is fenced as repository text and
can direct nothing. Walking is deepest-first (closest to cwd wins on
total budget). Stops at the user's home directory or a filesystem root,
whichever comes first. A file that is a link out of its own directory
is skipped.

Callers: ``turn_context.project_instructions`` (the terminal session,
``-p``/run and the dashboard all compose a turn through it) and the
sub-agent spawn (``subagents._inherited_context``, smaller cap).

Each found file is rendered as

    --- Project memory: <relative-path> ---
    <content>

so the agent can tell which rules came from which scope.

The loader is safe to call cheaply on every turn — IO is small and
files are read at most once per call.
"""

from __future__ import annotations

import hashlib
import hmac
import secrets
from pathlib import Path
from typing import Iterable

_RECOGNISED_NAMES = ("DELFIN.MD", "DELFIN.md", "AGENTS.md")
#: Read in a workspace root only (the session workspace and each extra
#: workspace directory), never in a parent: a CONTRIBUTING.md or CLAUDE.md
#: one level up belongs to some other project. CONTRIBUTING.md last, so
#: the total cap trims it first.
_WORKSPACE_ONLY_NAMES = ("CLAUDE.md", "CONTRIBUTING.md")
_DEFAULT_MAX_CHARS = 8000
_PER_FILE_CAP = 6000
#: Per-name caps below _PER_FILE_CAP. A contributor guide is mostly for
#: humans (CLA, issue templates, code of conduct); its head carries the
#: branch and commit rules.
_NAME_CAPS = {"CONTRIBUTING.md": 3000}


def _candidate_dirs(start: Path, *, stop_at_home: bool = True) -> list[Path]:
    """Yield ``start``, its parents, up to home or filesystem root."""
    out: list[Path] = []
    home = Path.home().resolve()
    cur = start.resolve()
    while True:
        out.append(cur)
        parent = cur.parent
        if parent == cur:
            break
        if stop_at_home and cur == home:
            break
        cur = parent
    return out


def discover_project_memory_files(
    cwd: Path | str | None = None,
    *,
    extra_roots: Iterable[Path] | None = None,
) -> list[Path]:
    """Return paths of project-memory files from cwd upward.

    Closest-to-cwd files come first. Each path appears at most once.
    ``extra_roots`` lets callers add additional search anchors (e.g.
    KitToolPermissions.extra_workspace_dirs) so multi-repo sessions
    pick up the right rules.
    """
    cwd = Path(cwd or Path.cwd()).resolve()
    seen: set[tuple[int, int]] = set()
    out: list[Path] = []
    anchors: list[Path] = list(_candidate_dirs(cwd))
    roots: list[Path] = [cwd]
    for r in extra_roots or ():
        try:
            r = Path(r).resolve()
        except Exception:
            continue
        roots.append(r)
        for d in _candidate_dirs(r):
            if d not in anchors:
                anchors.append(d)
    for d in anchors:
        names = _RECOGNISED_NAMES + (
            _WORKSPACE_ONLY_NAMES if d in roots else ())
        for name in names:
            p = d / name
            try:
                if not p.is_file():
                    continue
                # A file that is a link out of its directory is not that
                # directory's instructions -- and AGENTS.md pointing at a
                # key file would put the key into the prompt.
                try:
                    p.resolve().relative_to(d)
                except ValueError:
                    continue
                # Dedupe by inode so DELFIN.MD/DELFIN.md don't load twice
                # on a case-insensitive filesystem.
                st = p.stat()
                ident = (st.st_dev, st.st_ino)
                if ident in seen:
                    continue
                seen.add(ident)
                out.append(p)
            except OSError:
                continue
    return out


#: Per-process key for the fence nonce. The nonce has to be (a) stable
#: between prompt builds, or the cached prompt prefix changes every turn,
#: and (b) unknowable to the file's author, or the file can close its own
#: fence. HMAC(process key, content) is both; a content hash alone is (a)
#: without (b) for any reader who can hash. A new process draws a new
#: key, so the fence (and the cache) changes once per process.
_FENCE_KEY = secrets.token_bytes(32)


def _fence_nonce(body: str) -> str:
    return hmac.new(_FENCE_KEY, body.encode("utf-8", "surrogateescape"),
                    hashlib.sha256).hexdigest()[:16]


def fence_project_text(body: str, sources: Iterable[str]) -> str:
    """Wrap repository instruction text as reference, not as the user.

    Input: the concatenated file blocks and the names they came from.
    Output: one ``untrusted.wrap_reference`` fence, byte-stable within a
    process for the same content. "" for empty input.
    """
    if not (body or "").strip():
        return ""
    from . import untrusted
    names = ", ".join(dict.fromkeys(sources)) or "project files"
    return untrusted.wrap_reference(
        f"project instructions: {names}", body.strip(),
        nonce=_fence_nonce(body.strip()))


def load_project_memory(
    cwd: Path | str | None = None,
    *,
    extra_roots: Iterable[Path] | None = None,
    max_chars: int = _DEFAULT_MAX_CHARS,
    per_file_cap: int = _PER_FILE_CAP,
) -> str:
    """Concatenate discovered project-memory files into one fenced block.

    Returns "" if no files are found. Truncates to ``max_chars`` total
    and ``per_file_cap`` per file so a misplaced wall-of-text doesn't
    blow the system-prompt budget. The result is fenced by
    :func:`fence_project_text`: the files are written by whoever wrote
    the repository, so they reach the model as reference text, never as
    a message from the user, and a reader that cuts at the untrusted
    marker (the action classifier) never takes them for intent.
    """
    files = discover_project_memory_files(cwd, extra_roots=extra_roots)
    if not files:
        return ""
    parts: list[str] = []
    names: list[str] = []
    used = 0
    home = str(Path.home())
    base = Path(cwd or Path.cwd()).resolve()
    for p in files:
        cap = min(per_file_cap, _NAME_CAPS.get(p.name, per_file_cap))
        try:
            # Read no more than the cap needs: a repository decides the
            # file's size.
            with open(p, "r", encoding="utf-8", errors="replace") as fh:
                text = fh.read(cap + 1)
        except OSError:
            continue
        if not text.strip():
            continue
        if len(text) > cap:
            text = text[:cap].rstrip() + "\n[... truncated ...]"
        # Relative to the workspace where it lies inside it: the absolute
        # path is the same bytes in every block and says nothing more.
        try:
            rel = p.relative_to(base).as_posix()
        except ValueError:
            rel = str(p)
            if rel.startswith(home):
                rel = "~" + rel[len(home):]
        header = f"--- Project memory: {rel} ---"
        block = f"{header}\n{text.strip()}\n"
        if used + len(block) > max_chars and parts:
            break
        parts.append(block)
        names.append(p.name)
        used += len(block)
    return fence_project_text("\n".join(parts), names)


__all__ = ["discover_project_memory_files", "fence_project_text",
           "load_project_memory"]
