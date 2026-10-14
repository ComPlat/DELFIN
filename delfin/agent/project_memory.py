"""Project-memory file loader.

Implements ``DELFIN.MD`` auto-loading: walks up from the
agent's working directory and concatenates every project-memory file
it finds. Recognised filenames: ``DELFIN.MD`` / ``DELFIN.md``
(DELFIN's own project-instructions file) and ``AGENTS.md`` (the
tool-neutral cross-agent standard). Instruction files of other agent
frameworks are deliberately NOT loaded — DELFIN is standalone and
must not inherit another tool's rules. Walking is deepest-first
(closest to cwd wins on total budget). Stops at the user's home
directory or a filesystem root, whichever comes first.

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
_DEFAULT_MAX_CHARS = 8000
_PER_FILE_CAP = 6000


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
    for r in extra_roots or ():
        try:
            r = Path(r).resolve()
        except Exception:
            continue
        for d in _candidate_dirs(r):
            if d not in anchors:
                anchors.append(d)
    for d in anchors:
        for name in _RECOGNISED_NAMES:
            p = d / name
            try:
                if not p.is_file():
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
    for p in files:
        try:
            text = p.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        if not text.strip():
            continue
        if len(text) > per_file_cap:
            text = text[:per_file_cap].rstrip() + "\n[... truncated ...]"
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
