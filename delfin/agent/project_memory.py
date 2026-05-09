"""Project-memory file loader.

Mirrors Claude Code's ``CLAUDE.md`` auto-loading: walks up from the
agent's working directory and concatenates every project-memory file
it finds. Recognised filenames: ``CLAUDE.md``, ``AGENTS.md``,
``DELFIN.md``. Walking is deepest-first (closest to cwd wins on
total budget). Stops at the user's home directory or a filesystem
root, whichever comes first.

Each found file is rendered as

    --- Project memory: <relative-path> ---
    <content>

so the agent can tell which rules came from which scope.

The loader is safe to call cheaply on every turn — IO is small and
files are read at most once per call.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

_RECOGNISED_NAMES = ("CLAUDE.md", "AGENTS.md", "DELFIN.md")
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
    seen: set[Path] = set()
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
                if p.is_file() and p not in seen:
                    seen.add(p)
                    out.append(p)
            except OSError:
                continue
    return out


def load_project_memory(
    cwd: Path | str | None = None,
    *,
    extra_roots: Iterable[Path] | None = None,
    max_chars: int = _DEFAULT_MAX_CHARS,
    per_file_cap: int = _PER_FILE_CAP,
) -> str:
    """Concatenate discovered project-memory files into a single string.

    Returns "" if no files are found. Truncates to ``max_chars`` total
    and ``per_file_cap`` per file so a misplaced wall-of-text doesn't
    blow the system-prompt budget.
    """
    files = discover_project_memory_files(cwd, extra_roots=extra_roots)
    if not files:
        return ""
    parts: list[str] = []
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
        used += len(block)
    return "\n".join(parts).strip()


__all__ = ["discover_project_memory_files", "load_project_memory"]
