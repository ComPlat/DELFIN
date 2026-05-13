"""Persistent memory store for the DELFIN Agent.

Stores facts, preferences, and project context that persist across
agent cycles and dashboard sessions.

Files live exclusively in ``~/.delfin/`` (per-user, per-machine) and
are written with 0600 permissions — they MUST NOT be committed to
the repo.
"""

from __future__ import annotations

import json
import os
import re
import time
from pathlib import Path
from typing import Any


_DEFAULT_PATH = Path.home() / ".delfin" / "agent_memory.json"
_STOPWORDS = {
    "the", "and", "for", "with", "this", "that", "from", "into", "your",
    "about", "task", "user", "agent", "project", "repo", "code", "file",
}


def _tokenize(text: str) -> set[str]:
    return {
        token
        for token in re.split(r"[^a-z0-9]+", (text or "").lower())
        if len(token) >= 3 and token not in _STOPWORDS
    }


def _read(path: Path) -> dict[str, Any]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError):
        return {"facts": [], "updated_at": 0}


def _write(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    data["updated_at"] = time.time()
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")
    try:
        os.chmod(path, 0o600)
    except OSError:
        pass


def load_memories(path: Path | None = None) -> list[dict[str, Any]]:
    """Return list of stored memory facts."""
    return _read(path or _DEFAULT_PATH).get("facts", [])


def save_memory(
    text: str,
    source: str = "user",
    path: Path | None = None,
) -> int:
    """Append a memory fact. Returns the index of the new fact."""
    p = path or _DEFAULT_PATH
    data = _read(p)
    facts = data.get("facts", [])
    facts.append({
        "text": text.strip(),
        "source": source,
        "created_at": time.time(),
    })
    data["facts"] = facts
    _write(p, data)
    return len(facts) - 1


def delete_memory(index: int, path: Path | None = None) -> bool:
    """Remove a memory by index. Returns True if removed."""
    p = path or _DEFAULT_PATH
    data = _read(p)
    facts = data.get("facts", [])
    if 0 <= index < len(facts):
        facts.pop(index)
        data["facts"] = facts
        _write(p, data)
        return True
    return False


def format_memory_context(
    path: Path | None = None,
    *,
    task_text: str = "",
    max_entries: int = 4,
) -> str:
    """Format relevant memories into a compact string for the system prompt."""
    facts = load_memories(path)
    if not facts:
        return ""

    selected: list[tuple[int, dict[str, Any]]] = []
    task_tokens = _tokenize(task_text)
    if task_tokens:
        for i, fact in enumerate(facts):
            text = str(fact.get("text", ""))
            score = len(task_tokens & _tokenize(text))
            if score > 0:
                selected.append((score * 100 + i, fact))
        selected = sorted(selected, key=lambda item: item[0], reverse=True)[:max_entries]
    else:
        selected = [
            (i, fact)
            for i, fact in enumerate(facts[-max_entries:])
        ]

    if not selected:
        return ""

    lines = []
    for _, fact in selected:
        src = fact.get("source", "?")
        lines.append(f"({src}) {fact['text']}")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Typed Claude-Code-style memory writer
# ---------------------------------------------------------------------------
#
# In addition to the flat ``~/.delfin/agent_memory.json`` store we mirror
# memories into the per-project layout that the prompt loader already reads
# back at startup (``~/.claude/projects/<slug>/memory/MEMORY.md`` plus
# typed sidecar files). This gives the user a richer, browseable, version-
# controllable memory while keeping the legacy JSON for backwards-compat.

_TYPE_LABELS = {
    "user": "User (background)",
    "feedback": "Feedback (how to work)",
    "project": "Project (current work)",
    "reference": "Reference (read only when relevant)",
}

_TYPE_PREFIX_RE = re.compile(
    r"^(user|feedback|project|reference)\s*[:\-]\s*", re.IGNORECASE
)


def parse_memory_type(text: str) -> tuple[str, str]:
    """Strip a leading ``<type>:`` prefix and return ``(memory_type, body)``.

    Defaults to ``"user"`` when no prefix is present.
    """
    m = _TYPE_PREFIX_RE.match(text or "")
    if m:
        return m.group(1).lower(), (text[m.end():]).strip()
    return "user", (text or "").strip()


def _slugify(text: str, max_len: int = 60) -> str:
    s = re.sub(r"[^a-z0-9\s\-]", "", (text or "").lower())
    s = re.sub(r"\s+", "-", s).strip("-")
    return (s or "memory")[:max_len].rstrip("-") or "memory"


def _claude_memory_dir(repo_root: Path) -> Path:
    slug = "-" + str(repo_root.resolve()).replace("/", "-").lstrip("-")
    return Path.home() / ".claude" / "projects" / slug / "memory"


def save_typed_memory(
    text: str,
    *,
    repo_root: Path | str,
    memory_type: str | None = None,
    title: str | None = None,
) -> tuple[Path, str, str]:
    """Persist a typed memory in the Claude-Code project memory layout.

    Writes:
    - ``~/.claude/projects/<slug>/memory/<type>_<kebab-slug>.md`` with
      ``name:``, ``description:`` and ``metadata.type`` frontmatter
    - Prepends a one-line pointer to that file under the matching section
      of ``MEMORY.md``, creating the section if missing.

    Returns ``(file_path, slug, memory_type)``.
    """
    body = (text or "").strip()
    if memory_type is None:
        memory_type, body = parse_memory_type(body)
    memory_type = (memory_type or "user").lower()
    if memory_type not in _TYPE_LABELS:
        memory_type = "user"

    first_line = (body.splitlines() or [""])[0].strip()
    display_title = (title or first_line)[:80].strip() or "memory"
    slug = _slugify(display_title)

    memory_dir = _claude_memory_dir(Path(repo_root))
    memory_dir.mkdir(parents=True, exist_ok=True)

    fname = f"{memory_type}_{slug}.md"
    fpath = memory_dir / fname
    if fpath.exists():
        # Don't clobber — disambiguate with a monotonic suffix.
        fname = f"{memory_type}_{slug}_{int(time.time())}.md"
        fpath = memory_dir / fname

    description = first_line[:160] or display_title
    front = (
        "---\n"
        f"name: {slug}\n"
        f"description: {description}\n"
        f"metadata:\n"
        f"  type: {memory_type}\n"
        "---\n\n"
        f"{body}\n"
    )
    fpath.write_text(front, encoding="utf-8")

    _update_memory_index(
        memory_dir,
        memory_type=memory_type,
        title=display_title,
        filename=fname,
        hook=description,
    )
    return fpath, slug, memory_type


def _update_memory_index(
    memory_dir: Path,
    *,
    memory_type: str,
    title: str,
    filename: str,
    hook: str,
) -> None:
    """Insert a one-line link under the matching ``## <Label>`` section of
    ``MEMORY.md``. Creates the index file or the section header if missing."""
    index = memory_dir / "MEMORY.md"
    label = _TYPE_LABELS[memory_type]
    section_header = f"## {label}"
    line = f"- [{title}]({filename}) — {hook}"

    if index.exists():
        content = index.read_text(encoding="utf-8")
    else:
        content = "# DELFIN Project Memory\n\n"

    if section_header in content:
        idx = content.index(section_header)
        try:
            eol = content.index("\n", idx) + 1
        except ValueError:
            eol = len(content)
        content = content[:eol] + line + "\n" + content[eol:]
    else:
        content = content.rstrip() + f"\n\n{section_header}\n{line}\n"

    index.write_text(content, encoding="utf-8")
