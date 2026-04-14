"""Persistent memory store for the DELFIN Agent.

Stores facts, preferences, and project context that persist across
agent cycles and dashboard sessions.
"""

from __future__ import annotations

import json
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
