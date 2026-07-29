"""Storage helpers for the toy bookmark keeper.

Plain standard library on purpose: this workspace stands in for an ordinary
user project, so nothing here may depend on the surrounding checkout.
"""

from __future__ import annotations

import json
from pathlib import Path

DEFAULT_STORE = Path(__file__).resolve().parent / "bookmarks.json"


def load_bookmarks(path: str | Path | None = None) -> list[dict]:
    """Return the bookmarks stored in ``path`` (defaults to the JSON file)."""
    store = Path(path) if path else DEFAULT_STORE
    with store.open(encoding="utf-8") as handle:
        data = json.load(handle)
    return [item for item in data if isinstance(item, dict)]


def collect_tags(bookmarks: list[dict]) -> list[str]:
    """Return every tag used by ``bookmarks``, de-duplicated and sorted."""
    seen: set[str] = set()
    for entry in bookmarks:
        seen.update(entry.get("tags") or [])
    return sorted(seen)
