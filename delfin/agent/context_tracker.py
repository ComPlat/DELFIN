"""Track which injected context sections the agent actually references.

Records per-section hit rates so the prompt loader can skip sections that
are consistently ignored — reducing token cost without losing quality.

Storage: ``~/.delfin/context_usage.jsonl`` (append-only, max 500 entries).
"""

from __future__ import annotations

import json
import re
import threading
import time
from collections import Counter
from pathlib import Path
from typing import Any


_DEFAULT_PATH = Path.home() / ".delfin" / "context_usage.jsonl"
_MAX_ENTRIES = 500

# Minimum data points before we trust the hit rate enough to skip a section.
_MIN_SAMPLES = 8
# Skip threshold: if a section is referenced less than this fraction, skip it.
_SKIP_THRESHOLD = 0.15

# Section-specific detection keywords.  If the response contains any of these
# terms (case-insensitive) we count it as a "hit" for that section.
_SECTION_MARKERS: dict[str, list[str]] = {
    "playbook": [
        "playbook", "step 1", "step 2", "key_invariant",
        "key invariant", "according to the playbook",
    ],
    "repo_map": [
        "repo map", "repo_map", "delfin/",
        "tests/test_", "build_up_complex.py", "smiles_converter.py",
    ],
    "briefing": [
        "briefing", "pass rate", "outcome history",
        "denied commands", "task class:",
    ],
    "profile": [
        "provider profile", "success rate", "common failure",
        "provider failures", "denied pattern",
    ],
    "memory": [
        "project memory", "memory fact", "remembered",
    ],
    "prior_outputs": [
        "prior output", "session_manager plan",
        "critic review", "builder output",
    ],
}


def _detect_references(
    response_text: str,
    sections_injected: list[str],
) -> dict[str, bool]:
    """Detect which sections the response actually references."""
    lower = (response_text or "").lower()
    hits: dict[str, bool] = {}
    for section in sections_injected:
        markers = _SECTION_MARKERS.get(section, [])
        hits[section] = any(m in lower for m in markers)
    return hits


class ContextUsageTracker:
    """Thread-safe tracker for context section usage."""

    def __init__(self, path: Path | None = None):
        self._path = path or _DEFAULT_PATH
        self._lock = threading.Lock()
        self._cache: list[dict[str, Any]] | None = None

    def record_usage(
        self,
        sections_injected: list[str],
        response_text: str,
        role_id: str = "",
        provider: str = "",
    ) -> dict[str, bool]:
        """Record which sections were referenced in the response.

        Returns the hit map for transparency logging.
        """
        hits = _detect_references(response_text, sections_injected)
        entry = {
            "ts": time.time(),
            "role": role_id,
            "provider": provider,
            "injected": sections_injected,
            "hits": {k: v for k, v in hits.items()},
        }
        with self._lock:
            self._append(entry)
            self._cache = None  # invalidate
        return hits

    def get_hit_rates(
        self,
        role_id: str = "",
        provider: str = "",
    ) -> dict[str, float]:
        """Return hit rates per section.  Optionally filtered by role/provider."""
        entries = self._load()
        counts: Counter[str] = Counter()
        totals: Counter[str] = Counter()
        for e in entries:
            if role_id and e.get("role") != role_id:
                continue
            if provider and e.get("provider") != provider:
                continue
            for section in e.get("injected", []):
                totals[section] += 1
                if e.get("hits", {}).get(section, False):
                    counts[section] += 1
        return {
            section: counts[section] / totals[section]
            for section in totals
            if totals[section] > 0
        }

    def should_skip(
        self,
        section_name: str,
        role_id: str = "",
        provider: str = "",
    ) -> bool:
        """True if the section has enough data and a low hit rate."""
        entries = self._load()
        injected_count = 0
        hit_count = 0
        for e in entries:
            if role_id and e.get("role") != role_id:
                continue
            if provider and e.get("provider") != provider:
                continue
            if section_name in e.get("injected", []):
                injected_count += 1
                if e.get("hits", {}).get(section_name, False):
                    hit_count += 1
        if injected_count < _MIN_SAMPLES:
            return False
        return (hit_count / injected_count) < _SKIP_THRESHOLD

    # -- storage helpers ---------------------------------------------------

    def _load(self) -> list[dict[str, Any]]:
        with self._lock:
            if self._cache is not None:
                return self._cache
        entries = _read_jsonl(self._path)
        with self._lock:
            self._cache = entries
        return entries

    def _append(self, entry: dict[str, Any]) -> None:
        self._path.parent.mkdir(parents=True, exist_ok=True)
        with open(self._path, "a", encoding="utf-8") as f:
            f.write(json.dumps(entry, ensure_ascii=False) + "\n")
        _trim_if_needed(self._path)


def _read_jsonl(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    try:
        lines = path.read_text(encoding="utf-8").strip().splitlines()
        entries = []
        for line in lines[-_MAX_ENTRIES:]:
            try:
                entries.append(json.loads(line))
            except json.JSONDecodeError:
                continue
        return entries
    except Exception:
        return []


def _trim_if_needed(path: Path) -> None:
    try:
        lines = path.read_text(encoding="utf-8").strip().splitlines()
        if len(lines) > _MAX_ENTRIES:
            trimmed = lines[-_MAX_ENTRIES:]
            path.write_text("\n".join(trimmed) + "\n", encoding="utf-8")
    except Exception:
        pass
