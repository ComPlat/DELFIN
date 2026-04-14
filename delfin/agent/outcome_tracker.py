"""Outcome tracking for DELFIN agent cycles.

Records task outcomes (pass/fail, cost, retries, errors) in an append-only
JSONL file.  Used by the provider profile system for adaptive routing.
"""

from __future__ import annotations

import json
import time
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path


_DEFAULT_PATH = Path.home() / ".delfin" / "outcome_history.jsonl"
_MAX_ENTRIES = 1000


@dataclass
class CycleOutcome:
    """Single cycle outcome record."""

    task: str = ""                      # User task (truncated to 200 chars)
    provider: str = ""                  # "claude" | "openai" | "kit"
    model: str = ""                     # "opus" | "sonnet" | "gpt-5.4" | ...
    mode: str = ""                      # "solo" | "quick" | "reviewed" | ...
    verdict: str = ""                   # "PASS" | "FAIL" | "PARTIAL" | ""
    cost_usd: float = 0.0
    duration_s: float = 0.0
    retries: int = 0
    denied_commands: list[str] = field(default_factory=list)
    error_type: str | None = None       # "timeout" | "crash" | "permission" | None
    task_class: str = ""                # "chemistry" | "coding" | "dashboard" | ""
    timestamp: str = ""


def append_outcome(
    outcome: CycleOutcome,
    path: Path | None = None,
) -> None:
    """Append a cycle outcome to the JSONL history file."""
    p = path or _DEFAULT_PATH
    p.parent.mkdir(parents=True, exist_ok=True)

    if not outcome.timestamp:
        outcome.timestamp = datetime.now().isoformat()

    with open(p, "a", encoding="utf-8") as f:
        f.write(json.dumps(asdict(outcome), ensure_ascii=False) + "\n")

    # Trim to max entries (keep newest)
    _trim_if_needed(p)


def load_outcomes(
    path: Path | None = None,
    max_entries: int = 100,
) -> list[CycleOutcome]:
    """Load recent outcomes from the JSONL file."""
    p = path or _DEFAULT_PATH
    if not p.exists():
        return []
    try:
        lines = p.read_text(encoding="utf-8").strip().splitlines()
        recent = lines[-max_entries:]
        results = []
        for line in recent:
            try:
                data = json.loads(line)
                results.append(CycleOutcome(**{
                    k: v for k, v in data.items()
                    if k in CycleOutcome.__dataclass_fields__
                }))
            except (json.JSONDecodeError, TypeError):
                continue
        return results
    except Exception:
        return []


def outcomes_for_provider(
    provider: str,
    path: Path | None = None,
    max_entries: int = 50,
) -> list[CycleOutcome]:
    """Load recent outcomes filtered by provider."""
    all_outcomes = load_outcomes(path, max_entries=max_entries * 3)
    return [o for o in all_outcomes if o.provider == provider][-max_entries:]


def _trim_if_needed(path: Path) -> None:
    """Keep only the newest _MAX_ENTRIES lines."""
    try:
        lines = path.read_text(encoding="utf-8").strip().splitlines()
        if len(lines) > _MAX_ENTRIES:
            trimmed = lines[-_MAX_ENTRIES:]
            path.write_text("\n".join(trimmed) + "\n", encoding="utf-8")
    except Exception:
        pass
