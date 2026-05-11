"""Outcome tracking for DELFIN agent cycles.

Records task outcomes (pass/fail, cost, retries, errors) in an append-only
JSONL file.  Used by the provider profile system for adaptive routing.

The history file is per-user/per-machine (``~/.delfin/outcome_history.jsonl``)
and is written with 0600 permissions.  It MUST NOT be committed to the repo.
"""

from __future__ import annotations

import json
import os
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
    cost_usd: float = 0.0               # Cumulative engine.cost_usd at outcome time
    cost_usd_delta: float = 0.0         # Δ for THIS cycle (cost_usd - prev cycle's
                                        # cost_usd). The honest per-turn spend; sum
                                        # across cycles equals real session total.
                                        # 0.0 on legacy entries written before A6.
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

    new_file = not p.exists()
    with open(p, "a", encoding="utf-8") as f:
        f.write(json.dumps(asdict(outcome), ensure_ascii=False) + "\n")
    if new_file:
        try:
            os.chmod(p, 0o600)
        except OSError:
            pass

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


def update_last_outcome(
    *,
    retries: int | None = None,
    verdict: str | None = None,
    cost_usd: float | None = None,
    cost_usd_delta: float | None = None,
    duration_s: float | None = None,
    denied_commands: list[str] | None = None,
    error_type: str | None = None,
    path: Path | None = None,
) -> bool:
    """Mutate the most recent JSONL outcome in-place.

    Used by ``/retry`` (B3): instead of appending a fresh outcome for
    a re-attempted task, bump the existing entry's ``retries`` counter
    and optionally swap the verdict to reflect the retry result.

    Only the keyword args that are passed (non-None) get rewritten.
    Returns True on success, False if the file is missing/empty/unreadable.
    """
    p = path or _DEFAULT_PATH
    if not p.exists():
        return False
    try:
        text = p.read_text(encoding="utf-8")
    except Exception:
        return False
    lines = text.splitlines()
    # Walk from the end to find the last non-empty line (some writers
    # leave a trailing newline that produces a blank tail entry).
    last_idx = -1
    for i in range(len(lines) - 1, -1, -1):
        if lines[i].strip():
            last_idx = i
            break
    if last_idx < 0:
        return False
    try:
        record = json.loads(lines[last_idx])
    except json.JSONDecodeError:
        return False
    if not isinstance(record, dict):
        return False

    updates: dict = {}
    if retries is not None:
        updates["retries"] = int(retries)
    if verdict is not None:
        updates["verdict"] = str(verdict)
    if cost_usd is not None:
        updates["cost_usd"] = float(cost_usd)
    if cost_usd_delta is not None:
        updates["cost_usd_delta"] = float(cost_usd_delta)
    if duration_s is not None:
        updates["duration_s"] = float(duration_s)
    if denied_commands is not None:
        updates["denied_commands"] = list(denied_commands)
    if error_type is not None:
        updates["error_type"] = str(error_type)

    if not updates:
        return False

    record.update(updates)
    lines[last_idx] = json.dumps(record, ensure_ascii=False)
    try:
        p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    except Exception:
        return False
    return True
