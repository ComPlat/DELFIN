"""What a backend has actually accepted, per model: a floor under the window.

The static window table is a hand table, and it understates. For
``kit.deepseek-v4-flash`` it says 131,072 while the endpoint accepted a
single request of 166,379 tokens (field metrics, 2026-10-08); for
``kit.glm-5.3-flash`` it had no entry at all and the heuristic gave
32,768 against requests of 85,535 that went through. Every mechanism that
protects a long session -- the sliding trim, compaction, the fresh start
-- fires at a share of the window, so a window that is too small is not
conservative: it throws context away that the model would have held.

The live ``/v1/models`` probe is the first source and reads the true
``max_model_len`` when it can authenticate; this is the second, for
hosts where it cannot. A request the backend answered with a token count
is proof the window is at least that big. The floor only ever raises the
window and never claims more than was accepted.

Persisted per model in the user's state directory, so a floor learned in
one session holds for the next. Never raises: a floor that cannot be
written is simply not learned.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Optional

#: Where the floors live. A module attribute, not a computed path, so the
#: test fixture that points every user-state sink into a test's own
#: directory can redirect this one too (tests/conftest.py, state_paths).
_PATH: Path = Path.home() / ".delfin" / "observed_windows.json"


def _path() -> Optional[Path]:
    try:
        return Path(_PATH)
    except Exception:
        return None


def _load() -> dict:
    p = _path()
    if p is None or not p.exists():
        return {}
    try:
        data = json.loads(p.read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else {}
    except Exception:
        return {}


def floor(model: str) -> int:
    """The largest input the backend has accepted for *model*, or 0."""
    model = str(model or "").strip()
    if not model:
        return 0
    try:
        return int(_load().get(model) or 0)
    except Exception:
        return 0


def note(model: str, input_tokens: int) -> None:
    """Record that *model* accepted a request of *input_tokens*.

    Writes only when the number is a new maximum, so the per-request
    cost is one dictionary lookup almost always. Owner-only file.
    """
    model = str(model or "").strip()
    try:
        n = int(input_tokens or 0)
    except (TypeError, ValueError):
        return
    if not model or n <= 0:
        return
    try:
        data = _load()
        if n <= int(data.get(model) or 0):
            return
        data[model] = n
        p = _path()
        if p is None:
            return
        p.parent.mkdir(parents=True, exist_ok=True)
        tmp = p.with_name(p.name + ".tmp")
        tmp.write_text(json.dumps(data, indent=1, sort_keys=True), encoding="utf-8")
        try:
            os.chmod(tmp, 0o600)
        except OSError:
            pass
        tmp.replace(p)
    except Exception:
        pass


def at_least(model: str, window: int) -> int:
    """*window* raised to the observed floor for *model*, never lowered."""
    try:
        return max(int(window or 0), floor(model))
    except Exception:
        return int(window or 0)
