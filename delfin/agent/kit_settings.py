"""Persistent settings for the KIT-Toolbox coding agent.

Implements two-tier settings model:

- **User settings** at ``~/.delfin/settings.json`` — global per user.
- **Repo settings** at ``<repo>/.delfin/settings.json`` — project overrides.

Both are merged on load, repo wins on conflict. All keys live under the
``"kit"`` namespace so future non-KIT settings can coexist.

Schema (under ``settings["kit"]``)::

    {
      "default_mode":          "default" | "plan" | "acceptEdits" | "bypassPermissions",
      "extra_workspace_dirs":  ["/abs/path/1", ...],
      "allow_patterns":        ["^pytest", "^ruff", ...],
      "deny_patterns":         ["rm -rf /", ...]
    }

The settings system is intentionally minimal: a missing or unreadable
file is treated as empty config (no exceptions bubble up). The user
can edit the JSON by hand, or the agent updates it in response to
"always allow X" interactions through ``persist_*`` helpers.
"""

from __future__ import annotations

import json
import os
import threading
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional


USER_SETTINGS_PATH = Path("~/.delfin/settings.json").expanduser()
REPO_SETTINGS_RELPATH = Path(".delfin/settings.json")

_KIT_KEY = "kit"
_VALID_MODES = {"plan", "default", "acceptEdits", "bypassPermissions"}


def _fresh_defaults() -> dict[str, Any]:
    """Return a fully-fresh defaults dict (lists are new each call).

    Built per-call rather than as a module constant so callers can mutate
    the returned dict without poisoning future loads.
    """
    return {
        "default_mode": "default",
        "extra_workspace_dirs": [],
        "allow_patterns": [],
        "deny_patterns": [],
    }

_lock = threading.Lock()


@dataclass
class KitSettings:
    """In-memory view of the merged KIT settings."""

    default_mode: str = "default"
    extra_workspace_dirs: list[str] = field(default_factory=list)
    allow_patterns: list[str] = field(default_factory=list)
    deny_patterns: list[str] = field(default_factory=list)
    user_path: Path = USER_SETTINGS_PATH
    repo_path: Optional[Path] = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "default_mode": self.default_mode,
            "extra_workspace_dirs": list(self.extra_workspace_dirs),
            "allow_patterns": list(self.allow_patterns),
            "deny_patterns": list(self.deny_patterns),
        }


# ---------------------------------------------------------------------------
# Low-level JSON I/O
# ---------------------------------------------------------------------------

def _read_json(path: Path) -> dict[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as f:
            data = json.load(f)
        if isinstance(data, dict):
            return data
    except FileNotFoundError:
        return {}
    except (OSError, json.JSONDecodeError):
        return {}
    return {}


def _atomic_write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, sort_keys=True)
        f.write("\n")
    os.replace(tmp, path)


def _normalize_kit_block(block: Any) -> dict[str, Any]:
    """Coerce a loaded ``settings["kit"]`` block to a clean dict."""
    if not isinstance(block, dict):
        return _fresh_defaults()
    out = _fresh_defaults()
    mode = block.get("default_mode")
    if isinstance(mode, str) and mode in _VALID_MODES:
        out["default_mode"] = mode
    for key in ("extra_workspace_dirs", "allow_patterns", "deny_patterns"):
        v = block.get(key)
        if isinstance(v, list):
            out[key] = [str(x) for x in v if isinstance(x, (str, os.PathLike))]
    return out


def _merge(user: dict[str, Any],
           repo: Optional[dict[str, Any]]) -> dict[str, Any]:
    """Merge user + repo settings.

    Repo wins on scalars but only when it has an explicit value (i.e. the
    repo block was actually present on disk). Lists are unioned with repo
    entries kept first; ordering preserved, duplicates dropped.
    """
    out = _fresh_defaults()
    if repo is not None:
        out["default_mode"] = repo.get("default_mode") or user.get("default_mode") or "default"
    else:
        out["default_mode"] = user.get("default_mode") or "default"
    for key in ("extra_workspace_dirs", "allow_patterns", "deny_patterns"):
        merged: list[str] = []
        sources = (repo.get(key, []) if repo else [], user.get(key, []))
        for src in sources:
            for item in src:
                if item not in merged:
                    merged.append(item)
        out[key] = merged
    return out


# ---------------------------------------------------------------------------
# Public load / save API
# ---------------------------------------------------------------------------

def repo_settings_path(repo_dir: Optional[Path | str]) -> Optional[Path]:
    if not repo_dir:
        return None
    return Path(repo_dir).expanduser().resolve() / REPO_SETTINGS_RELPATH


def load(repo_dir: Optional[Path | str] = None,
         user_path: Optional[Path] = None) -> KitSettings:
    """Read user + repo settings and return the merged view."""
    if user_path is None:
        user_path = USER_SETTINGS_PATH
    with _lock:
        user_raw = _read_json(user_path)
        user_block = _normalize_kit_block(user_raw.get(_KIT_KEY)) \
            if _KIT_KEY in user_raw else _fresh_defaults()
        repo_path = repo_settings_path(repo_dir)
        repo_block: Optional[dict[str, Any]] = None
        if repo_path is not None and repo_path.exists():
            repo_raw = _read_json(repo_path)
            if _KIT_KEY in repo_raw:
                repo_block = _normalize_kit_block(repo_raw.get(_KIT_KEY))
        merged = _merge(user_block, repo_block)
    return KitSettings(
        default_mode=merged["default_mode"],
        extra_workspace_dirs=merged["extra_workspace_dirs"],
        allow_patterns=merged["allow_patterns"],
        deny_patterns=merged["deny_patterns"],
        user_path=user_path,
        repo_path=repo_path,
    )


def _save_block(path: Path, block: dict[str, Any]) -> None:
    """Write the kit block back to ``path`` while preserving non-kit keys."""
    with _lock:
        existing = _read_json(path)
        existing[_KIT_KEY] = block
        _atomic_write_json(path, existing)


# ---------------------------------------------------------------------------
# Conversational mutators (the agent calls these via remember_permission)
# ---------------------------------------------------------------------------

def _target_path(scope: str, repo_dir: Optional[Path | str],
                 user_path: Optional[Path] = None) -> Path:
    if user_path is None:
        user_path = USER_SETTINGS_PATH
    if scope == "repo":
        rp = repo_settings_path(repo_dir)
        if rp is None:
            raise ValueError("repo scope requires a repo_dir")
        return rp
    if scope == "user":
        return user_path
    raise ValueError(f"unknown scope: {scope!r} (use 'user' or 'repo')")


def _mutate(scope: str, repo_dir: Optional[Path | str],
            user_path: Optional[Path],
            mutator) -> KitSettings:
    if user_path is None:
        user_path = USER_SETTINGS_PATH
    path = _target_path(scope, repo_dir, user_path)
    raw = _read_json(path)
    block = _normalize_kit_block(raw.get(_KIT_KEY))
    mutator(block)
    _save_block(path, block)
    return load(repo_dir=repo_dir, user_path=user_path)


def persist_extra_dir(directory: str | os.PathLike, *,
                      scope: str = "user",
                      repo_dir: Optional[Path | str] = None,
                      user_path: Optional[Path] = None) -> KitSettings:
    """Add ``directory`` to the persisted ``extra_workspace_dirs`` list."""
    resolved = str(Path(directory).expanduser().resolve())

    def m(block: dict[str, Any]) -> None:
        if resolved not in block["extra_workspace_dirs"]:
            block["extra_workspace_dirs"].append(resolved)

    return _mutate(scope, repo_dir, user_path, m)


def remove_extra_dir(directory: str | os.PathLike, *,
                     scope: str = "user",
                     repo_dir: Optional[Path | str] = None,
                     user_path: Optional[Path] = None) -> KitSettings:
    """Remove ``directory`` from the persisted list (safe if absent)."""
    resolved = str(Path(directory).expanduser().resolve())

    def m(block: dict[str, Any]) -> None:
        block["extra_workspace_dirs"] = [
            d for d in block["extra_workspace_dirs"] if d != resolved
        ]

    return _mutate(scope, repo_dir, user_path, m)


def persist_pattern(pattern: str, *, kind: str = "allow",
                    scope: str = "user",
                    repo_dir: Optional[Path | str] = None,
                    user_path: Optional[Path] = None) -> KitSettings:
    """Append a regex ``pattern`` to allow_patterns or deny_patterns."""
    if kind not in {"allow", "deny"}:
        raise ValueError(f"kind must be 'allow' or 'deny', got {kind!r}")
    if not pattern:
        raise ValueError("pattern must be non-empty")
    key = "allow_patterns" if kind == "allow" else "deny_patterns"

    def m(block: dict[str, Any]) -> None:
        if pattern not in block[key]:
            block[key].append(pattern)

    return _mutate(scope, repo_dir, user_path, m)


def remove_pattern(pattern: str, *, kind: str = "allow",
                   scope: str = "user",
                   repo_dir: Optional[Path | str] = None,
                   user_path: Optional[Path] = None) -> KitSettings:
    if kind not in {"allow", "deny"}:
        raise ValueError(f"kind must be 'allow' or 'deny', got {kind!r}")
    key = "allow_patterns" if kind == "allow" else "deny_patterns"

    def m(block: dict[str, Any]) -> None:
        block[key] = [p for p in block[key] if p != pattern]

    return _mutate(scope, repo_dir, user_path, m)


def persist_default_mode(mode: str, *,
                         scope: str = "user",
                         repo_dir: Optional[Path | str] = None,
                         user_path: Optional[Path] = None) -> KitSettings:
    if mode not in _VALID_MODES:
        raise ValueError(f"mode must be one of {sorted(_VALID_MODES)}")

    def m(block: dict[str, Any]) -> None:
        block["default_mode"] = mode

    return _mutate(scope, repo_dir, user_path, m)


# ---------------------------------------------------------------------------
# Convenience: derive sensible bash auto-allow regex for a one-shot command
# ---------------------------------------------------------------------------

def suggest_pattern_for_command(cmd: str) -> str:
    """Translate a shell command into a re-usable auto-allow regex.

    Examples::

        'pytest -xvs tests/foo.py'   -> '^\\s*pytest\\b'
        'git status'                 -> '^\\s*git status\\b'
        'python3 -m delfin.cli x'    -> '^\\s*python3\\s+-m\\s+delfin\\.cli\\b'
    """
    parts = (cmd or "").strip().split()
    if not parts:
        return ""
    head = parts[0]
    if head in ("git", "uv", "pip", "conda", "npm", "yarn", "cargo") and len(parts) >= 2:
        verb = parts[1]
        return rf"^\s*{head}\s+{verb}\b"
    if head in ("python", "python3") and len(parts) >= 3 and parts[1] == "-m":
        mod = parts[2].replace(".", r"\.")
        return rf"^\s*{head}\s+-m\s+{mod}\b"
    return rf"^\s*{head}\b"
