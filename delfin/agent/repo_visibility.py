"""Whether a repository is public, for the share path's rule set.

Input: a repository root. Output: ``"public"``, ``"private"``,
``"internal"`` or ``"unknown"``. :func:`treated_public` maps everything
but a confirmed private/internal answer to public: a share into a
repository whose visibility could not be established gets the public
rules and the public warning.

Source: ``gh repo view --json visibility`` in the repository, run only
when ``probe=True`` (the share dialog), with a 15 s timeout; the answer
is cached per repository root in ``~/.delfin/knowledge/visibility.json``
(0600) for the end-of-turn hint and the /memory label, which never run
a subprocess. No gh, no login, no network: ``unknown``.

Rejected: parsing the remote URL (a github.com remote may be private,
and a mirror may be public); asking the person (the answer is a fact
the forge knows, and a guess here decides which rules apply).
"""
from __future__ import annotations

import json
import subprocess
from pathlib import Path
from typing import Callable, Optional

PUBLIC = "public"
PRIVATE = "private"
INTERNAL = "internal"
UNKNOWN = "unknown"

#: Tests replace this; production runs gh.
_RUNNER: Optional[Callable[[Path], str]] = None


def _cache_path() -> Path:
    return Path.home() / ".delfin" / "knowledge" / "visibility.json"


def _read_cache() -> dict:
    try:
        data = json.loads(_cache_path().read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else {}
    except (OSError, ValueError):
        return {}


def _write_cache(key: str, value: str) -> None:
    try:
        from .memory_store import _atomic_write, _private_dir
        p = _cache_path()
        _private_dir(p.parent)
        data = _read_cache()
        data[key] = value
        _atomic_write(p, json.dumps(data, sort_keys=True))
    except Exception:
        pass


def _gh(repo_root: Path) -> str:
    try:
        out = subprocess.run(
            ["gh", "repo", "view", "--json", "visibility", "-q",
             ".visibility"], cwd=str(repo_root), capture_output=True,
            text=True, timeout=15)
    except (OSError, subprocess.SubprocessError):
        return UNKNOWN
    if out.returncode != 0:
        return UNKNOWN
    return (out.stdout or "").strip().lower() or UNKNOWN


def visibility(repo_root: Path | str, *, probe: bool = False) -> str:
    key = str(Path(repo_root).resolve())
    if probe:
        runner = _RUNNER or _gh
        try:
            value = str(runner(Path(key)) or UNKNOWN).strip().lower()
        except Exception:
            value = UNKNOWN
        if value not in (PUBLIC, PRIVATE, INTERNAL):
            value = UNKNOWN
        _write_cache(key, value)
        return value
    value = str(_read_cache().get(key) or UNKNOWN)
    return value if value in (PUBLIC, PRIVATE, INTERNAL) else UNKNOWN


def treated_public(value: str) -> bool:
    return value not in (PRIVATE, INTERNAL)


def label(value: str) -> str:
    return "public" if treated_public(value) else value
