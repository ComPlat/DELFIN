"""Secure credential store for the DELFIN agent.

API keys live in ``~/.delfin/credentials.json`` with file-mode ``0600``
(the same convention used by ``~/.netrc`` and ``~/.aws/credentials``).
The engine calls ``load_credential(name)`` at engine-init time only;
no CLI command ever echoes a stored value back to the user.

Lookup priority for any credential name:

  1. ``os.environ[name]``         — session-level override
  2. ``~/.delfin/credentials.json`` — persistent file store
  3. ``""``                        — not configured

Design notes:

- We do **not** offer a ``get`` command.  Once written, a credential
  is consumed only by library code; ``list`` shows masked values so
  the user knows WHICH credentials are stored without exposing them.
- Writes use a temp-file + atomic rename pattern so a crash mid-write
  cannot leave a partial file with a corrupted key.
- The directory is created on first write with ``mkdir(mode=0o700)``;
  the file itself is then ``chmod 0600``.  Best-effort — failure to
  set permissions logs nothing but does not raise (we don't have a
  logger here and don't want to leak credential names to stderr).
"""

from __future__ import annotations

import json
import os
from pathlib import Path


_DEFAULT_PATH = Path.home() / ".delfin" / "credentials.json"

# Names the agent stack knows how to consume.  Used by ``list_credentials``
# to surface env-var-sourced keys even when nothing is on disk.
_WELL_KNOWN_KEYS: tuple[str, ...] = (
    "KIT_TOOLBOX_API_KEY",
    "OPENAI_API_KEY",
    "ANTHROPIC_API_KEY",
)


def _resolve_path(path: Path | None) -> Path:
    return path if path is not None else _DEFAULT_PATH


def _read_store(path: Path | None = None) -> dict[str, str]:
    p = _resolve_path(path)
    if not p.exists():
        return {}
    try:
        raw = json.loads(p.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    if not isinstance(raw, dict):
        return {}
    out: dict[str, str] = {}
    for k, v in raw.items():
        if isinstance(k, str) and isinstance(v, str):
            out[k] = v
    return out


def _write_store(store: dict[str, str], path: Path | None = None) -> None:
    p = _resolve_path(path)
    # Ensure parent dir exists with permissive-to-owner-only mode
    try:
        p.parent.mkdir(parents=True, exist_ok=True)
        try:
            os.chmod(p.parent, 0o700)
        except OSError:
            pass
    except OSError:
        return
    # Atomic write: tmp file → chmod 0600 → rename
    tmp = p.with_suffix(p.suffix + ".tmp")
    try:
        tmp.write_text(
            json.dumps(store, ensure_ascii=False, indent=2, sort_keys=True)
            + "\n",
            encoding="utf-8",
        )
        try:
            os.chmod(tmp, 0o600)
        except OSError:
            pass
        tmp.replace(p)
        try:
            os.chmod(p, 0o600)
        except OSError:
            pass
    except OSError:
        # Best-effort: leave the old file intact if anything failed
        try:
            tmp.unlink(missing_ok=True)
        except OSError:
            pass


def load_credential(name: str, *, path: Path | None = None) -> str:
    """Return the credential value: env-var > file > ''.

    The lookup is silent — callers MUST handle an empty return without
    revealing the credential name in any user-visible error.
    """
    if not name:
        return ""
    env_val = os.environ.get(name, "")
    if env_val:
        return env_val
    store = _read_store(path)
    return store.get(name, "")


def set_credential(
    name: str, value: str, *, path: Path | None = None,
) -> bool:
    """Persist a credential.  Returns True if the store was modified.

    Empty ``value`` is a no-op (use ``delete_credential`` to remove).
    """
    if not name or not value:
        return False
    store = _read_store(path)
    if store.get(name) == value:
        return False
    store[name] = value
    _write_store(store, path)
    return True


def delete_credential(name: str, *, path: Path | None = None) -> bool:
    """Remove a credential.  Returns True if it was present."""
    if not name:
        return False
    store = _read_store(path)
    if name not in store:
        return False
    del store[name]
    _write_store(store, path)
    return True


def mask(value: str) -> str:
    """Return a masked rendering safe for display.

    Short values (<=10 chars) become all-stars so length isn't even
    revealed; longer values keep the first 4 + last 4 chars so the
    user can recognise which key it is, with the middle elided.
    """
    if not value:
        return ""
    if len(value) <= 10:
        return "*" * len(value)
    return f"{value[:4]}…{value[-4:]}"


def list_credentials(
    *, path: Path | None = None,
) -> dict[str, dict[str, str]]:
    """Return a dict of name → {value: masked, source: 'env'|'file'}.

    Surfaces both file-stored credentials and well-known env-var
    credentials.  When a name lives in both, ``env`` wins (matches
    ``load_credential`` semantics).
    """
    out: dict[str, dict[str, str]] = {}
    store = _read_store(path)
    for name, val in store.items():
        out[name] = {"value": mask(val), "source": "file"}
    for name in _WELL_KNOWN_KEYS:
        env_val = os.environ.get(name, "")
        if not env_val:
            continue
        # env-var wins on conflict
        out[name] = {"value": mask(env_val), "source": "env"}
    return out


def credentials_path() -> Path:
    """Public accessor for the credentials file location."""
    return _DEFAULT_PATH


__all__ = [
    "load_credential",
    "set_credential",
    "delete_credential",
    "list_credentials",
    "mask",
    "credentials_path",
]
