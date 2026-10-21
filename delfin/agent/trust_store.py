"""The account's trust store: a key and pins that no agent route reaches.

Input: a purpose name (``session_messages``, ``audit_log``) or a record
name. Output: a 32-byte key derived for that purpose, an HMAC under it,
or a JSON record kept beside the key. Location: ``~/.delfin/trust_store``
(0700), resolved per call (``store_dir``; the suite and a benchmark
attempt redirect it like every user-state sink).

What is kept here:

* ``account.key`` -- 32 bytes from ``secrets``, 0600, created once with
  ``O_CREAT | O_EXCL`` by the first DELFIN process that needs it. Every
  key DELFIN uses is derived from it per purpose
  (``HMAC-SHA256(account key, "delfin:" + purpose)``), so a key that
  leaks for one purpose does not sign for another;
* ``mcp_pins.json`` -- the hash each MCP tool's description had when it
  was first seen (``mcp_client.pin_verdict``);
* ``audit-anchor-*.json`` -- where each audit log ended at its last write
  (``audit_log``).

Who can reach it. It is a protected location on every agent route: the
path deny globs of the file tools and of the bash read scan
(``api_client._DEFAULT_PATH_DENY_GLOBS``, by the name ``trust_store``),
the CLI backend's deny rules, and the cages (``sandbox._HOME_SECRET_DIRS``
and ``sandbox.runtime_masked_dirs``: a bubblewrap tmpfs after the binds, a
Seatbelt hide, a Landlock hide unless a write root contains it). A process
that carries the agent-session marker (``approval_guard.ENV_MARKER``) gets
no key from here and creates none, so DELFIN code an agent's command
imports cannot sign as the account.

A key file is used only when it is a regular file of this account, not a
link, with no group or other permission bits and exactly 32 bytes; any
other file in its place is refused (recorded once per process as
``trust_store_refused``) and the callers fall back to their unkeyed form,
which they report as unverified.

Limits. A same-account process outside every cage that reads the file
directly has the key: the account is outside the boundary
(docs/SECURITY_MODEL.md, "Against whom"). Where no isolation holds, the
protection is the read gates, the same level as the approval-code
directory. Deleting the key makes every earlier signature unverifiable;
a new key is created at the next use and the earlier records then read
as not verified.

Rejected: the OS keyring (not present on most cluster nodes, and a
keyring item is readable by the same interpreter without a prompt on
macOS); a key per purpose in separate files (more files to protect, same
boundary); a key in the environment (every child process inherits it).
Cost: one ``open``/``read`` of 32 bytes per process and purpose (cached),
one HMAC-SHA256 per signed or verified item (about 1 us for 300 bytes).
"""

from __future__ import annotations

import hashlib
import hmac
import json
import os
import secrets
import stat as _stat
import threading
from pathlib import Path
from typing import Any, Optional

#: The directory name, also the name the path deny globs match.
STORE_NAME = "trust_store"
KEY_FILE = "account.key"
KEY_BYTES = 32

_NOFOLLOW = getattr(os, "O_NOFOLLOW", 0)
_LOCK = threading.Lock()
#: {(key path, purpose): derived key}; keyed by path so a redirected store
#: (a test, a benchmark attempt) never reuses another store's key.
_DERIVED: "dict[tuple[str, str], bytes]" = {}
_REFUSAL_NOTED = False


def home_store_dir() -> Path:
    """``~/.delfin/trust_store``, never redirected: the canonical location
    the cages mask as well."""
    return Path.home() / ".delfin" / STORE_NAME


def store_dir() -> Path:
    """The store directory. Resolved per call (tests move HOME), and a
    user-state resolver (``state_paths.USER_STATE_RESOLVERS``)."""
    return home_store_dir()


def ensure_store_dir() -> Optional[Path]:
    """Create the directory 0700 so a cage built now masks it. None when
    it cannot be made or is not a private directory of this account."""
    d = store_dir()
    try:
        d.mkdir(parents=True, exist_ok=True, mode=0o700)
        info = os.lstat(d)
        if not _stat.S_ISDIR(info.st_mode) or (
                hasattr(os, "getuid") and info.st_uid != os.getuid()):
            return None
        if _stat.S_IMODE(info.st_mode) != 0o700:
            os.chmod(d, 0o700)
        return d
    except OSError:
        return None


def _inside_agent_session() -> bool:
    try:
        from .approval_guard import inside_agent_session
        return inside_agent_session()
    except Exception:
        return bool(os.environ.get("DELFIN_AGENT_SESSION"))


def _note_refusal(detail: str) -> None:
    global _REFUSAL_NOTED
    if _REFUSAL_NOTED:
        return
    _REFUSAL_NOTED = True
    try:
        from . import security_events
        security_events.record("trust_store_refused", "trust_store",
                               detail[:300], blocked=True)
    except Exception:
        pass


def _read_key(path: Path) -> Optional[bytes]:
    """The key in ``path``, or None (absent). Raises ValueError for a file
    that is there and must not be used."""
    try:
        fd = os.open(str(path), os.O_RDONLY | _NOFOLLOW)
    except FileNotFoundError:
        return None
    except OSError as exc:
        raise ValueError(f"{path} cannot be opened: {exc}") from exc
    try:
        info = os.fstat(fd)
        if (not _stat.S_ISREG(info.st_mode) or info.st_nlink != 1
                or (hasattr(os, "getuid") and info.st_uid != os.getuid())
                or info.st_mode & 0o077):
            raise ValueError(f"{path} is not a private regular file of "
                             "this account")
        data = os.read(fd, KEY_BYTES + 1)
    finally:
        os.close(fd)
    if len(data) != KEY_BYTES:
        raise ValueError(f"{path} does not hold a {KEY_BYTES}-byte key")
    return data


def _create_key(path: Path) -> Optional[bytes]:
    """Write a new key with O_EXCL; the loser of a race reads the winner's."""
    data = secrets.token_bytes(KEY_BYTES)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.partial")
    try:
        fd = os.open(str(tmp), os.O_WRONLY | os.O_CREAT | os.O_EXCL
                     | _NOFOLLOW, 0o600)
    except OSError:
        return None
    try:
        view = memoryview(data)
        while view:
            view = view[os.write(fd, view):]
        os.fsync(fd)
    finally:
        os.close(fd)
    try:
        os.link(str(tmp), str(path))      # fails if a key appeared meanwhile
    except FileExistsError:
        pass
    except OSError:
        try:
            os.unlink(str(tmp))
        except OSError:
            pass
        return None
    try:
        os.unlink(str(tmp))
    except OSError:
        pass
    return _read_key(path)


def account_key(*, create: bool = True) -> Optional[bytes]:
    """The account key, or None: inside an agent session, where the store
    cannot be made, or where the file in its place is refused."""
    if _inside_agent_session():
        return None
    d = ensure_store_dir() if create else store_dir()
    if d is None:
        _note_refusal(f"{store_dir()} is not a private directory of this "
                      "account; nothing is signed with the account key")
        return None
    path = d / KEY_FILE
    try:
        key = _read_key(path)
        if key is None and create:
            key = _create_key(path)
        return key
    except ValueError as exc:
        _note_refusal(f"{exc}; nothing is signed with the account key")
        return None


def derived_key(purpose: str, *, create: bool = True) -> Optional[bytes]:
    """The key for ``purpose``, derived from the account key, or None."""
    if _inside_agent_session():
        return None
    cache = (str(store_dir() / KEY_FILE), str(purpose))
    with _LOCK:
        got = _DERIVED.get(cache)
    if got is not None:
        return got
    base = account_key(create=create)
    if base is None:
        return None
    got = hmac.new(base, b"delfin:" + str(purpose).encode("utf-8"),
                   hashlib.sha256).digest()
    with _LOCK:
        _DERIVED[cache] = got
    return got


def mac(purpose: str, data: bytes, *, create: bool = True) -> str:
    """HMAC-SHA256 of ``data`` under the key for ``purpose`` (hex), or ""
    when there is no key."""
    key = derived_key(purpose, create=create)
    if key is None:
        return ""
    return hmac.new(key, data, hashlib.sha256).hexdigest()


def check_mac(purpose: str, data: bytes, tag: Any, *,
              key: Optional[bytes] = None) -> Optional[bool]:
    """True when ``tag`` is the HMAC of ``data`` for ``purpose``, False when
    it is not, None when there is no key to check it with."""
    if key is None:
        key = derived_key(purpose, create=False)
    if key is None:
        return None
    want = hmac.new(key, data, hashlib.sha256).hexdigest()
    return hmac.compare_digest(want, str(tag or ""))


def forget_cached_keys() -> None:
    """Drop the derived keys this process holds (tests; a deleted key)."""
    with _LOCK:
        _DERIVED.clear()


# -- records beside the key ---------------------------------------------------

def _record_path(name: str) -> Optional[Path]:
    safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in name)[:120]
    if not safe or safe.startswith(".") or safe == KEY_FILE:
        return None
    return store_dir() / safe


def read_record(name: str) -> Optional[dict]:
    """A JSON record of the store, or None (absent, unreadable, not a
    private regular file of this account, not an object)."""
    path = _record_path(name)
    if path is None:
        return None
    try:
        fd = os.open(str(path), os.O_RDONLY | _NOFOLLOW)
    except OSError:
        return None
    try:
        info = os.fstat(fd)
        if (not _stat.S_ISREG(info.st_mode)
                or (hasattr(os, "getuid") and info.st_uid != os.getuid())
                or info.st_mode & 0o077 or info.st_size > 8 << 20):
            return None
        chunks = []
        while True:
            chunk = os.read(fd, 1 << 20)
            if not chunk:
                break
            chunks.append(chunk)
    finally:
        os.close(fd)
    try:
        data = json.loads(b"".join(chunks).decode("utf-8"))
    except ValueError:
        return None
    return data if isinstance(data, dict) else None


def write_record(name: str, data: dict) -> bool:
    """Write a JSON record whole (temp file 0600, then rename). Refused
    inside an agent session, like the key. True when written."""
    if _inside_agent_session():
        return False
    path = _record_path(name)
    if path is None or ensure_store_dir() is None:
        return False
    import tempfile
    try:
        fd, tmp = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp",
                                   dir=str(path.parent))
    except OSError:
        return False
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            json.dump(data, fh, ensure_ascii=False, sort_keys=True)
        os.chmod(tmp, 0o600)
        os.replace(tmp, str(path))
        return True
    except (OSError, TypeError, ValueError):
        try:
            os.unlink(tmp)
        except OSError:
            pass
        return False


__all__ = ["STORE_NAME", "home_store_dir", "store_dir", "ensure_store_dir",
           "account_key", "derived_key", "mac", "check_mac",
           "forget_cached_keys", "read_record", "write_record"]
