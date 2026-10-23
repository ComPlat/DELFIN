"""Secure credential store for the DELFIN agent.

Two stores. The operating system's credential store (macOS Keychain,
Linux Secret Service, Windows Credential Manager), reached through the
optional ``keyring`` package, and ``~/.delfin/credentials.json`` with
file-mode ``0600`` (the convention of ``~/.netrc`` and
``~/.aws/credentials``). The engine calls ``load_credential(name)`` at
engine-init time only; no CLI command ever echoes a stored value back.

Lookup priority for any credential name:

  1. ``os.environ[name]``         — session-level override
  2. the keyring, service ``delfin`` — only for names recorded in
     ``~/.delfin/credentials.keyring.json`` (names, never values)
  3. ``~/.delfin/credentials.json`` — persistent file store
  4. ``""``                        — not configured

Where a write goes is ``agent.credentials.store`` in the user's own
settings file (``auto`` | ``keyring`` | ``file``; a workspace settings
file does not choose it). ``auto`` takes the keyring when a secure
backend answers within the timeout, else the file. Every keyring call
runs in a daemon thread bounded by a timeout (3 s without a terminal and
display, 60 s with one, so an unlock prompt can be answered); a call that
does not return marks the keyring unusable for the rest of the process,
so a headless run waits at most once. Backends below priority 1
(``fail``, ``null``, ``keyrings.alt`` encrypted/plaintext files) are not
a store. Rejected: probing ``secretstorage`` collections for a lock flag
-- backend-specific, and a timeout answers the same question on every
backend.

Nothing moves between the stores on its own: ``migrate_to_keyring``
runs on request and removes a file entry only after the keyring returned
the same value.

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
import re
import sys
import threading
from pathlib import Path


_DEFAULT_PATH = Path.home() / ".delfin" / "credentials.json"

# Names the agent stack knows how to consume.  Used by ``list_credentials``
# to surface env-var-sourced keys even when nothing is on disk.
_WELL_KNOWN_KEYS: tuple[str, ...] = (
    "KIT_TOOLBOX_API_KEY",
    "OPENAI_API_KEY",
    "ANTHROPIC_API_KEY",
    # A GitHub token opens the repository this agent pushes to.
    "GITHUB_TOKEN",
    "GH_TOKEN",
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


class CredentialStoreNotPrivate(OSError):
    """The key would not be private to the user where it would be stored."""


def _write_store(store: dict[str, str], path: Path | None = None) -> None:
    """Write the store so that nobody but the user can read it, or not at all.

    The temporary file was created under the umask and only then chmodded:
    a window in which it was world-readable on a default umask, and for
    good on a filesystem that ignores modes (CIFS, FAT, some network
    mounts). It is now created 0600 with O_EXCL and O_NOFOLLOW, and its
    mode and owner are checked on the open descriptor before a byte of the
    key is written. Where that check fails the key is not stored.
    """
    p = _resolve_path(path)
    try:
        p.parent.mkdir(parents=True, exist_ok=True)
        try:
            os.chmod(p.parent, 0o700)
        except OSError:
            pass
    except OSError:
        return
    tmp = p.with_suffix(p.suffix + ".tmp")
    try:
        tmp.unlink()
    except FileNotFoundError:
        pass
    except OSError:
        return
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0)
    try:
        fd = os.open(tmp, flags, 0o600)
    except OSError:
        return
    try:
        st = os.fstat(fd)
        if (st.st_mode & 0o077) or (hasattr(os, "getuid") and st.st_uid != os.getuid()):
            raise CredentialStoreNotPrivate(
                f"{p.parent} does not keep files private to you (mode "
                f"{oct(st.st_mode & 0o777)}); the key was not stored. Keep "
                "it in an environment variable instead, or move the store "
                "to a filesystem with POSIX permissions.")
        data = (json.dumps(store, ensure_ascii=False, indent=2, sort_keys=True)
                + "\n").encode("utf-8")
        os.write(fd, data)
        os.fsync(fd)
    except CredentialStoreNotPrivate:
        os.close(fd)
        try:
            tmp.unlink()
        except OSError:
            pass
        raise
    except OSError:
        os.close(fd)
        try:
            tmp.unlink()
        except OSError:
            pass
        return
    os.close(fd)
    try:
        tmp.replace(p)
    except OSError:
        try:
            tmp.unlink()
        except OSError:
            pass


# ---------------------------------------------------------------------------
# The operating system's credential store, through the optional keyring
# ---------------------------------------------------------------------------

#: Service name every DELFIN entry is filed under; the username is the
#: credential name.
KEYRING_SERVICE = "delfin"
STORE_MODES: tuple[str, ...] = ("auto", "keyring", "file")
#: Seconds a keyring call may take without a terminal and display, where
#: an unlock prompt has nobody to answer it.
HEADLESS_TIMEOUT_S = 3.0
#: Seconds with a terminal and display: long enough to type a password
#: into the unlock dialog.
INTERACTIVE_TIMEOUT_S = 60.0
_PROBE_USER = "__delfin_probe__"


class CredentialStoreUnavailable(RuntimeError):
    """The keyring cannot take or give this key; the message says why and
    never carries a value."""


class _KeyringTimeout(Exception):
    pass


def configured_store() -> str:
    """``agent.credentials.store`` from the user's settings file.

    Read through ``user_settings.load_settings`` -- the file in the home
    directory, outside every writable root and on the protected list -- and
    never from a workspace's ``.delfin/settings.json``. Anything
    unrecognised reads as ``auto``.
    """
    try:
        from delfin.user_settings import load_settings
        agent = (load_settings() or {}).get("agent") or {}
        section = agent.get("credentials") if isinstance(agent, dict) else None
        raw = section.get("store") if isinstance(section, dict) else ""
    except Exception:
        return "auto"
    value = str(raw or "").strip().lower()
    return value if value in STORE_MODES else "auto"


def _interactive() -> bool:
    """A person is at a terminal who can answer an unlock dialog.

    Both standard streams on a tty, and on Linux/BSD a display as well:
    the Secret Service unlock prompt is a graphical window, so an ssh
    session to a cluster node is headless whatever its tty says.
    """
    try:
        tty = bool(sys.stdin and sys.stdin.isatty()
                   and sys.stdout and sys.stdout.isatty())
    except Exception:
        tty = False
    if not tty:
        return False
    if sys.platform == "darwin" or sys.platform.startswith("win"):
        return True
    return bool(os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY"))


def _keyring_timeout() -> float:
    return INTERACTIVE_TIMEOUT_S if _interactive() else HEADLESS_TIMEOUT_S


def _default_keyring_factory():
    """The backend ``keyring`` selects (honours PYTHON_KEYRING_BACKEND)."""
    import keyring
    return keyring.get_keyring()


def _disabled_keyring_factory():
    raise ImportError("keyring disabled in this process")


#: Replaced by the test suite with ``_disabled_keyring_factory``.
_keyring_factory = _default_keyring_factory

_KR_LOCK = threading.RLock()
_KR: dict = {"probed": False, "backend": None, "name": "", "reason": ""}
_KR_CACHE: dict[str, str] = {}
_KR_TRIED: set[str] = set()


def reset_keyring_state() -> None:
    """Forget the probe verdict and every value read in this process."""
    with _KR_LOCK:
        _KR.update(probed=False, backend=None, name="", reason="")
        _KR_CACHE.clear()
        _KR_TRIED.clear()


def _bounded(fn, *args, timeout: float):
    """``fn(*args)`` in a daemon thread; raises _KeyringTimeout when it has
    not returned after *timeout* seconds. The thread is left to finish on
    its own -- a blocked D-Bus call cannot be cancelled, and a daemon
    thread does not hold the process open at exit."""
    box: dict = {}

    def run():
        try:
            box["value"] = fn(*args)
        except BaseException as exc:  # noqa: BLE001 - re-raised below
            box["error"] = exc

    t = threading.Thread(target=run, name="delfin-keyring", daemon=True)
    t.start()
    t.join(timeout)
    if t.is_alive():
        raise _KeyringTimeout()
    if "error" in box:
        raise box["error"]
    return box.get("value")


def _backend_label(backend) -> str:
    cls = type(backend)
    return f"{cls.__module__}.{cls.__name__}"


def _insecure_reason(backend) -> str:
    """Why *backend* is not a store for a key, or "".

    keyring's own ranking: priority >= 1 is a recommended backend; the
    ``fail`` (0) and ``null`` (-1) backends and the ``keyrings.alt`` file
    backends (< 1, and the encrypted one asks for a password on stdin)
    are below it. A class named Plaintext is refused at any priority.
    """
    label = _backend_label(backend)
    if "plaintext" in label.lower():
        return f"{label} keeps keys in clear text"
    try:
        priority = float(getattr(backend, "priority", 0))
    except Exception:
        priority = 0.0
    if priority < 1:
        return f"{label} is not a secure store (priority {priority:g})"
    return ""


def _mark_unusable(reason: str) -> None:
    with _KR_LOCK:
        _KR.update(probed=True, backend=None, name="", reason=reason)


def keyring_status() -> dict:
    """``{usable, backend, reason}`` for this process; probed once.

    Usable means: the package imports, the selected backend is a secure
    one, and a read of a probe entry returned within the timeout. The
    verdict is kept, so a host where the keyring hangs pays the timeout
    once per process, not once per key.
    """
    with _KR_LOCK:
        if not _KR["probed"]:
            timeout = _keyring_timeout()
            reason = ""
            backend = None
            try:
                backend = _bounded(_keyring_factory, timeout=timeout)
            except _KeyringTimeout:
                reason = (f"the keyring did not load within {timeout:g} s "
                          "(timed out)")
            except ImportError:
                reason = ("the keyring package is not installed "
                          "(pip install 'delfin-complat[keyring]')")
            except Exception as exc:
                reason = f"the keyring could not be loaded ({type(exc).__name__})"
            if not reason and backend is None:
                reason = "no keyring backend"
            if not reason:
                reason = _insecure_reason(backend)
            if not reason:
                try:
                    _bounded(backend.get_password, KEYRING_SERVICE,
                             _PROBE_USER, timeout=timeout)
                except _KeyringTimeout:
                    reason = (f"{_backend_label(backend)} timed out after "
                              f"{timeout:g} s (locked, or waiting for an "
                              "unlock prompt nobody can answer here)")
                except Exception as exc:
                    reason = (f"{_backend_label(backend)} refused a read "
                              f"({type(exc).__name__})")
            if reason:
                _KR.update(probed=True, backend=None, name="", reason=reason)
            else:
                _KR.update(probed=True, backend=backend,
                           name=_backend_label(backend), reason="")
        return {"usable": _KR["backend"] is not None,
                "backend": _KR["name"], "reason": _KR["reason"]}


def _usable_backend():
    status = keyring_status()
    if not status["usable"]:
        raise CredentialStoreUnavailable(status["reason"])
    return _KR["backend"]


def _kr_call(fn, *args):
    timeout = _keyring_timeout()
    try:
        return _bounded(fn, *args, timeout=timeout)
    except _KeyringTimeout:
        reason = (f"{_KR['name'] or 'the keyring'} timed out after "
                  f"{timeout:g} s (timed out)")
        _mark_unusable(reason)
        raise CredentialStoreUnavailable(reason) from None
    except CredentialStoreUnavailable:
        raise
    except Exception as exc:
        raise CredentialStoreUnavailable(
            f"{_KR['name'] or 'the keyring'} refused the call "
            f"({type(exc).__name__})") from None


def _kr_get(name: str, *, cached: bool = True) -> str:
    if cached and name in _KR_CACHE:
        return _KR_CACHE[name]
    backend = _usable_backend()
    value = _kr_call(backend.get_password, KEYRING_SERVICE, name)
    value = value if isinstance(value, str) else ""
    if value:
        _KR_CACHE[name] = value
    else:
        _KR_CACHE.pop(name, None)
    return value


def _kr_set(name: str, value: str) -> None:
    """Store and read back; raises when the keyring returns anything else.
    A copy that did not read back is removed again (best effort)."""
    backend = _usable_backend()
    _kr_call(backend.set_password, KEYRING_SERVICE, name, value)
    _KR_CACHE.pop(name, None)
    back = _kr_get(name, cached=False)
    if back != value:
        try:
            _kr_call(backend.delete_password, KEYRING_SERVICE, name)
        except CredentialStoreUnavailable:
            pass
        _KR_CACHE.pop(name, None)
        raise CredentialStoreUnavailable(
            f"{_KR['name']} returned a different value on read-back; "
            "the key was not stored there")


def _kr_delete(name: str) -> None:
    backend = _usable_backend()
    try:
        _bounded(backend.delete_password, KEYRING_SERVICE, name,
                 timeout=_keyring_timeout())
    except _KeyringTimeout:
        reason = f"{_KR['name']} timed out (timed out)"
        _mark_unusable(reason)
        raise CredentialStoreUnavailable(reason) from None
    except Exception as exc:
        # Already absent is the outcome asked for.
        if type(exc).__name__ != "PasswordDeleteError":
            raise CredentialStoreUnavailable(
                f"{_KR['name']} refused the delete ({type(exc).__name__})"
            ) from None
    _KR_CACHE.pop(name, None)


def index_path(path: Path | None = None) -> Path:
    """The file naming which credentials are in the keyring. Names only:
    keyring backends cannot list their entries, and without this record
    every lookup of an unset name would have to ask the keyring."""
    p = _resolve_path(path)
    return p.with_name(p.stem + ".keyring.json")


def _index_names(path: Path | None = None) -> list[str]:
    p = index_path(path)
    try:
        raw = json.loads(p.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return []
    names = raw.get("names") if isinstance(raw, dict) else None
    if not isinstance(names, list):
        return []
    return [n for n in names if isinstance(n, str) and n]


def _index_write(names: list[str], path: Path | None = None) -> None:
    _write_store({"service": KEYRING_SERVICE,
                  "names": sorted(set(names))}, index_path(path))


def _index_add(name: str, path: Path | None = None) -> None:
    names = _index_names(path)
    if name not in names:
        _index_write(names + [name], path)


def _index_remove(name: str, path: Path | None = None) -> None:
    names = _index_names(path)
    if name in names:
        _index_write([n for n in names if n != name], path)


def _write_target(store: str | None) -> str:
    """"keyring" or "file" for a write under *store* (default: the
    setting). An explicit keyring setting with no usable keyring raises:
    the user asked for the key not to be on disk."""
    mode = store if store in STORE_MODES else configured_store()
    if mode == "file":
        return "file"
    status = keyring_status()
    if status["usable"]:
        return "keyring"
    if mode == "keyring":
        raise CredentialStoreUnavailable(
            f"agent.credentials.store is keyring and {status['reason']}; "
            "the key was not stored")
    return "file"


def stored_in(name: str, *, path: Path | None = None) -> list[str]:
    """Which stores hold *name*: "keyring" (per the index) and/or "file"."""
    where = []
    if name in _index_names(path):
        where.append("keyring")
    if name in _read_store(path):
        where.append("file")
    return where


def describe_store(name: str = "", *, path: Path | None = None) -> str:
    """Where *name* lives, for a message: the keyring backend or the file."""
    if name and name in _index_names(path):
        return f"the system keyring ({_KR['name'] or 'keyring'})"
    return f"{_resolve_path(path)} (0600)"


def load_credential(name: str, *, path: Path | None = None) -> str:
    """Return the credential value: env-var > keyring > file > ''.

    The keyring is asked only for names the index records, so a host
    without one never imports it. When it cannot answer, the file copy (if
    any) is returned. The lookup is silent — callers MUST handle an empty
    return without revealing the credential name in any user-visible error.
    """
    if not name:
        return ""
    env_val = os.environ.get(name, "")
    if env_val:
        return env_val
    if name in _index_names(path):
        try:
            value = _kr_get(name)
        except CredentialStoreUnavailable:
            value = ""
        if value:
            return value
    store = _read_store(path)
    return store.get(name, "")


def set_credential(
    name: str, value: str, *, path: Path | None = None,
    store: str | None = None,
) -> bool:
    """Persist a credential.  Returns True if the store was modified.

    *store* overrides the setting for this call ("auto", "keyring",
    "file"). Under auto a keyring that fails the write leaves the key in
    the file; under keyring it raises CredentialStoreUnavailable and
    nothing is written. A file copy of the same name is left in place --
    the keyring copy is read first; ``migrate_to_keyring`` or ``delete``
    remove it on request. Empty ``value`` is a no-op (use
    ``delete_credential`` to remove).
    """
    if not name or not value:
        return False
    mode = store if store in STORE_MODES else configured_store()
    if _write_target(mode) == "keyring":
        try:
            if name in _index_names(path) and _kr_get(name) == value:
                return False
            _kr_set(name, value)
            _index_add(name, path)
            return True
        except CredentialStoreUnavailable:
            if mode == "keyring":
                raise
    current = _read_store(path)
    if current.get(name) == value:
        return False
    current[name] = value
    _write_store(current, path)
    return True


def delete_credential(name: str, *, path: Path | None = None) -> bool:
    """Remove a credential from every store.  Returns True if it was present.

    Raises CredentialStoreUnavailable when the keyring holds it and cannot
    be reached; the file copy is removed first either way.
    """
    if not name:
        return False
    removed = False
    current = _read_store(path)
    if name in current:
        del current[name]
        _write_store(current, path)
        removed = True
    if name in _index_names(path):
        _kr_delete(name)
        _index_remove(name, path)
        removed = True
    return removed


def migrate_to_keyring(
    names: "tuple[str, ...] | list[str]" = (), *, path: Path | None = None,
) -> list[dict]:
    """Move file entries into the keyring. One row per name.

    Row: ``{name, action, reason}``; action is "migrated" (the keyring
    returned the same value and the file entry was removed), "differs"
    (the keyring holds another value for the name -- nothing changed),
    "failed" (the write or read-back failed -- the file entry stays) or
    "absent" (not in the file). Raises CredentialStoreUnavailable, with
    nothing changed, when the setting is ``file`` or no secure keyring is
    usable. The file entry is removed only if it still holds the value
    that was verified.
    """
    if configured_store() == "file":
        raise CredentialStoreUnavailable(
            "agent.credentials.store is file; set it to auto or keyring "
            "to move keys into the keyring")
    _usable_backend()
    current = _read_store(path)
    wanted = list(names) or sorted(current)
    rows: list[dict] = []
    for name in wanted:
        value = current.get(name, "")
        if not value:
            rows.append({"name": name, "action": "absent", "reason": ""})
            continue
        try:
            held = _kr_get(name, cached=False)
            if held and held != value:
                rows.append({"name": name, "action": "differs",
                             "reason": "the keyring holds a different value"})
                continue
            if held != value:
                _kr_set(name, value)
            if _kr_get(name, cached=False) != value:
                raise CredentialStoreUnavailable("read-back differs")
            _index_add(name, path)
        except CredentialStoreUnavailable as exc:
            rows.append({"name": name, "action": "failed", "reason": str(exc)})
            continue
        latest = _read_store(path)
        if latest.get(name) == value:
            del latest[name]
            _write_store(latest, path)
        rows.append({"name": name, "action": "migrated", "reason": ""})
    return rows


def held_secret_values(*, path: Path | None = None) -> list[str]:
    """Every stored credential value, both stores, for exact-match
    redaction and egress checks. Never raises. Each keyring name is asked
    at most once per process; afterwards the cached value is used."""
    out = [v for v in _read_store(path).values() if v]
    try:
        names = _index_names(path)
    except Exception:
        names = []
    for name in names:
        if name in _KR_CACHE:
            out.append(_KR_CACHE[name])
            continue
        if name in _KR_TRIED:
            continue
        _KR_TRIED.add(name)
        try:
            value = _kr_get(name)
        except Exception:
            continue
        if value:
            out.append(value)
    return out


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
    """Return name → {value: masked, source, stores}.

    ``source`` is where ``load_credential`` takes the value from: "env",
    "keyring" or "file" (env > keyring > file); ``stores`` lists every
    place holding it. A keyring entry that cannot be read right now is
    shown as "(unavailable)", never as an error.
    """
    out: dict = {}
    store = _read_store(path)
    for name, val in store.items():
        out[name] = {"value": mask(val), "source": "file", "stores": ["file"]}
    for name in _index_names(path):
        try:
            val = _kr_get(name)
            shown = mask(val) if val else "(missing)"
        except CredentialStoreUnavailable:
            shown = "(unavailable)"
        stores = ["keyring"] + (["file"] if name in store else [])
        out[name] = {"value": shown, "source": "keyring", "stores": stores}
    for name in _WELL_KNOWN_KEYS:
        env_val = os.environ.get(name, "")
        if not env_val:
            continue
        # env-var wins on conflict
        prior = (out.get(name) or {}).get("stores") or []
        out[name] = {"value": mask(env_val), "source": "env",
                     "stores": ["env"] + list(prior)}
    return out


#: Where a login shell usually exports things, most specific first.
SHELL_FILES: tuple[str, ...] = (
    ".bashrc", ".bash_profile", ".bash_login", ".profile", ".zshrc",
    ".zshenv", ".zprofile", ".kshrc",
    # Not shells: the places a LOGIN SESSION is given variables. A key
    # that reaches every process of the session without appearing in any
    # rc file comes from one of these (seen on a cluster, 2026-09-17,
    # where no shell file named it and every process still had it).
    ".ssh/environment", ".pam_environment",
    ".config/environment.d/*.conf",
)


def _sets(name: str, line: str) -> bool:
    """Whether *line* sets or exports *name* -- not merely reads it.

    ``OPENAI_API_KEY=$KIT_TOOLBOX_API_KEY`` names the second key and sets
    only the first. Taking every line that held the name and an ``=``
    anywhere for an export commented out an alias that read the key from
    a 0600 file (2026-09-20).
    """
    n = re.escape(name)
    return bool(
        re.search(rf"(?<![\w${{]){n}=", line)
        or re.search(rf"\b(?:export|declare\s+-\w*x\w*|typeset\s+-\w*x\w*)"
                     rf"\s+(?:[^\s;&|]+\s+)*{n}\b", line)
        or re.search(rf"\bsetenv\s+{n}\b", line)
        or re.search(rf"^\s*{n}\s+(?:DEFAULT|OVERRIDE)=", line))


def _written_value(name: str, line: str) -> str:
    """The value *line* writes for *name*, as it stands in the file.

    Empty when there is none to find there: a bare ``export NAME``, an
    empty value, or one that is a ``$REFERENCE`` or a command -- those put
    no key into the file, whoever can read it.
    """
    n = re.escape(name)
    m = (re.search(rf"(?<![\w${{]){n}=(\S*)", line)
         or re.search(rf"\bsetenv\s+{n}\s+(\S+)", line)
         or re.search(rf"^\s*{n}\s+(?:DEFAULT|OVERRIDE)=(\S*)", line))
    if not m:
        return ""
    value = m.group(1).lstrip("'\"")
    if not value or value.startswith(("$", "`")):
        return ""
    return value


def exported_in_shell_files(name: str, *, home: Path | None = None
                            ) -> list[tuple[Path, int, str]]:
    """Every (file, line number, line) that exports ``name``.

    Read-only: this finds the line, it does not touch it. A shell file is
    the user's own, and a key that a tool silently edits out of it is a
    tool that edits dotfiles.
    """
    root = Path(home) if home is not None else Path.home()
    found: list[tuple[Path, int, str]] = []
    if not name:
        return found
    candidates: list[Path] = []
    for rel in SHELL_FILES:
        if "*" in rel:
            try:
                candidates.extend(sorted(root.glob(rel)))
            except OSError:
                continue
        else:
            candidates.append(root / rel)
    for path in candidates:
        try:
            lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        except OSError:
            continue
        for number, line in enumerate(lines, start=1):
            stripped = line.strip()
            if stripped.startswith("#"):
                continue
            if name in stripped and _sets(name, stripped):
                found.append((path, number, line))
    return found


def _path_traversable(directory: Path, bit: int) -> bool:
    """Can somebody holding *bit* walk from the root down to *directory*?

    Every directory on the way needs its execute bit for that class; one
    without it shuts the door for everything below, however open the rest
    of the path is. Unreadable or vanished is answered NO -- a check that
    cannot see is not a check that found nothing, and for this question
    the safe side is to stay quiet rather than to raise an alarm about a
    path it could not walk.
    """
    try:
        current = Path(directory).resolve()
    except OSError:
        return False
    while True:
        try:
            if not (current.stat().st_mode & bit):
                return False
        except OSError:
            return False
        if current.parent == current:
            return True
        current = current.parent


def keys_in_readable_shell_files(
    names: "tuple[str, ...] | list[str]" = (),
    *,
    home: Path | None = None,
) -> list[dict]:
    """Provider keys sitting in shell files that somebody else can read.

    The environment is the wrong place to start looking. ``adopt_from_
    environment`` searches the shell files only for keys it already found
    in ``os.environ``, so a key that never reaches the environment is
    never searched for -- and that is the worse case, because it stays in
    the file instead of living only as long as a shell. The shape that
    hides this way is an alias:

        alias codex-kit='OPENAI_API_KEY="..." ... codex'

    The value is set for the duration of one command and is absent from
    every environment DELFIN inspects, while the line itself is readable
    by everyone who may read the file.

    The permission is the finding, not the presence: a key in a 0600
    ``~/.bashrc`` is where its owner put it and nobody else sees it. Each
    row is ``{name, file, line_no, mode, readable_by}`` and never carries
    the value. Read-only, always -- the value is not in the store, and a
    tool that removes the only copy of a key is worse than the leak.
    """
    rows: list[dict] = []
    for name in (tuple(names) or _WELL_KNOWN_KEYS):
        for path, number, line in exported_in_shell_files(name, home=home):
            if not _written_value(name, line):
                continue
            try:
                mode = path.stat().st_mode & 0o777
            except OSError:
                continue
            others = bool(mode & 0o004)
            group = bool(mode & 0o040)
            if not (group or others):
                continue
            # A file is reached by walking the directories above it, and
            # each of those needs its execute bit for the asker. A 0644
            # file inside a 0700 home is open to nobody, whatever its own
            # bits say -- and a 0700 home is the ordinary configuration,
            # so reading the file's mode alone raised this on almost every
            # host while naming a door that is shut. A security check
            # people learn to scroll past has stopped being one.
            group = group and _path_traversable(path.parent, 0o010)
            others = others and _path_traversable(path.parent, 0o001)
            if not (group or others):
                continue
            who = ("group, others" if group and others
                   else "others" if others else "group")
            rows.append({
                "name": name,
                "file": str(path),
                "line_no": number,
                "mode": f"{mode:04o}",
                "readable_by": who,
            })
    return rows


#: Suffix of the copy kept beside a shell file this touches.
SHELL_BACKUP_SUFFIX = ".delfin-backup"


def comment_out_exports(name: str, *, home: Path | None = None,
                        store_path: Path | None = None) -> list[dict]:
    """Comment out every line exporting ``name``. Returns what was done.

    A copy of the file is written beside it first, so the change is one
    ``mv`` away from undone, and the line stays in place as a comment
    rather than disappearing: a user reading their own shell file must be
    able to see what happened and why.
    """
    done: list[dict] = []
    for path, number, line in exported_in_shell_files(name, home=home):
        # An alias sets a variable for one command and exports nothing
        # into the shell, so commenting it out stops no export -- it only
        # breaks the alias. A key written in one is reported by
        # keys_in_readable_shell_files instead, and never edited.
        if line.lstrip().startswith("alias "):
            continue
        try:
            text = path.read_text(encoding="utf-8")
            lines = text.splitlines(keepends=True)
            index = number - 1
            if index >= len(lines) or name not in lines[index]:
                continue
            backup = path.with_suffix(path.suffix + SHELL_BACKUP_SUFFIX)
            if not backup.exists():
                # 0600, NOT the original file's mode: the copy still holds
                # the export, and inheriting a group-readable shell file's
                # permissions would move the key rather than protect it.
                fd = os.open(str(backup),
                             os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
                try:
                    with os.fdopen(fd, "w", encoding="utf-8") as handle:
                        handle.write(text)
                except Exception:
                    try:
                        os.close(fd)
                    except OSError:
                        pass
                    raise
                try:
                    os.chmod(backup, 0o600)
                except OSError:
                    pass
            note = (f"# disabled by DELFIN: {name} now lives in "
                    f"{describe_store(name, path=store_path)}. An exported key can be "
                    "read by every process of yours through /proc.\n")
            lines[index] = note + "# " + lines[index]
            path.write_text("".join(lines), encoding="utf-8")
            done.append({"file": str(path), "line": number,
                         "backup": str(backup)})
        except OSError:
            continue
    return done


def systemd_user_environment(name: str) -> bool:
    """Whether the systemd user session hands ``name`` to its children."""
    try:
        from . import contained_run as _contained
        proc = _contained.run(["systemctl", "--user", "show-environment"],
                              timeout=10)
    except Exception:
        return False
    if proc.returncode != 0:
        return False
    for line in str(proc.stdout or "").splitlines():
        if line.split("=", 1)[0].strip() == name:
            return True
    return False


def unset_in_systemd_user_environment(name: str) -> bool:
    """Take ``name`` out of the systemd user environment. Reversible.

    Only the user's own session is touched, and only for processes it
    starts from now on -- the one already running keeps what it was given
    at exec, which is why this is paired with a new-shell notice rather
    than offered as a cure for the current process.
    """
    if not systemd_user_environment(name):
        return False
    try:
        from . import contained_run as _contained
        proc = _contained.run(
            ["systemctl", "--user", "unset-environment", name], timeout=10)
        return proc.returncode == 0
    except Exception:
        return False


def secure_exported_keys(
    names: "tuple[str, ...] | list[str]" = (),
    *,
    path: Path | None = None,
    env: "dict | None" = None,
    home: Path | None = None,
    clean_shell: bool = True,
) -> list[dict]:
    """Put exported provider keys where they belong, and say what was done.

    Telling a user that their key is exported, and then leaving them four
    manual steps, is a warning that repeats every start. This takes the
    key into the 0600 store and comments out the line that exports it,
    keeping a copy of the file. What it will not do is touch a key whose
    stored value differs from the exported one -- that is a question, not
    a chore.
    """
    rows = adopt_from_environment(names, path=path, env=env, home=home)
    for row in rows:
        row["cleaned"] = []
        if not clean_shell:
            continue
        row["systemd"] = False
        if row["action"] in ("stored", "already"):
            if row["exports"]:
                row["cleaned"] = comment_out_exports(
                    row["name"], home=home, store_path=path)
            # And the session's own environment, where nothing in the
            # home directory mentions the key and every process has it.
            if env is None:
                row["systemd"] = unset_in_systemd_user_environment(row["name"])
    return rows


def adopt_from_environment(
    names: "tuple[str, ...] | list[str]" = (),
    *,
    path: Path | None = None,
    env: "dict | None" = None,
    home: Path | None = None,
) -> list[dict]:
    """Take exported keys into the store. Returns one row per name.

    Each row: ``{name, action, exports}`` where action is "stored" (it
    was not in the store), "already" (the store holds the same value),
    "differs" (the store holds a DIFFERENT value — nothing is
    overwritten) or "absent" (not exported at all), and ``exports`` lists
    the shell lines still exporting it.

    Nothing is deleted and no shell file is written: the point is to make
    the store the place the key lives, and then to say exactly which line
    the user should remove.
    """
    import os as _os

    environment = _os.environ if env is None else env
    wanted = tuple(names) or _WELL_KNOWN_KEYS
    rows: list[dict] = []
    for name in wanted:
        value = str(environment.get(name, "") or "")
        exports = [(str(f), n, line)
                   for f, n, line in exported_in_shell_files(name, home=home)]
        if not value:
            rows.append({"name": name, "action": "absent", "exports": exports})
            continue
        held = _stored_value(name, path)
        if held == value:
            action = "already"
        elif held:
            action = "differs"
        else:
            set_credential(name, value, path=path)
            action = "stored"
        rows.append({"name": name, "action": action, "exports": exports})
    return rows


def _stored_value(name: str, path: Path | None = None) -> str:
    """The value a store holds for *name*, ignoring the environment."""
    if name in _index_names(path):
        try:
            value = _kr_get(name)
        except CredentialStoreUnavailable:
            value = ""
        if value:
            return value
    return _read_store(path).get(name, "")


def credentials_path() -> Path:
    """Public accessor for the credentials file location."""
    return _DEFAULT_PATH


__all__ = [
    "adopt_from_environment",
    "unset_in_systemd_user_environment",
    "secure_exported_keys",
    "comment_out_exports",
    "exported_in_shell_files",
    "SHELL_FILES",
    "load_credential",
    "set_credential",
    "delete_credential",
    "list_credentials",
    "mask",
    "credentials_path",
    "index_path",
    "configured_store",
    "keyring_status",
    "reset_keyring_state",
    "stored_in",
    "describe_store",
    "migrate_to_keyring",
    "held_secret_values",
    "CredentialStoreUnavailable",
    "CredentialStoreNotPrivate",
    "KEYRING_SERVICE",
    "STORE_MODES",
]
