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
            if name in stripped and ("export" in stripped or "=" in stripped):
                found.append((path, number, line))
    return found


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
                    f"{_resolve_path(store_path)} (0600). An exported key can be "
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
    store = _read_store(path)
    rows: list[dict] = []
    for name in wanted:
        value = str(environment.get(name, "") or "")
        exports = [(str(f), n, line)
                   for f, n, line in exported_in_shell_files(name, home=home)]
        if not value:
            rows.append({"name": name, "action": "absent", "exports": exports})
            continue
        held = store.get(name, "")
        if held == value:
            action = "already"
        elif held:
            action = "differs"
        else:
            set_credential(name, value, path=path)
            action = "stored"
        rows.append({"name": name, "action": action, "exports": exports})
    return rows


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
]
