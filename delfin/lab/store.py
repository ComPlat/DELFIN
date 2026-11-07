"""Where a lab installation keeps its state, and who is acting.

Input: the environment (``DELFIN_LAB_DIR``, ``HOME``, ``DELFIN_PERSON``).
Output: the installation's lab directory, the acting person, and a
cross-process transaction over the shared state file.

The lab directory belongs to the INSTALLATION, not to a person: in the
launcher's people mode every person has ``<account>/people/<name>`` as HOME
on one account, and all of them must see one set of locks, one e-stop and
one manifest. So a HOME of that shape resolves to the account's own
``~/.delfin/lab``; any other HOME (a workstation, a test) to
``$HOME/.delfin/lab``. ``DELFIN_LAB_DIR`` overrides both (the systemd units
set it).

The person is ``DELFIN_PERSON`` when set, else ``<name>`` of a people-mode
HOME, else the account name. Same account means same Unix user: the person
is an attribution and a rule the kernel applies (a lock or an approval of
one person is never another's), not an operating-system boundary between
them. The boundary between the agent and the person is the agent-session
marker (``approval_guard.ENV_MARKER``) and the answer-room protection.

State (``state.json``) is one JSON document changed under ``flock`` on
``state.lock`` and replaced atomically (write to a temp file in the same
directory, fsync, rename), so a kill mid-write leaves the previous state.
"""

from __future__ import annotations

import contextlib
import json
import os
import re
import socket
import tempfile
import time
from pathlib import Path
from typing import Iterator, Optional

try:
    import fcntl as _fcntl
except ImportError:                 # pragma: no cover - non-POSIX
    _fcntl = None

_PERSON_RE = re.compile(r"^[a-z0-9_-]{1,32}$")


def lab_dir() -> Path:
    env = os.environ.get("DELFIN_LAB_DIR", "").strip()
    if env:
        return Path(env).expanduser()
    home = Path.home()
    if home.parent.name == "people":
        return home.parent.parent / ".delfin" / "lab"
    return home / ".delfin" / "lab"


def ensure_lab_dir() -> Path:
    d = lab_dir()
    d.mkdir(parents=True, exist_ok=True, mode=0o700)
    try:
        os.chmod(d, 0o700)
    except OSError:
        pass
    return d


def manifest_path() -> Path:
    env = os.environ.get("DELFIN_LAB_MANIFEST", "").strip()
    return Path(env).expanduser() if env else lab_dir() / "manifest.yaml"


def lab_configured() -> bool:
    """True when this installation has a lab manifest (tools advertised)."""
    try:
        return manifest_path().is_file()
    except OSError:
        return False


def person() -> str:
    env = os.environ.get("DELFIN_PERSON", "").strip().lower()
    if env and _PERSON_RE.match(env):
        return env
    home = Path.home()
    if home.parent.name == "people" and _PERSON_RE.match(home.name):
        return home.name
    try:
        import pwd
        return pwd.getpwuid(os.getuid()).pw_name
    except Exception:
        return os.environ.get("USER", "unknown") or "unknown"


def host() -> str:
    try:
        return socket.gethostname()
    except Exception:
        return ""


def pid_alive(pid: int) -> bool:
    try:
        os.kill(int(pid), 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    except (OSError, ValueError, TypeError):
        return False
    try:
        # A zombie (killed, not yet reaped) is dead for our purposes.
        with open(f"/proc/{int(pid)}/stat", "rb") as fh:
            return fh.read().split(b") ", 1)[-1][:1] != b"Z"
    except OSError:
        return True


def _write_atomic(path: Path, data: dict) -> None:
    fd, tmp = tempfile.mkstemp(prefix=".state-", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            json.dump(data, fh, sort_keys=True)
            fh.flush()
            os.fsync(fh.fileno())
        os.chmod(tmp, 0o600)
        os.replace(tmp, path)
    except BaseException:
        with contextlib.suppress(OSError):
            os.unlink(tmp)
        raise


def read_json(path: Path, default: dict) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return json.loads(json.dumps(default))


@contextlib.contextmanager
def locked_json(name: str, default: Optional[dict] = None) -> Iterator[dict]:
    """``with locked_json("state.json") as doc:`` -- doc is written back."""
    d = ensure_lab_dir()
    path = d / name
    lock = d / (name + ".lock")
    fd = os.open(str(lock), os.O_RDWR | os.O_CREAT, 0o600)
    try:
        if _fcntl is not None:
            _fcntl.flock(fd, _fcntl.LOCK_EX)
        doc = read_json(path, default or {})
        yield doc
        _write_atomic(path, doc)
    finally:
        if _fcntl is not None:
            with contextlib.suppress(OSError):
                _fcntl.flock(fd, _fcntl.LOCK_UN)
        os.close(fd)


_STATE_DEFAULT = {"locks": {}, "module_locks": {}, "active": {},
                  "estop": {"global": None, "modules": {}}, "latched": {},
                  "approvals": {}, "last_set": {}, "nonces": {}}


@contextlib.contextmanager
def state() -> Iterator[dict]:
    with locked_json("state.json", _STATE_DEFAULT) as doc:
        for k, v in _STATE_DEFAULT.items():
            doc.setdefault(k, json.loads(json.dumps(v)))
        yield doc


def snapshot() -> dict:
    doc = read_json(lab_dir() / "state.json", _STATE_DEFAULT)
    for k, v in _STATE_DEFAULT.items():
        doc.setdefault(k, json.loads(json.dumps(v)))
    return doc


def now() -> float:
    return time.time()
