"""The pairing state: which contact is the user's phone.

Stored 0600 in a 0700 folder under $HOME/.delfin/messenger/pair.json:
the SimpleX contact id of the one allowed sender and, after ``pair``,
nothing else (no address, no name -- those live in the daemon's own
database, 0700 as well). ``revoke`` removes the file: the phone may no
longer command the agent.
"""

from __future__ import annotations

import json
import os
import stat
from pathlib import Path


def state_dir() -> Path:
    return Path(os.environ.get("DELFIN_MESSENGER_DIR", "")
                or (Path.home() / ".delfin" / "messenger"))


def pair_file() -> Path:
    return state_dir() / "pair.json"


def load() -> dict:
    """The pairing: {"contact_id": "..."}; {} when not paired."""
    path = pair_file()
    try:
        data = json.loads(path.read_text())
        return data if isinstance(data, dict) else {}
    except (OSError, ValueError):
        return {}


def save(contact_id: str) -> None:
    d = state_dir()
    d.mkdir(parents=True, exist_ok=True)
    os.chmod(d, 0o700)
    path = pair_file()
    path.write_text(json.dumps({"contact_id": str(contact_id)}) + "\n")
    os.chmod(path, 0o600)


def revoke() -> bool:
    """Unpair. True when something was there to remove."""
    path = pair_file()
    try:
        path.unlink()
        return True
    except FileNotFoundError:
        return False


def paired_contact() -> str | None:
    return str(load().get("contact_id") or "") or None


def _ensure_private(d: Path) -> None:
    try:
        mode = stat.S_IMODE(d.stat().st_mode)
        if mode & 0o077:
            os.chmod(d, 0o700)
    except OSError:
        pass
