"""Who a phone can belong to: the account owner or one person of people mode.

People mode (the Windows launcher, ``~/.config/delfin/launcher-people.conf``
in the account's home) gives each person ``~/people/<name>`` as HOME: own
sessions, agent service (``<home>/.run``), messenger activation records and
pairing requests. Everyone is the same Unix user, so the separation drawn
here is an application-level boundary: the gateway and the bridge never
route one person's phone into another person's sessions, but a process of
the account can read every person's files.

Input: a person name as the bridge reports it (the directory a pairing
request was found in). Output: that person's home, or None when the name is
not a person of this account. The owner is the empty name ``""``: the
account's own home and the sessions outside ``~/people``.
"""
from __future__ import annotations

import os
import re
import stat
from pathlib import Path

#: The account owner: the person an admin pairs with ``pi_messenger.sh --pair``.
OWNER = ""

#: The launcher's person names (tools/windows-launcher/remote_bootstrap.sh:
#: lower-case letters, digits, ``_`` and ``-``, not starting with ``-``).
_NAME = re.compile(r"[a-z0-9_][a-z0-9_-]{0,31}")

PEOPLE_CONF = Path(".config") / "delfin" / "launcher-people.conf"


def valid(name: object) -> bool:
    return isinstance(name, str) and (name == OWNER or bool(_NAME.fullmatch(name)))


def account_home() -> Path:
    """The account's home: the HOME of the bridge (the service user)."""
    return Path.home()


def mode_on(account: Path | None = None) -> bool:
    return ((account or account_home()) / PEOPLE_CONF).is_file()


def home(name: str, account: Path | None = None) -> Path | None:
    """The person's home, or None. A person's home must be a real directory
    (never a link) of this uid under ``<account>/people``, and people mode
    must be on."""
    account = account or account_home()
    if name == OWNER:
        return account
    if not valid(name) or not mode_on(account):
        return None
    path = account / "people" / name
    try:
        info = os.lstat(path)
    except OSError:
        return None
    if not stat.S_ISDIR(info.st_mode) or info.st_uid != os.getuid():
        return None
    return path


def names(account: Path | None = None) -> list[str]:
    """The owner, then every person of people mode, sorted."""
    account = account or account_home()
    found = [OWNER]
    if mode_on(account):
        try:
            entries = sorted(os.listdir(account / "people"))
        except OSError:
            entries = []
        found += [n for n in entries if n != OWNER and home(n, account) is not None]
    return found


def runtime_socket(name: str, account: Path | None = None) -> Path | None:
    """A person's agent-service socket (``<home>/.run/delfin-<uid>/agent.sock``,
    the launcher's DELFIN_RUNTIME_DIR); None for the owner (the bridge's
    default) and for an unknown name."""
    if name == OWNER:
        return None
    base = home(name, account)
    return None if base is None else base / ".run" / f"delfin-{os.getuid()}" / "agent.sock"


def runtime_dir(name: str, account: Path | None = None) -> Path | None:
    sock = runtime_socket(name, account)
    return None if sock is None else sock.parent


def activation_dir(name: str, account: Path | None = None) -> Path | None:
    base = home(name, account)
    return None if base is None else base / ".delfin" / "messenger_sessions"


def pairing_dir(name: str, account: Path | None = None) -> Path | None:
    base = home(name, account)
    return None if base is None else base / ".delfin" / "messenger_pairing"


def workspaces_root(name: str, owner_root: str, account: Path | None = None) -> Path | None:
    """Where ``neu NAME`` puts a person's session folders: the owner's
    configured root (DELFIN_MESSENGER_WORKSPACES), or ``<home>/DELFIN-Handy``
    of the person when the owner's root is configured at all."""
    if not owner_root:
        return None
    if name == OWNER:
        return Path(owner_root)
    base = home(name, account)
    return None if base is None else base / "DELFIN-Handy"
