"""The OS account this process runs as, for decisions and records.

Input: none. Output: :func:`name` -- the login name of the real user id
(``pwd.getpwuid(os.getuid()).pw_name``); the decimal uid when the
password database has no entry (a container with an unlisted uid);
``getpass.getuser()`` only where the OS has no user ids (Windows), where
it reads %USERNAME%.

Semantics: the kernel's uid cannot be changed by the process; the
environment can. ``getpass.getuser()`` returns ``$LOGNAME``, ``$USER``,
``$LNAME`` or ``$USERNAME`` first, so ``USER=<admin> delfin-agent ...``
claimed another account wherever it was used for an authorisation (the
knowledge admin check), an audit destination (``{user}`` in a managed
audit path) or a record of who answered. Those use :func:`name`;
``$USER`` stays acceptable where a name is only shown.

Rejected: ``os.getlogin()`` (the controlling terminal's login, absent
without a terminal and not the effective account under su/sudo);
caching (a test or a long process may change uid; the lookup is cheap).
"""
from __future__ import annotations

import os

__all__ = ["name"]


def name() -> str:
    getuid = getattr(os, "getuid", None)
    if getuid is None:
        try:
            import getpass
            return getpass.getuser()
        except Exception:
            return ""
    uid = getuid()
    try:
        import pwd
        return pwd.getpwuid(uid).pw_name
    except (ImportError, KeyError):
        return str(uid)
