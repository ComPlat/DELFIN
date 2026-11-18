"""Named accounts of the dashboard, beside the single token.

Each account is a name and its own secret, stored only as a hash in the
format jupyter-server's own authentication already relies on
(:func:`jupyter_server.auth.passwd` / :func:`passwd_check`, argon2 by
default) -- the store never holds a secret in clear, and no secret or
hash reaches a log line, an error message or the page.

The accounts file is a JSON list of
``{"name": ..., "secret_hash": ..., "head": true|false}`` objects; its
path is ``DELFIN_DASHBOARD_ACCOUNTS``. An empty or absent file means no
accounts: the dashboard then behaves exactly as without this module
(jupyter-server's single token). Exactly one account must be flagged as
head, and no name may be repeated; anything else is refused with
:class:`AccountsError` -- the server starts without accounts rather than
with an accounts file that cannot be audited.

This module does not touch who may reach the dashboard over the network:
the server keeps binding to the loopback.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Dict, List, Optional

#: Environment variable holding the path of the accounts file.
ACCOUNTS_ENV = "DELFIN_DASHBOARD_ACCOUNTS"

#: The field names of one account entry.
_FIELDS = ("name", "secret_hash", "head")


class AccountsError(Exception):
    """The accounts file cannot be served: it is malformed or ambiguous."""


def hash_secret(secret: str) -> str:
    """A secret as jupyter-server's hash format (never the secret itself)."""
    from jupyter_server.auth import passwd

    return passwd(secret)


def load_path() -> Optional[Path]:
    """The configured accounts file, or ``None`` when none is configured."""
    value = os.environ.get(ACCOUNTS_ENV, "")
    if not value.strip():
        return None
    return Path(value)


class Accounts:
    """The named accounts of one dashboard server.

    Truthy exactly when accounts are configured, so callers can branch on
    ``if accounts:`` to keep the single-token path untouched.
    """

    def __init__(self, source: Any = None):
        self._by_name: Dict[str, Dict[str, Any]] = {}
        self._head: Optional[str] = None
        path = source if source is not None else load_path()
        if path is None:
            return
        text = Path(path).read_text(encoding="utf-8")
        if not text.strip():
            return
        try:
            entries = json.loads(text)
        except ValueError as exc:
            raise AccountsError(f"{ACCOUNTS_ENV}: not valid JSON") from exc
        if not isinstance(entries, list):
            raise AccountsError(f"{ACCOUNTS_ENV}: expected a list of accounts")
        for entry in entries:
            if not isinstance(entry, dict) \
                    or any(name not in entry for name in _FIELDS):
                raise AccountsError(
                    f"{ACCOUNTS_ENV}: every account needs the fields "
                    + ", ".join(_FIELDS))
            name = entry["name"]
            if not isinstance(name, str) or not name or len(name) > 64 \
                    or not all(ch.isprintable() for ch in name):
                raise AccountsError(
                    f"{ACCOUNTS_ENV}: an account name must be a short, "
                    "printable, non-empty string")
            if name in self._by_name:
                raise AccountsError(f"{ACCOUNTS_ENV}: the name is repeated")
            if not isinstance(entry["secret_hash"], str) \
                    or ":" not in entry["secret_hash"]:
                raise AccountsError(
                    f"{ACCOUNTS_ENV}: account {name}: secret_hash must be a "
                    "jupyter-server hash (type:salt:hash or argon2:hash)")
            self._by_name[name] = {"secret_hash": entry["secret_hash"]}
            if entry["head"]:
                if self._head is not None:
                    raise AccountsError(
                        f"{ACCOUNTS_ENV}: two accounts are flagged as head")
                self._head = name
        if self._head is None:
            raise AccountsError(
                f"{ACCOUNTS_ENV}: no account is flagged as head")

    def __bool__(self) -> bool:
        """True when accounts are configured."""
        return bool(self._by_name)

    def names(self) -> List[str]:
        """The account names, in file order."""
        return list(self._by_name)

    def head(self) -> Optional[str]:
        """The head account's name, or ``None`` with no accounts."""
        return self._head

    def is_head(self, name: str) -> bool:
        return name == self._head

    def verify(self, name: str, secret: str) -> bool:
        """True when *secret* is the secret of account *name*.

        Never raises for a wrong name or secret; an unusable hash counts
        as a mismatch, not as a server error.
        """
        entry = self._by_name.get(str(name))
        if entry is None:
            return False
        try:
            from jupyter_server.auth import passwd_check

            return bool(passwd_check(entry["secret_hash"], str(secret)))
        except Exception:
            return False
