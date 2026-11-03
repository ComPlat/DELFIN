"""Host hardening turns SSH passwords off only for a login user with a key.

Input: scripts/pi/pi_harden_host.sh. It writes PasswordAuthentication no and
removes the earlier password exception (60-delfin-password.conf). A box set
up by hand with a password and no authorized_keys for the login user was then
unreachable over SSH after the next reload.
"""
from __future__ import annotations

import pathlib

_SH = (pathlib.Path(__file__).resolve().parents[1] / "scripts" / "pi"
       / "pi_harden_host.sh").read_text(encoding="utf-8")


def test_the_key_is_checked_before_passwords_are_turned_off():
    check = _SH.find("authorized_keys")
    assert check != -1, "no check for the login user's key"
    assert check < _SH.index("PasswordAuthentication no")
    assert check < _SH.index("60-delfin-password.conf")
