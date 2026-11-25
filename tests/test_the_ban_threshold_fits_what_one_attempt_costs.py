"""One click must not spend the whole ban budget.

fail2ban counts SSH authentication failures. The Windows launcher makes up
to THREE connections per click -- the script check, the copy, and the
dashboard -- and OpenSSH asks for the password once per connection. So a
single attempt with a wrong password costs three failures, not one.

At the jail's earlier ``maxretry = 5`` that meant a newcomer was banned on
their SECOND click, for an hour, with ``bantime.increment`` lengthening
every further ban. Measured on the lab Pi on 2026-11-19: 21 failures, one
banned address, and an afternoon spent looking at passwords, consoles and
clipboards, because a ban and a wrong password are the same sentence to
whoever is typing -- so the person retries, which extends the ban.

This pins the relationship between the two numbers, which live in two
files that nothing else connects: the launcher's connection count and the
jail's threshold. Raising one without the other is the defect, and it is
invisible in either file alone.
"""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HARDEN = (ROOT / "scripts" / "pi" / "pi_harden_host.sh").read_text(
    encoding="utf-8")
APP = (ROOT / "tools" / "windows-launcher" / "DELFIN.ps1").read_text(
    encoding="utf-8")

#: Connections one click makes at most, each with its own password prompt:
#: on a new launcher version the dashboard connection finds no script copy,
#: the copy (Copy-ServerScripts) follows, then the dashboard again. Every
#: later click is one connection.
CONNECTIONS_PER_CLICK = 3


def _jail_value(name: str) -> str:
    match = re.search(rf"^{name}\s*=\s*(\S+)", HARDEN, re.M)
    assert match, f"{name} is not set in the jail written by the script"
    return match.group(1)


def _seconds(value: str) -> int:
    unit = {"s": 1, "m": 60, "h": 3600, "d": 86400}
    if value[-1] in unit:
        return int(value[:-1]) * unit[value[-1]]
    return int(value)


def test_one_click_really_does_cost_three_prompts():
    """The number above, read from the launcher rather than asserted."""
    assert "Invoke-OpenSSH $scp" in APP and "Invoke-OpenSSH $ssh" not in APP, (
        "the copy is the one connection besides the dashboard; if that "
        "changed, so must CONNECTIONS_PER_CLICK")
    assert "Copy-ServerScripts ([string]$p.Target)" in APP
    assert "$scriptsCopied" in APP, "the dashboard is started again at most once"
    assert "Get-RemoteCommand $p 'dashboard'" in APP


def test_the_threshold_allows_more_than_one_attempt():
    maxretry = int(_jail_value("maxretry"))

    assert maxretry >= 3 * CONNECTIONS_PER_CLICK, (
        f"maxretry={maxretry} with {CONNECTIONS_PER_CLICK} prompts per "
        "click bans a person before they have tried three times; they "
        "cannot tell that from a refused password, so they keep trying")


def test_a_ban_is_a_pause_not_a_wall():
    """Long enough to make guessing hopeless, short enough that a person
    waits it out instead of reporting the machine as broken."""
    bantime = _seconds(_jail_value("bantime"))

    assert bantime <= 15 * 60, f"bantime={bantime}s is a wall"
    assert bantime >= 60, f"bantime={bantime}s stops nothing"


def test_a_persistent_address_is_still_punished_harder():
    """Forgiving the first ban must not forgive the thousandth."""
    assert _jail_value("bantime.increment") == "true"
    assert _jail_value("enabled") == "true"
    assert _jail_value("banaction") == "nftables-multiport"


def test_the_readme_says_to_look_at_the_jail_first():
    """The whole cost here was not the setting, it was not knowing to
    look. A ban and a wrong password read identically."""
    readme = (ROOT / "scripts" / "pi" / "README.md").read_text(
        encoding="utf-8")

    assert "looks exactly like a wrong password" in readme
    assert "fail2ban-client status sshd" in readme
    assert "unbanip" in readme
