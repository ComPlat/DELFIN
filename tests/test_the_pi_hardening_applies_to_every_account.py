"""The Pi's SSH hardening is global, not scoped to one account.

Measured on a live box on 2026-11-18. ``/etc/ssh/sshd_config.d`` held,
in sorted order:

    30-temp-password.conf   Match User delfin
                                PasswordAuthentication yes
    40-delfin-hardening.conf    PasswordAuthentication no
                                PermitRootLogin no
                                MaxAuthTries 3
                                AllowUsers pi delfin

sshd reads the drop-ins in sorted order and concatenates them, and a
``Match`` block stays open across a file boundary until the next
``Match``. So every directive of the hardening file sat INSIDE
``Match User delfin``: the whole block, ``AllowUsers`` included, applied
to that one account, and every other account on the machine kept the
OpenSSH defaults -- six auth tries, a two-minute grace, no allow-list.

``pi_harden_host.sh`` already removed a password exception at ``60-``,
which sorts after the hardening and would simply override it. The one
that does the damage sorts BEFORE it, which the script did not consider.

The fix is one line, ``Match all`` at the top of the hardening file, and
these tests pin the property rather than the line: the check below
composes the drop-ins the way sshd does and asserts which directives are
global. A test that only grepped for the string would pass for a file
that happens to contain it somewhere harmless.

The earlier exception is NOT removed by the script: it is how a new
person's PC registers its key (the launcher sends the key on the first
password login). It is reported instead, because an exception nobody can
see is an exception nobody withdraws.
"""

from __future__ import annotations

import re
from pathlib import Path

HARDEN = Path(__file__).resolve().parents[1] / "scripts" / "pi" / \
    "pi_harden_host.sh"

# The onboarding exception a people-mode setup writes, and which the live
# box had.
TEMP_PASSWORD = "Match User delfin\n    PasswordAuthentication yes\n"


def _hardening_body() -> str:
    """The heredoc the script writes to 40-delfin-hardening.conf."""
    text = HARDEN.read_text(encoding="utf-8")
    m = re.search(
        r"cat > /etc/ssh/sshd_config\.d/40-delfin-hardening\.conf <<EOF\n"
        r"(.*?)\nEOF\n", text, re.S)
    assert m, "the hardening heredoc is no longer where this test reads it"
    return m.group(1) + "\n"


def _scope_of_each_directive(conf: str) -> dict[str, str]:
    """Which Match block each directive falls in, as sshd parses it.

    Returns {directive_lowercased: "global" | the Match criteria}. The
    first occurrence wins, which is sshd's own rule ("the first obtained
    value for each parameter is used").
    """
    scope = "global"
    out: dict[str, str] = {}
    for raw in conf.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        word, _, rest = line.partition(" ")
        if word.lower() == "match":
            crit = rest.strip()
            scope = "global" if crit.lower() == "all" else crit
            continue
        key = word.lower()
        if key not in out:
            out[key] = scope
    return out


def _composed(*files: str) -> str:
    """What sshd reads: the drop-ins concatenated in sorted order."""
    return "".join(files)


def test_the_hardening_is_global_behind_an_earlier_match_block():
    conf = _composed(TEMP_PASSWORD, _hardening_body())
    scope = _scope_of_each_directive(conf)

    for directive in ("permitrootlogin", "maxauthtries", "logingracetime",
                      "x11forwarding", "kbdinteractiveauthentication",
                      "allowusers"):
        assert scope.get(directive) == "global", (
            f"{directive} must apply to every account, not only to "
            f"{scope.get(directive)!r} -- an earlier Match block in a "
            "lower-numbered drop-in would otherwise swallow it")


def test_the_onboarding_exception_still_wins_for_its_own_account():
    """The password window is deliberate: it is how a new PC registers its
    key. sshd keeps the FIRST value, so the exception must survive the
    hardening that follows it."""
    conf = _composed(TEMP_PASSWORD, _hardening_body())
    scope = _scope_of_each_directive(conf)
    assert scope.get("passwordauthentication") == "User delfin", (
        "the first PasswordAuthentication read is the exception's; the "
        "hardening's later 'no' must not become the effective value for "
        "that account, or no new PC could ever register its key")


def test_the_hardening_alone_is_global_too():
    """Without any exception file, nothing changes for the ordinary box."""
    scope = _scope_of_each_directive(_hardening_body())
    assert scope.get("passwordauthentication") == "global"
    assert scope.get("allowusers") == "global"


def test_match_all_is_the_first_directive_of_the_file():
    """Not a style point: a directive placed above it would still be
    captured by an earlier block."""
    body = [ln.strip() for ln in _hardening_body().splitlines()
            if ln.strip() and not ln.strip().startswith("#")]
    assert body[0].lower() == "match all", body[:3]


def test_the_script_names_an_exception_it_does_not_remove():
    text = HARDEN.read_text(encoding="utf-8")
    assert "rm -f /etc/ssh/sshd_config.d/60-delfin-password.conf" in text, (
        "the exception that sorts AFTER the hardening is still removed")
    assert "opens a Match block" in text, (
        "an exception that sorts BEFORE it must be reported: it is needed "
        "for onboarding and must not be silent")
    assert "sudo rm" in text, "say how to withdraw it"


def test_the_script_still_validates_before_reloading():
    text = HARDEN.read_text(encoding="utf-8")
    assert text.index("sshd -t") < text.index("systemctl reload ssh"), (
        "a broken config must be caught before the reload that would "
        "apply it -- sshd is the only way back into the box")
