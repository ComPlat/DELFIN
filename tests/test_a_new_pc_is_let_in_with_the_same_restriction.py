"""A new person's PC is let in under the restriction the launcher uses.

On a server that accepts passwords the launcher registers a PC by itself:
the first login uses the password and ``remote_bootstrap.sh`` appends the
PC's key to ``authorized_keys``, restricted to what the launcher actually
needs. Every later login uses the key.

Where ``PasswordAuthentication no`` is set -- a hardened box, which is
what the lab Raspberry Pi is; measured 2026-11-19 -- that first login
cannot happen. The key is added once you are in, and you get in with the
key. A new person is locked out, and nothing in the tree did that one
step: ``delfin-people`` was read-only.

So it was done by hand, and a hand-written ``ssh-ed25519 ...`` line is an
UNRESTRICTED key: agent forwarding, X11 and ``~/.ssh/rc`` all come with
it, none of which the launcher uses. That is the failure this guards
against -- not a missing feature, a quietly wider grant.

The check and the prefix here are the bootstrap's own, so the two routes
cannot accept different things.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "pi" / "delfin_people.sh"
KEY = ("ssh-ed25519 "
       "AAAAC3NzaC1lZDI1NTE5AAAAILc42IKut51Wm5lcpANYgPVi5I0xCBfgkBRQjeSBpnA7"
       " delfin-launcher-example")


def _run(home: Path, *args: str, stdin: str = "") -> subprocess.CompletedProcess:
    return subprocess.run(
        ["bash", str(SCRIPT), *args], input=stdin, capture_output=True,
        text=True, env={"HOME": str(home), "PATH": "/usr/bin:/bin"})


@pytest.fixture
def home(tmp_path):
    (tmp_path / ".ssh").mkdir()
    return tmp_path


def _keys(home: Path) -> str:
    path = home / ".ssh" / "authorized_keys"
    return path.read_text(encoding="utf-8") if path.exists() else ""


def test_the_key_is_added_restricted_not_bare(home):
    result = _run(home, "--add-key", "simone", stdin=KEY + "\n")

    assert result.returncode == 0, result.stderr
    line = _keys(home).strip()
    assert line.startswith("restrict,port-forwarding,pty "), (
        "a bare key line carries agent forwarding, X11 and ~/.ssh/rc; the "
        "launcher uses none of them. " + line)
    assert "delfin-launcher-simone" in line, "say whose PC it is, to revoke it"


def test_the_file_is_private(home):
    _run(home, "--add-key", "simone", stdin=KEY + "\n")

    assert (home / ".ssh" / "authorized_keys").stat().st_mode & 0o077 == 0


def test_it_works_when_there_is_no_ssh_directory_yet(tmp_path):
    result = _run(tmp_path, "--add-key", "simone", stdin=KEY + "\n")

    assert result.returncode == 0, result.stderr
    assert (tmp_path / ".ssh" / "authorized_keys").exists()
    assert (tmp_path / ".ssh").stat().st_mode & 0o077 == 0


def test_the_same_key_is_not_added_twice(home):
    _run(home, "--add-key", "simone", stdin=KEY + "\n")
    _run(home, "--add-key", "simone", stdin=KEY + "\n")

    assert _keys(home).count("AAAAC3Nza") == 1


def test_anything_that_is_not_one_ed25519_line_is_refused(home):
    for bad, why in (
            ("ssh-rsa AAAAB3Nza... someone", "another key type"),
            ("not a key at all", "prose"),
            ("", "nothing"),
            (KEY + "\nssh-ed25519 AAAAsecond extra", "a second line"),
            ("restrict,port-forwarding,pty " + KEY, "an options prefix"),
    ):
        _run(home, "--add-key", "simone", stdin=bad + "\n")
        assert "AAAAsecond" not in _keys(home), why

    assert "ssh-rsa" not in _keys(home)
    assert "not a key" not in _keys(home)


def test_a_second_line_cannot_ride_in_on_the_first(home):
    """The bootstrap refuses a multi-line value for this reason: the check
    reads a line, so a second one would be appended unchecked."""
    _run(home, "--add-key", "simone",
         stdin=KEY + "\nssh-ed25519 AAAAsnuck in-unchecked\n")

    assert "snuck" not in _keys(home)


def test_a_bad_person_name_adds_nothing(home):
    for name in ("../../root", "Simone Meier", "-rf", ""):
        result = _run(home, "--add-key", name, stdin=KEY + "\n")
        assert result.returncode != 0, name

    assert _keys(home) == ""


def test_an_account_that_refuses_launcher_keys_is_honoured(home):
    (home / ".config" / "delfin").mkdir(parents=True)
    (home / ".config" / "delfin" / "launcher-no-keys.conf").write_text("")

    result = _run(home, "--add-key", "simone", stdin=KEY + "\n")

    assert result.returncode != 0
    assert _keys(home) == "", (
        "the account said no launcher keys; this route must not be the way "
        "around that")


def test_the_listing_names_a_key_that_carries_no_restriction(home):
    """What the listing is for: an unrestricted key already in the file is
    the thing you want to find, so it must not be shown as ordinary."""
    (home / ".ssh" / "authorized_keys").write_text(
        f"restrict,port-forwarding,pty {KEY}\n"
        "ssh-ed25519 AAAAolder someone-elses-pc\n", encoding="utf-8")

    result = _run(home, "--keys")

    assert result.returncode == 0, result.stderr
    assert "NONE - full access" in result.stdout, result.stdout
    assert "someone-elses-pc" in result.stdout
    assert "delfin-launcher-example" in result.stdout
