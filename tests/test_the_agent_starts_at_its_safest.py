"""Isolation is on wherever this host can hold a command, not only in bypass.

The shipped setting ``agent.bash_isolation`` is "auto", and "auto" used to
mean: isolate in the unattended (bypass) profile, and for a locked scope,
and nowhere else. So an ordinary attended session -- the one a person
actually sits in front of -- ran its commands through a plain
``/bin/bash -c``, and the banner said so honestly:

    isolation  off — a command the agent runs can still write outside
               the workspace (--isolate)

Honest is not the same as safe. The reasoning behind the old default was
that a human approves each command in the attended modes, but approval is
given on the command TEXT, and the text is not the act: an interpreter, a
symlink or a base64 round-trip walks past any reading of it. That is the
same argument the locked-scope branch already makes two screens further
up, and it does not stop being true because somebody is watching.

So the agent starts at the safest posture this host can actually provide
-- bubblewrap, else Landlock, else Seatbelt -- and a person who needs
unrestricted writes turns it off deliberately, with ``--no-isolate`` or
the setting. Safe by default, disarmed on purpose.

What must NOT happen is a host with no mechanism refusing every command.
That would be "secure" and useless, and it is a different promise from
the one a locked scope makes, where refusing IS the right answer.
"""

from __future__ import annotations

import pytest

from delfin.agent import api_client as A


@pytest.fixture()
def perms(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    return A.KitToolPermissions(workspace=ws, mode="default")


@pytest.fixture(autouse=True)
def _no_override():
    A.set_bash_isolation_override("")
    yield
    A.set_bash_isolation_override("")


def _argv(perms, tmp_path, monkeypatch, *, bwrap=False, landlock=False,
          seatbelt=False):
    # shutil.which too: the branch below the decision refuses a forced
    # "bwrap" when the binary is absent, so a fake that says "functional"
    # while the binary is missing tests the refusal, not the choice.
    monkeypatch.setattr(A.shutil, "which",
                        lambda name, *a, **k: "/usr/bin/bwrap" if (
                            name == "bwrap" and bwrap) else None)
    monkeypatch.setattr(A, "_bwrap_functional", lambda: bwrap)
    monkeypatch.setattr(A, "_landlock_functional", lambda: landlock)
    monkeypatch.setattr(A, "_seatbelt_functional", lambda: seatbelt)
    return A._bash_isolation_argv("echo hi", str(tmp_path), perms)


def test_an_attended_session_is_isolated_where_bwrap_works(
        perms, tmp_path, monkeypatch):
    argv = _argv(perms, tmp_path, monkeypatch, bwrap=True)
    assert any("bwrap" in str(part) for part in argv), argv


def test_landlock_is_taken_when_there_is_no_bwrap(
        perms, tmp_path, monkeypatch):
    argv = _argv(perms, tmp_path, monkeypatch, landlock=True)
    joined = " ".join(str(p) for p in argv)
    assert "landlock" in joined.lower(), argv


def test_seatbelt_is_taken_on_a_mac(perms, tmp_path, monkeypatch):
    argv = _argv(perms, tmp_path, monkeypatch, seatbelt=True)
    joined = " ".join(str(p) for p in argv)
    assert "sandbox-exec" in joined or "seatbelt" in joined.lower(), argv


def test_a_host_with_nothing_still_runs_the_command(
        perms, tmp_path, monkeypatch):
    """Secure and useless is not a posture. The banner is what says the
    host cannot hold it; refusing here would be a different promise from
    the one a locked scope makes."""
    argv = _argv(perms, tmp_path, monkeypatch)
    joined = " ".join(str(p) for p in argv)
    assert "bash" in joined, argv
    assert "refus" not in joined.lower(), argv


def test_it_can_be_disarmed(perms, tmp_path, monkeypatch):
    A.set_bash_isolation_override("off")
    argv = _argv(perms, tmp_path, monkeypatch, bwrap=True)
    assert not any("bwrap" in str(part) for part in argv), (
        "'off' is the escape hatch for a setup that needs unrestricted "
        "writes, and it has to actually let go")


def test_the_cli_offers_the_disarm():
    """A weakness a user cannot act on is half the truth; so is a
    protection they cannot switch off when it blocks real work."""
    import inspect

    from delfin.agent import cli
    src = inspect.getsource(cli.build_parser)
    assert '"--no-isolate"' in src, "there must be a way to turn it off"
