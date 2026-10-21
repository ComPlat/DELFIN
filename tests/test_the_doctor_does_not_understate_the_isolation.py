"""The doctor's isolation row agrees with the resolver that decides.

The row said `auto — isolated in bypassPermissions only` and offered
`agent.bash_isolation = "bwrap"` as the fix. That was true of an older
resolver: `auto` used to wall the unattended mode and a locked scope and
nothing else.

`_bash_isolation_argv` changed -- an approval is given on the command
TEXT, and the text is not the act, so an interpreter, a symlink or a
base64 round trip walks past it whether or not somebody is watching --
and the row did not. On a host with bubblewrap, all four of default,
acceptEdits, plan and bypassPermissions come back walled.

Understating protection is not the harmless direction: it tells the user
attended sessions are unguarded when they are not, and sends them to
change a setting already in force.

So the row is asserted AGAINST the resolver rather than against a
sentence, in both directions -- a host that can isolate and one that
cannot. That is the pin: the next change to the resolver has to bring
this row with it.
"""

from __future__ import annotations

import pytest

from delfin.agent import api_client as A
from delfin.agent import doctor as D


_ATTENDED = ("default", "acceptEdits", "plan")
_MODES = _ATTENDED + ("bypassPermissions",)


@pytest.fixture
def host(monkeypatch, tmp_path):
    """Everything the resolver and the doctor ask about, answered here."""
    def _supply(*, isolation="auto", mechanism="bwrap"):
        import delfin.user_settings as us
        monkeypatch.setattr(A, "_BASH_ISOLATION_OVERRIDE", "")
        monkeypatch.setattr(
            us, "load_settings",
            lambda *a, **k: {"agent": {"bash_isolation": isolation}})
        monkeypatch.setattr(A.shutil, "which",
                            lambda name: "/usr/bin/bwrap"
                            if mechanism == "bwrap" else None)
        for name, which in (("_bwrap_functional", "bwrap"),
                            ("_landlock_functional", "landlock"),
                            ("_seatbelt_functional", "seatbelt")):
            monkeypatch.setattr(A, name,
                                (lambda w=which: lambda: mechanism == w)())
        monkeypatch.setattr(A, "_process_cage_functional", lambda: False)
        monkeypatch.setattr(A, "_PROCESS_CAGE_FUNCTIONAL", False)
        monkeypatch.delenv(A._PROCESS_CAGE_ENV, raising=False)
        return tmp_path
    return _supply


def _walled(mode: str, workspace) -> bool:
    """What the resolver actually builds for *mode*."""
    perms = A.KitToolPermissions(workspace=workspace, mode=mode)
    parts = [str(p) for p in
             A._bash_isolation_argv("true", str(workspace), perms)]
    if "--ro-bind" in parts:
        return True
    if "--fs" in parts:
        i = parts.index("--fs")
        return i + 1 < len(parts) and parts[i + 1] != "0"
    return False


def _row(ctx=None):
    rows = D._check_bash_isolation(ctx or {})
    assert len(rows) == 1, rows
    return rows[0]


class TestWhereTheHostCanIsolate:
    @pytest.mark.parametrize("mode", _MODES)
    def test_every_permission_mode_is_walled(self, host, mode):
        """The fact the row has to report. Asserted first, because if
        this ever stops being true the row below must change with it."""
        ws = host(isolation="auto", mechanism="bwrap")
        assert _walled(mode, ws), f"{mode} is not confined under auto"

    def test_the_row_says_every_mode_and_passes(self, host):
        host(isolation="auto", mechanism="bwrap")
        row = _row()
        assert row["status"] == "PASS", row
        assert "every permission mode" in row["detail"], row
        assert "bypassPermissions only" not in row["detail"]

    def test_it_does_not_propose_a_setting_already_in_force(self, host):
        """The old row's fix was to switch on what was already on."""
        host(isolation="auto", mechanism="bwrap")
        row = _row()
        assert not row.get("fix"), row
        assert not row.get("setting"), row

    def test_landlock_is_named_where_it_is_what_holds(self, host):
        """A cluster login node rarely allows the user namespace bwrap
        needs, and the same kernel usually offers Landlock."""
        host(isolation="auto", mechanism="landlock")
        row = _row()
        assert row["status"] == "PASS", row
        assert "Landlock" in row["detail"], row


class TestWhereItCannot:
    @pytest.mark.parametrize("mode", _ATTENDED)
    def test_nothing_walls_an_attended_command(self, host, mode):
        ws = host(isolation="auto", mechanism="none")
        assert not _walled(mode, ws)

    def test_the_row_warns_and_names_what_is_left(self, host):
        host(isolation="auto", mechanism="none")
        row = _row()
        assert row["status"] == "WARN", row
        assert "nothing here can isolate" in row["detail"], row
        # Which protection is actually holding, not only that one is not.
        assert "write-target gate" in row["detail"], row
        assert "socket guard" in row["detail"], row

    def test_it_does_not_propose_bwrap_with_no_bwrap(self, host):
        """Proposing it would turn a warning into a refusal of every
        shell command, which is worse than the state being fixed."""
        host(isolation="auto", mechanism="none")
        row = _row()
        assert not row.get("setting"), row
        assert "install bubblewrap" in str(row.get("fix") or ""), row


class TestTheOtherSettings:
    def test_off_still_warns(self, host):
        host(isolation="off", mechanism="bwrap")
        row = _row()
        assert row["status"] == "WARN"
        assert "explicitly off" in row["detail"]

    def test_bwrap_with_a_mechanism_passes(self, host):
        host(isolation="bwrap", mechanism="bwrap")
        assert _row()["status"] == "PASS"

    def test_bwrap_with_none_fails_and_says_commands_are_refused(self, host):
        host(isolation="bwrap", mechanism="none")
        row = _row()
        assert row["status"] == "FAIL"
        assert "refused" in row["detail"]
