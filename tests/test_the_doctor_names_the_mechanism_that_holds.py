"""What holds a command here has three names, and the doctor knew one.

``delfin-agent doctor`` asked ``_bwrap_functional()`` and nothing else.
The product does not: where bubblewrap cannot run it uses Landlock, and
on macOS Seatbelt, and a forced isolation with none of the three REFUSES
the command rather than running it unheld. So on a host without
bubblewrap the row read

    auto — isolated in bypassPermissions only; bwrap unusable here, so
    never isolated
    fix: install bwrap for real containment

on a machine that was in fact isolating every unattended command through
Landlock -- and a session configured with bash_isolation = "bwrap" was
reported as FAIL "bwrap does not work here" while it was running held.

That is a security surface reporting on the wrong mechanism, which is
the one kind of wrong report that ends with somebody trusting less, or
more, than they should. The host-level ``delfin doctor`` already knew all
of this; the two had drifted, and the last test here holds them
together.
"""

from __future__ import annotations

import pytest

from delfin.agent import doctor as AD


@pytest.fixture()
def probes(monkeypatch):
    """Set what this host can do, one mechanism at a time."""
    import delfin.agent.api_client as A

    def _set(bwrap=False, landlock=False, seatbelt=False, mode="auto"):
        monkeypatch.setattr(A, "_bwrap_functional", lambda: bwrap)
        monkeypatch.setattr(A, "_landlock_functional", lambda: landlock)
        monkeypatch.setattr(A, "_seatbelt_functional", lambda: seatbelt)
        monkeypatch.setattr(
            "delfin.user_settings.load_settings",
            lambda *a, **k: {"agent": {"bash_isolation": mode}})
    return _set


def _row(ctx=None):
    rows = AD._check_bash_isolation(ctx or {})
    assert len(rows) == 1
    return rows[0]


# -- the default setting ----------------------------------------------------

def test_bubblewrap_is_named_when_it_is_what_runs(probes):
    probes(bwrap=True)
    row = _row()
    # PASS, not WARN: under auto a host that can hold a command holds it
    # in every permission mode. The WARN here was the older resolver's.
    assert row["status"] == "PASS"
    assert "bwrap" in row["detail"]


def test_landlock_is_named_where_bubblewrap_cannot_run(probes):
    """A cluster login node rarely allows the user namespace bwrap needs,
    and offers Landlock on the same kernel. The old row called that host
    unprotected."""
    probes(bwrap=False, landlock=True)
    row = _row()
    assert "Landlock" in row["detail"]
    assert "never isolated" not in row["detail"]


def test_seatbelt_is_named_on_a_mac(probes):
    probes(bwrap=False, landlock=False, seatbelt=True)
    assert "Seatbelt" in _row()["detail"]


def test_a_host_with_none_of_them_is_told_so(probes):
    probes()
    row = _row()
    assert "nothing here can isolate" in row["detail"]
    assert "Landlock" in row["fix"], (
        "the way out is not only bubblewrap; a newer kernel does it too")


# -- isolation switched on by name ------------------------------------------

def test_a_forced_isolation_passes_on_landlock(probes):
    """bash_isolation = "bwrap" means "hold every command". Landlock does
    that where bubblewrap is absent, and the product uses it; reporting
    FAIL there told the user their session was unprotected when it was
    not."""
    probes(bwrap=False, landlock=True, mode="bwrap")
    row = _row()
    assert row["status"] == "PASS"
    assert "Landlock" in row["detail"]


def test_a_forced_isolation_without_any_mechanism_fails(probes):
    probes(mode="bwrap")
    row = _row()
    assert row["status"] == "FAIL"
    assert "refused" in row["detail"], (
        "the product does not run the command unheld -- it refuses it, "
        "and a user reading this row needs to know which of the two")


def test_switching_it_off_still_says_so(probes):
    probes(bwrap=True, mode="off")
    row = _row()
    assert row["status"] == "WARN"
    assert "explicitly off" in row["detail"]


# -- the two doctors do not contradict each other ---------------------------

@pytest.mark.parametrize("mechanism", ["bwrap", "landlock", "none"])
def test_the_host_doctor_and_the_agent_doctor_agree(probes, monkeypatch,
                                                    mechanism):
    """Both answer "what holds a command here". They are separate
    functions in separate modules, and they drift.

    This asserted only that the agent doctor does not DENY isolation on a
    host the host doctor calls isolated. That is a narrower question than
    the one it is named for, and the next drift went straight through it:
    the agent doctor called an isolated host isolated "in
    bypassPermissions only" and warned about it, naming a mechanism all
    the while. Both halves of the old assertion held while the two
    doctors gave the user opposite impressions of the same machine.

    So the VERDICTS are compared, in both directions and for each
    mechanism. The wordings may differ -- they are written for different
    readers -- but a machine cannot be fine in one report and a warning
    in the other.
    """
    from delfin import doctor as HD
    import delfin.agent.socket_guard as SG

    probes(bwrap=(mechanism == "bwrap"), landlock=(mechanism == "landlock"))
    # The host doctor's OK needs Landlock AND the socket guard; the agent
    # doctor's needs only a filesystem mechanism. Supplied, so the
    # comparison is between the two verdicts and not between two hosts.
    monkeypatch.setattr(SG, "available", lambda: True)

    host = HD.check_command_isolation()
    row = _row()
    isolated = host.status == HD.OK

    assert (row["status"] == "PASS") == isolated, (
        f"with {mechanism}, the host doctor says {host.status} and the "
        f"agent doctor says {row['status']}: one report calls this "
        f"machine fine and the other does not.\n"
        f"  host:  {host.detail}\n  agent: {row['detail']}")
    if isolated:
        assert "nothing here can isolate" not in row["detail"], row
        assert "Landlock" in row["detail"] or "bwrap" in row["detail"], row
