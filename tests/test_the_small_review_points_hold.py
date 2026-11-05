"""Review finding 10: three small points.

* store.expire() (every 5 s) rewrote pair.json each time -- SD-card wear.
  Now a store write happens only on a change (control: an expiry that
  removes something writes).
* ``revoke owner`` meant the account owner, so a person named ``owner``
  could not be revoked by name. The owner is ``revoke --owner`` now.
* delfin-person-start.service ran the person's service inside its own
  oneshot cgroup, unhardened. The starter is hardened (ProtectSystem,
  ProtectHome, PrivateTmp, RestrictSUIDSGID, LockPersonality, MemoryMax,
  NoNewPrivileges) and starts each person's service as a transient unit of
  its own (systemd-run --uid, its own cgroup) -- the supervisor, as the
  account's delfin-agent-service unit runs it, in the person's environment.
"""
from __future__ import annotations

import json
import os
import time

import pytest

from delfin.messenger import pairing, people, store


@pytest.fixture
def state(tmp_path, monkeypatch):
    monkeypatch.setenv("DELFIN_MESSENGER_DIR", str(tmp_path / "pairing"))
    monkeypatch.setattr(store, "_REFUSED_SEEN", set())
    return tmp_path / "pairing"


def test_an_expiry_with_nothing_to_expire_does_not_write(state):
    store.save("7")
    before = os.stat(state / "pair.json")
    time.sleep(0.01)
    for _ in range(3):
        store.expire()
    after = os.stat(state / "pair.json")
    assert (after.st_ino, after.st_mtime_ns) == (before.st_ino, before.st_mtime_ns)
    # control: an invitation past its time is removed, and that is written
    store.add_invitation("41", person="anna", request="r" * 32, now=1000.0)
    mid = os.stat(state / "pair.json")
    store.expire(now=1000.0 + pairing.INVITATION_TTL_S + 1)
    assert os.stat(state / "pair.json").st_ino != mid.st_ino
    assert not store.invitations()


def test_a_person_named_owner_is_revocable_by_name(state):
    from delfin.messenger import __main__ as admin
    store.save("7")                              # the account owner
    store.save("40", person="owner")             # a person called "owner"
    assert admin.main(["revoke", "owner"]) == 0
    assert store.contact("40") is None and store.contact("7") is not None
    assert admin.main(["revoke", "--owner"]) == 0
    assert store.contact("7") is None


def test_the_starter_unit_is_hardened(tmp_path):
    from test_the_messenger_provisioning_covers_every_pinned_server import run_script
    done = run_script(tmp_path)
    assert done.returncode == 0, done.stderr
    unit = (tmp_path / "units" / "delfin-person-start.service").read_text()
    lines = [line.strip() for line in unit.splitlines()]
    for want in ("NoNewPrivileges=yes", "ProtectSystem=strict", "ProtectHome=read-only",
                 "PrivateTmp=yes", "RestrictSUIDSGID=yes", "LockPersonality=yes"):
        assert want in lines, want
    assert any(line.startswith("MemoryMax=") for line in lines)
    assert "KillMode=process" not in lines           # nothing of a person's service stays in it


def test_a_persons_service_starts_as_a_transient_unit_of_its_own(tmp_path, monkeypatch):
    home = tmp_path / "acct"
    (home / ".config" / "delfin").mkdir(parents=True)
    (home / ".config" / "delfin" / "launcher-people.conf").write_text("on\n")
    (home / "people" / "anna").mkdir(parents=True)
    monkeypatch.setattr(people, "account_home", lambda: home)
    argv = people.transient_argv("anna")
    assert argv[0].endswith("systemd-run") and "--collect" in argv
    assert f"--uid={os.getuid()}" in argv and "--unit=delfin-person-anna" in argv
    assert "--property=KillMode=mixed" in argv
    joined = " ".join(argv)
    assert f"--setenv=HOME={home / 'people' / 'anna'}" in argv
    assert "--setenv=DELFIN_PERSON=anna" in argv
    assert "delfin.agent.service.supervisor" in joined
    assert people.transient_argv("ghost") is None and people.transient_argv("") is None
    json.dumps(argv)
