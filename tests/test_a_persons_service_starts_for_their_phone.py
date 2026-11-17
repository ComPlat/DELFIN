"""A person's phone works while no window of theirs is open.

Input: a person's agent service is down (it exits when idle, and is
otherwise started by their TUI/dashboard through
delfin.agent.service.launcher). Before this, every phone command of that
person ended in "service not running".

Output:

* the bridge starts the person's service through the launcher's own start
  command (launcher.service_command) with that person's HOME, runtime dir,
  DELFIN_PERSON and settings defaults -- the people-mode environment of
  tools/windows-launcher/remote_bootstrap.sh -- and the phone's
  ``sessions`` gets its answer;
* each start is audited (messenger_service_started); starts are
  rate-limited per person (a second one within a minute is refused and
  audited);
* a name that is not a real, owned ~/people/<name> directory is refused
  before anything runs;
* on a host with the starter unit (DELFIN_MESSENGER_PERSON_STARTS) the
  bridge only leaves a request file, and the starter runs the same start.
"""
from __future__ import annotations

import os

import pytest

from _messenger_people_rig import Account
from delfin.agent import security_events
from delfin.agent.service import launcher
from delfin.messenger import people, store


@pytest.fixture(autouse=True)
def _german_phone(monkeypatch):
    """These tests read the German phone texts (lang.DE); English is the
    default (test_the_phone_speaks_the_persons_language)."""
    monkeypatch.setenv("DELFIN_MESSENGER_LANG", "de")


@pytest.fixture
def acct(tmp_path, monkeypatch):
    security_events.clear()
    account = Account(tmp_path, monkeypatch)
    store.save("21", person="anna", label="anna phone")
    sub = account.bridge._person_bridge("anna")
    sub._service = None
    account.up = False
    account.starts = []

    def connect():
        if not account.up:
            raise FileNotFoundError("no socket")
        sub._service = account.services["anna"]
        return sub._service

    monkeypatch.setattr(sub, "_connect_service", connect)

    def start(name, **_):
        account.starts.append(name)
        account.up = True
        return True

    monkeypatch.setattr(people, "start_service", start)
    monkeypatch.setattr(people, "START_WAIT_S", 0.2)
    yield account
    security_events.clear()


def _kinds():
    return [e.kind for e in security_events.recent(100)]


def _settle(acct, text):
    """One phone message, then the gateway's answer after the start."""
    import asyncio

    async def main():
        await acct.rooms.handle({"chat": "@21", "sender": "21", "text": text,
                                 "id": str(len(acct.daemon.sent) + 1)})
        await asyncio.wait_for(asyncio.gather(*list(acct.rooms._retries)), 10)
    asyncio.run(main())


def test_a_phone_command_starts_the_persons_service_and_gets_its_answer(acct):
    _settle(acct, "sessions")
    assert acct.starts == ["anna"]
    assert "1. Butan" in acct.daemon.to("@21")[-1]
    assert "messenger_service_started" in _kinds()


def test_starts_are_rate_limited_per_person(acct):
    _settle(acct, "sessions")
    acct.up = False
    acct.bridge._person_bridge("anna")._service = None
    _settle(acct, "sessions")
    assert acct.starts == ["anna"]
    assert "läuft gerade nicht" in acct.daemon.to("@21")[-1]
    assert "messenger_service_start_refused" in _kinds()


def test_the_start_runs_the_launchers_command_in_the_persons_environment(tmp_path, monkeypatch):
    account = Account(tmp_path, monkeypatch)
    seen = {}

    def run(cmd, **kwargs):
        seen.update(cmd=cmd, **kwargs)

        class Done:
            returncode = 0
        return Done()

    monkeypatch.delenv(people.START_REQUESTS_ENV, raising=False)
    assert people.start_service("anna", run=run) is True
    home = account.home / "people" / "anna"
    env = seen["env"]
    assert env["HOME"] == str(home) and env["DELFIN_PERSON"] == "anna"
    assert env["DELFIN_RUNTIME_DIR"] == str(home / ".run")
    assert env["DELFIN_SETTINGS_DEFAULTS"] == str(account.home / ".delfin_settings.json")
    assert not any(k.startswith("DELFIN_MESSENGER_") for k in env)
    assert seen["cmd"][-len(launcher.service_command()):] == launcher.service_command()
    assert seen["cwd"] == str(home / ".run")
    assert seen["stdin"] is not None and seen.get("start_new_session") is True


@pytest.mark.parametrize("name", ["", "ghost", "../anna", "link"])
def test_a_start_for_anything_but_a_real_person_is_refused(tmp_path, monkeypatch, name):
    account = Account(tmp_path, monkeypatch)
    os.symlink(account.home / "people" / "anna", account.home / "people" / "link")
    ran = []
    assert people.start_service(name, run=lambda *a, **k: ran.append(a)) is False
    assert not ran


def test_with_the_starter_unit_the_bridge_leaves_a_request_and_the_starter_starts(tmp_path, monkeypatch):
    account = Account(tmp_path, monkeypatch)
    requests = tmp_path / "starts"
    requests.mkdir(mode=0o700)
    monkeypatch.setenv(people.START_REQUESTS_ENV, str(requests))
    ran = []
    assert people.request_start("anna", run=lambda cmd, **k: ran.append(cmd)) is True
    assert not ran and (requests / "start-anna").exists()
    (requests / "start-..%2Fx").write_text("")
    started = []
    people.start_requested(requests, start=lambda name: started.append(name) or True)
    assert started == ["anna"] and not list(requests.iterdir())
