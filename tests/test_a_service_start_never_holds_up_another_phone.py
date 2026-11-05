"""Starting one person's agent service never holds up anybody else's phone.

Review finding 5 (DoS): the start of a person's service ran on the
bridge's event loop -- time.sleep up to 30 s, subprocess.run up to 90 s --
so every other phone waited behind one start.

Output:

* the start runs in a thread of its own; the command that needed it is
  answered at once with "wird gestartet" (``starting``), and the gateway
  answers the phone again once the service is up;
* while one start hangs, another person's phone is served within 0.5 s
  (control: the same command for that person with no start going);
* a second command for the starting person does not start a second time.
"""
from __future__ import annotations

import asyncio
import threading
import time

import pytest

from _messenger_people_rig import Account
from delfin.messenger import people, store


@pytest.fixture
def acct(tmp_path, monkeypatch):
    account = Account(tmp_path, monkeypatch)
    store.save("21", person="anna", label="anna phone")
    store.save("31", person="ben", label="ben phone")
    sub = account.bridge._person_bridge("anna")
    sub._service = None
    account.up = threading.Event()
    account.release = threading.Event()
    account.starts = []

    def connect():
        if not account.up.is_set():
            raise FileNotFoundError("no socket")
        sub._service = account.services["anna"]
        return sub._service

    monkeypatch.setattr(sub, "_connect_service", connect)

    def start(name, **_):
        account.starts.append(name)
        account.release.wait(10)          # a start that hangs until released
        account.up.set()
        return True

    monkeypatch.setattr(people, "start_service", start)
    monkeypatch.setattr(people, "START_WAIT_S", 5.0)
    yield account
    account.release.set()


def test_another_persons_phone_is_served_while_a_start_hangs(acct):
    async def main():
        t0 = time.monotonic()
        first = asyncio.create_task(acct.bridge._run_command({"cmd": "sessions", "person": "anna"}))
        await asyncio.sleep(0.05)
        ben = await acct.bridge._run_command({"cmd": "sessions", "person": "ben"})
        ben_at = time.monotonic() - t0
        anna = await first
        anna_at = time.monotonic() - t0
        again = await acct.bridge._run_command({"cmd": "sessions", "person": "anna"})
        return ben, ben_at, anna, anna_at, again

    ben, ben_at, anna, anna_at, again = asyncio.run(main())
    assert ben.get("ok") and ben_at < 0.5
    assert anna.get("starting") and anna_at < 0.5
    assert again.get("starting") and acct.starts == ["anna"]
    acct.release.set()


def test_the_phone_is_told_starting_then_gets_its_answer(acct):
    async def main():
        await acct.rooms.handle({"chat": "@21", "sender": "21", "text": "sessions", "id": "1"})
        told = list(acct.daemon.to("@21"))
        acct.release.set()
        await asyncio.wait_for(asyncio.gather(*list(acct.rooms._retries)), 10)
        return told

    told = asyncio.run(main())
    assert "wird gestartet" in told[-1]
    assert "1. Butan" in acct.daemon.to("@21")[-1]
