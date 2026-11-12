"""One person's flood of pairing requests never crowds out another person's.

Review finding 6: the bridge's pairing poll returned the first 20 requests
over all persons; 20 requests of one person hid everybody else's.

Output: at most 3 requests per person per poll, the oldest first; the
other persons' requests are in the same poll (control: a single request is
returned as it is).
"""
from __future__ import annotations

import asyncio
import time

from _messenger_people_rig import Account
from delfin.messenger import pairing


def _poll(acct):
    return asyncio.run(acct.bridge._run_command({"cmd": "pairing", "action": "poll"}))["requests"]


def test_a_flood_of_one_person_leaves_room_for_the_others(tmp_path, monkeypatch):
    acct = Account(tmp_path, monkeypatch)
    first = []
    for i in range(25):
        rid = pairing.submit("phones", {}, directory=acct.pairing_dir("anna"))
        if i < 3:
            first.append(rid)
        time.sleep(0.001)
    ben = pairing.submit("phones", {}, directory=acct.pairing_dir("ben"))
    got = _poll(acct)
    anna = [r["id"] for r in got if r["person"] == "anna"]
    assert anna == first
    assert [r["id"] for r in got if r["person"] == "ben"] == [ben]


def test_a_single_request_is_returned_control(tmp_path, monkeypatch):
    acct = Account(tmp_path, monkeypatch)
    rid = pairing.submit("phones", {}, directory=acct.pairing_dir("ben"))
    assert [r["id"] for r in _poll(acct)] == [rid]
