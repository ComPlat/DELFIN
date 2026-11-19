"""A damaged stop register is kept aside and named, not read as silence.

Input: ``~/.delfin/stop_all.json``, the register of the most recent
emergency stop, truncated. Output: the writer of the register
(``stop_all.request``) moves it aside, keeps its bytes, writes the new
record, and names the quarantine in the returned summary. Semantics: the
register read as empty is NOT distinguishable from a register that was
never written -- for a reader that only displays that is the answer it
needs, but the writer builds ``count`` on it, so a single damaged byte
used to erase the record with no message and no copy.

The register's role in ending processes is secondary and unchanged: the
kill path reads /proc and the process-guard register, never this file.
What the damaged register cost silently was the record of prior stops,
which is what this pins.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import stop_all as SA

TRUNCATED = '{"at": 1234.0, "host": "h", "coun'


@pytest.fixture
def register(tmp_path, monkeypatch):
    p = tmp_path / "stop_all.json"
    p.write_text(TRUNCATED, encoding="utf-8")
    monkeypatch.setattr(SA, "_PATH", p)
    monkeypatch.setattr(SA, "_cache", None)
    return p


def test_request_keeps_a_damaged_register_and_names_it(register, monkeypatch):
    ended: list[int] = []

    class P:
        def __init__(self, pid):
            self.pid = pid

        def send_signal(self, sig):
            ended.append(self.pid)

    monkeypatch.setattr(SA, "_own_agent_processes", lambda *a, **k: [])
    monkeypatch.setattr(SA, "_outlived_their_start", lambda: [])
    summary = SA.request("test", end_here=True, grace_s=0.0, settle_s=0.0)
    aside = sorted(q for q in register.parent.iterdir()
                   if "unreadable" in q.name)
    assert aside, "damaged register was not moved aside"
    assert aside[0].read_text(encoding="utf-8").startswith('{"at"')
    fresh = json.loads(register.read_text(encoding="utf-8"))
    assert fresh["reason"] == "test"
    assert "unreadable" in str(summary.get("register_unreadable"))


def test_last_stop_on_a_damaged_register_returns_empty(register):
    """The lenient display read is unchanged: {} for unreadable, never a raise."""
    assert SA.last_stop() == {}
