"""A refusal is remembered with the user's reason, and survives compaction.

A session already could not ask twice for the same refused path or
command -- but the second answer said only "already declined", without
the reason the user gave the first time, and a compaction dropped the
refusals the summary did not happen to quote. The permissions object now
keeps them (`refusal_memory`, in memory: a file in the workspace could be
edited by the agent it describes), the repeat answer carries the reason,
and the engine hands them to the working-state block.
"""

from __future__ import annotations

from types import SimpleNamespace

from delfin.agent import api_client
from delfin.agent.engine import AgentEngine


def test_a_remembered_refusal_answers_a_repeat_with_its_reason():
    perms = SimpleNamespace(refusal_memory=None)
    api_client._remember_refusal(perms, "read_file", "/srv/x/gate",
                                 "outside your workspace; use the gate")
    note = api_client._earlier_refusal_reason(
        perms, "read_file", {"path": "/srv/x/gate"})
    assert "outside your workspace; use the gate" in note
    # The same target through the shell is the same ask.
    assert api_client._earlier_refusal_reason(
        perms, "bash", {"command": "sed -n 1,50p /srv/x/gate"})
    # Another file in the same directory was never refused.
    assert api_client._earlier_refusal_reason(
        perms, "read_file", {"path": "/srv/x/lint"}) == ""


def test_no_memory_and_a_broken_memory_say_nothing():
    assert api_client._earlier_refusal_reason(
        SimpleNamespace(), "bash", {"command": "cat /x"}) == ""
    broken = SimpleNamespace(refusal_memory=object())
    assert api_client._earlier_refusal_reason(
        broken, "bash", {"command": "cat /x"}) == ""
    api_client._remember_refusal(object(), "bash", "cat /x", "r")  # no raise


def test_the_engine_hands_the_refusals_to_the_working_state():
    perms = SimpleNamespace(refusal_memory=None)
    api_client._remember_refusal(perms, "bash", "cat /tmp/spy.log",
                                 "logs belong in your own tree")
    engine = SimpleNamespace(client=SimpleNamespace(_permissions=perms))
    engine.kit_permissions = perms
    entries = AgentEngine._refusal_entries(engine)
    assert entries and entries[0]["reason"] == "logs belong in your own tree"
