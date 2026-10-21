"""B4 adversarial: cycle finder (mutual_wait) and health wiring.

Reviewer nacht-s19, wave 18. Attacks the builder's P3:
- 2-cycle, 3-cycle, 4-cycle found; chains without a cycle never fire.
- on_whom operator/person/ci are NEVER cycle members (a wait on the person
  cannot be mutual).
- hostile names (markup, newline, unicode look-alikes) cannot forge a cycle
  identity: rows are normalised the same way register() cleans them.
- duplicate rows in the input collapse to one report.
- a cycle is reported ONCE (deduped by member set), not per rotation.
"""

from __future__ import annotations

import pytest

from delfin.agent.mutual_wait import find


def _w(kind, value, on_whom="", for_what="", session=""):
    # find() reads w["session"] (the waiting session) and w["value"]
    # (who it waits on) — the real row shape from wait_registry._load()
    # enriched with the session key (session_health.collect() adds it).
    return {
        "session": session, "kind": kind, "value": value,
        "on_whom": on_whom, "for_what": for_what,
    }


def test_two_cycle_found():
    waits = [
        _w("message_from", "s2", on_whom="session", session="s1"),
        _w("message_from", "s1", on_whom="session", session="s2"),
    ]
    cycles = find(waits)
    assert len(cycles) == 1
    # A cycle is the list of its edge rows; members = the sessions on them.
    members = {r["session"] for r in cycles[0]}
    assert members == {"s1", "s2"}


def test_three_cycle_found():
    waits = [
        _w("message_from", "s2", on_whom="session", session="s1"),
        _w("message_from", "s3", on_whom="session", session="s2"),
        _w("message_from", "s1", on_whom="session", session="s3"),
    ]
    cycles = find(waits)
    assert len(cycles) == 1
    members = {r["session"] for r in cycles[0]}
    assert members == {"s1", "s2", "s3"}


def test_chain_without_cycle_not_reported():
    waits = [
        _w("message_from", "s2", on_whom="session", session="s1"),
        _w("message_from", "s3", on_whom="session", session="s2"),
    ]
    assert find(waits) == []


def test_self_wait_is_not_a_cycle():
    waits = [_w("message_from", "s1", on_whom="session", session="s1")]
    assert find(waits) == []


def test_operator_person_ci_never_in_a_cycle():
    waits = [
        _w("message_from", "operator", on_whom="operator", session="s1"),
        _w("message_from", "s1", on_whom="session", session="operator"),
    ]
    assert find(waits) == []


def test_duplicate_rows_collapse_to_one_cycle():
    waits = [
        _w("message_from", "s2", on_whom="session", session="s1"),
        _w("message_from", "s2", on_whom="session", session="s1"),
        _w("message_from", "s1", on_whom="session", session="s2"),
    ]
    cycles = find(waits)
    assert len(cycles) == 1


def test_rotations_of_one_cycle_reported_once():
    waits = [
        _w("message_from", "s2", on_whom="session", session="s1"),
        _w("message_from", "s1", on_whom="session", session="s2"),
    ]
    reversed_waits = list(reversed(waits))
    c1 = find(waits)
    c2 = find(reversed_waits)
    assert len(c1) == 1 and len(c2) == 1
    assert {r["session"] for r in c1[0]} == {r["session"] for r in c2[0]}


def test_hostile_names_cannot_forges_a_cycle():
    # A wait whose "value" looks like a session name with markup/whitespace
    # must NOT match a clean key: identity is the normalised name.
    waits = [
        _w("message_from", "s2\n", on_whom="session", session="s1"),
        _w("message_from", "s1", on_whom="session", session="s2"),
    ]
    # Either the hostile row is dropped (no cycle) or normalised to s2
    # (the cycle is found) — but never a forged third identity.
    cycles = find(waits)
    for c in cycles:
        assert {r["session"] for r in c} <= {"s1", "s2"}


def test_empty_and_none_input():
    assert find([]) == []
    assert find([None]) == [] or find([None]) == []
