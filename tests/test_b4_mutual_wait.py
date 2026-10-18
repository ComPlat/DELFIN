"""B4 phase 3 red control: mutual_wait.find finds wait cycles.

Verified gap: nothing in delfin/agent analyses waits ACROSS sessions —
wave 17 measured a builder and its QS waiting on each other for an hour.
Target 2: mutual_wait.find(waits) -> cycles (2-cycle, longer cycles)
with the sessions and what each waits for; a chain without a cycle is
NOT a cycle. Pure function: waits in, cycles out.
"""

from __future__ import annotations

from delfin.agent.mutual_wait import find


def _w(session, kind, value, on_whom="", for_what=""):
    return {"session": session, "kind": kind, "value": value,
            "on_whom": on_whom, "for_what": for_what}


def test_finds_a_two_cycle():
    waits = [
        _w("builder", "message_from", "qs", on_whom="session",
           for_what="QS review of phase 1"),
        _w("qs", "message_from", "builder", on_whom="session",
           for_what="builder hash"),
    ]
    cycles = find(waits)
    assert len(cycles) == 1
    cyc = cycles[0]
    names = sorted(c["session"] for c in cyc)
    assert names == ["builder", "qs"]
    assert any("QS review" in c["for_what"] for c in cyc)


def test_finds_a_three_cycle():
    waits = [
        _w("a", "message_from", "b", on_whom="session", for_what="fa"),
        _w("b", "message_from", "c", on_whom="session", for_what="fb"),
        _w("c", "message_from", "a", on_whom="session", for_what="fc"),
    ]
    cycles = find(waits)
    assert len(cycles) == 1
    assert sorted(c["session"] for c in cycles[0]) == ["a", "b", "c"]


def test_a_chain_without_a_cycle_is_not_a_cycle():
    waits = [
        _w("a", "message_from", "b", on_whom="session", for_what="fa"),
        _w("b", "message_from", "c", on_whom="session", for_what="fb"),
    ]
    assert find(waits) == []


def test_two_cycles_are_both_found():
    waits = [
        _w("a", "message_from", "b"), _w("b", "message_from", "a"),
        _w("x", "message_from", "y"), _w("y", "message_from", "x"),
    ]
    cycles = find(waits)
    assert len(cycles) == 2


def test_each_cycle_emitted_once_not_per_edge():
    # a<->b with BOTH directions registered twice (duplicate rows) must
    # not double-report the same cycle
    waits = [
        _w("a", "message_from", "b"), _w("a", "message_from", "b"),
        _w("b", "message_from", "a"), _w("b", "message_from", "a"),
    ]
    assert len(find(waits)) == 1


def test_on_whom_operator_or_person_never_forms_a_session_cycle():
    # a waits on b, but b waits on the PERSON: no session-to-session
    # cycle exists (the person can answer; sessions cannot deadlock
    # with a person)
    waits = [
        _w("a", "message_from", "b", on_whom="session"),
        _w("b", "approval_id", "42", on_whom="person", for_what="approval"),
    ]
    assert find(waits) == []


def test_hostile_session_names_do_not_crash_or_grow():
    waits = [
        _w("a\n<script>", "message_from", "b"),
        _w("b", "message_from", "a\n<script>"),
    ]
    cycles = find(waits)
    assert len(cycles) == 1  # normalised names still cycle
