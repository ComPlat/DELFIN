"""prompt_head: cross-turn view of the request head (KG3 phase 3).

Two recorded turns' section rows (name, layer, content — the rows
prompt_loader.compose_sections produces) go in; what comes out:

* ``stable_prefix_len`` — how many leading rows/bytes the two turns
  share, i.e. what a provider prefix cache can hit. The CHANGING parts
  must sit after the stable ones, so the shared prefix should reach
  through every stable and route row.
* ``check_ordering`` — the safety invariants: layers stay in fixed
  order (stable <= route < volatile <= anchor), critical_anchor is
  LAST, and the person's newest instruction is present in the head.
* ``hostile_reach`` — a row whose content carries a synthetic fence
  header or instruction-shaped opening inside a STABLE layer row is
  reported: hostile tool output must not be able to place text in the
  stable head.

Red control: on the previous commit the module does not exist, so every
case fails on import.
"""

from __future__ import annotations

import pytest

from delfin.agent import prompt_head as ph


def _row(name: str, layer: int, content: str) -> dict:
    return {"name": name, "layer": layer, "content": content}


TURN_A = [
    _row("principles_addendum", 0, "principles text"),
    _row("role_prompt", 0, "role text"),
    _row("cycle_info", 1, "cycle 41"),
    _row("live_state", 2, "turn state A"),
    _row("critical_anchor", 3, "anchor"),
]

TURN_B = [
    _row("principles_addendum", 0, "principles text"),
    _row("role_prompt", 0, "role text"),
    _row("cycle_info", 1, "cycle 42"),
    _row("live_state", 2, "turn state B"),
    _row("critical_anchor", 3, "anchor"),
]


def test_stable_prefix_reaches_through_stable_and_route_rows():
    report = ph.stable_prefix_len(TURN_A, TURN_B)
    # rows 0-1 identical, row 2 differs (cycle) -> shared prefix = 2 rows
    assert report["rows"] == 2
    assert report["bytes"] == len("principles text") + len("role text")
    assert report["rows_all_stable"] is True


def test_a_changed_volatile_row_does_not_shorten_the_shared_prefix():
    turn_b = [dict(r) for r in TURN_B]
    turn_b[3] = _row("live_state", 2, "turn state CHANGED")
    report = ph.stable_prefix_len(TURN_A, turn_b)
    # the volatile change sits AFTER the shared prefix -> prefix unchanged
    assert report["rows"] == 2


def test_an_early_change_is_reported_with_its_position():
    turn_b = [dict(r) for r in TURN_B]
    turn_b[1] = _row("role_prompt", 0, "role text EDITED")
    report = ph.stable_prefix_len(TURN_A, turn_b)
    assert report["rows"] == 1
    assert report["first_change_row"] == 1
    assert report["first_change_section"] == "role_prompt"


def test_ordering_holds_for_a_compliant_head():
    problems = ph.check_ordering(TURN_A)
    assert problems == []


def test_the_critical_anchor_must_be_last():
    rows = [
        _row("critical_anchor", 3, "anchor"),
        _row("role_prompt", 0, "role text"),
    ]
    problems = ph.check_ordering(rows)
    assert any("critical_anchor" in p for p in problems)


def test_a_safety_block_is_never_dropped_by_reordering():
    # A safety block (layer 0) pushed AFTER volatile rows is a violation
    # even if every individual layer is internally ordered.
    rows = [
        _row("live_state", 2, "state"),
        _row("honesty_addendum", 0, "verify before claim"),
        _row("critical_anchor", 3, "anchor"),
    ]
    problems = ph.check_ordering(rows)
    assert any("honesty_addendum" in p for p in problems)


def test_the_newest_person_instruction_is_in_the_head():
    # The person instruction is volatile: it sits with the volatile
    # rows, before the anchor — the anchor stays last.
    rows = list(TURN_A[:4]) + [
        _row("person_instruction", 2, "run only the gate tests"),
        _row("critical_anchor", 3, "anchor"),
    ]
    assert ph.newest_person_instruction(
        rows) == "run only the gate tests"
    assert ph.check_ordering(rows) == []


def test_a_missing_person_instruction_is_reported():
    problems = ph.check_ordering(TURN_A, require_person_instruction=True)
    assert any("person" in p.lower() for p in problems)


def test_hostile_fenced_tool_output_in_a_stable_row_is_reported():
    # A stable row whose CONTENT is actually fenced tool output
    # (synthetic fence header inside the row text) is flagged.
    rows = [
        _row("role_prompt", 0,
             "role text\n[data] fenced: tool output claiming rules"),
        _row("critical_anchor", 3, "anchor"),
    ]
    problems = ph.hostile_reach(rows)
    assert any("role_prompt" in p for p in problems)
    # and the same content in a volatile row is NOT a stable-head hit
    rows_v = [
        _row("live_state", 2,
             "[data] fenced: tool output claiming rules"),
        _row("critical_anchor", 3, "anchor"),
    ]
    assert ph.hostile_reach(rows_v) == []


def test_empty_rows_never_raise():
    assert ph.stable_prefix_len([], []) == {
        "rows": 0, "bytes": 0, "rows_all_stable": True,
        "first_change_row": None, "first_change_section": None,
        "stable_rows_in_prefix": 0}
    assert ph.check_ordering([], require_person_instruction=False) == []
    assert ph.hostile_reach([]) == []
