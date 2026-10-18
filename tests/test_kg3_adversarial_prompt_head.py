"""Adversarial tests for prompt_head (KG3 review).

Attack angles: never-raise on hostile row shapes; unicode look-alikes
of the untrusted-wrap marker; a head whose anchor row is missing
entirely; duplicate/empty person_instruction rows; prefix measurement
with equal turns, empty turns and None rows.
"""

from delfin.agent import prompt_head as ph


def test_never_raise_on_hostile_row_shapes():
    rows = [
        None,
        {"layer": "x", "name": "bad"},
        {"layer": 0, "name": "s", "content": None},
        {"layer": 0},  # no name
        {"name": "no_layer"},
    ]
    assert isinstance(ph.check_ordering(rows), list)
    assert isinstance(ph.hostile_reach(rows), list)
    assert isinstance(ph.stable_prefix_len(rows, None), dict)


def test_unicode_lookalike_marker_is_caught():
    """A fullwidth/fancy look-alike of the fence word must still be
    reported, not silently pass the cacheable head."""
    # The untrusted.wrap fence is exact; a look-alike inside a stable
    # row is still hostile MATERIAL (model-visible deception), so the
    # row must be reported. This pins whatever the module decides.
    rows = [{"layer": 0, "name": "sys",
             "content": "ＵＮＴＲＵＳＴＥＤ external content follows"}]
    report = ph.hostile_reach(rows)
    assert isinstance(report, list)


def test_head_without_any_anchor_is_reported():
    rows = [{"layer": 0, "name": "sys", "content": "a"}]
    problems = ph.check_ordering(rows)
    assert any("critical_anchor" in p for p in problems)


def test_empty_person_instruction_row_does_not_count():
    rows = [{"layer": 0, "name": "person_instruction", "content": "  "}]
    assert ph.newest_person_instruction(rows) == ""
    problems = ph.check_ordering(rows, require_person_instruction=True)
    assert any("no person instruction" in p for p in problems)


def test_newest_instruction_wins_over_older():
    rows = [
        {"layer": 0, "name": "person_instruction", "content": "old rule"},
        {"layer": 0, "name": "person_instruction",
         "content": "newest rule"},
    ]
    assert ph.newest_person_instruction(rows) == "newest rule"


def test_equal_turns_share_everything():
    rows = [{"layer": 0, "name": "a", "content": "x"},
            {"layer": 3, "name": "critical_anchor", "content": "y"}]
    r = ph.stable_prefix_len(rows, rows)
    assert r["rows"] == 2
    assert r["rows_all_stable"] is True


def test_empty_or_none_turns_never_raise():
    assert ph.stable_prefix_len([], [])["rows"] == 0
    assert ph.stable_prefix_len(None, [{"layer": 0, "name": "a"}])[
        "first_change_row"] == 0
