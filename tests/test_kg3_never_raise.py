"""Probe for nacht-s22's two IndexError triggers (never-raise contract)."""
from delfin.agent import prompt_head as ph


def test_none_and_odd_turns_never_raise():
    a = [{"name": f"r{i}", "layer": 0, "content": "x"} for i in range(5)]
    # trigger (a): turn_a = 5 rows, turn_b = None
    assert ph.stable_prefix_len(a, None)["rows"] == 0
    # trigger (b): turn_a = None, turn_b = 1 row
    b = [{"name": "only", "layer": 0, "content": "y"}]
    r = ph.stable_prefix_len(None, b)
    assert r["rows"] == 0
    assert r["first_change_row"] == 0
    assert r["first_change_section"] == "only"
    # hostile rows: non-dicts, missing keys
    h1 = [None, "text", {"name": "n"}]
    h2 = [None, {"layer": "two", "content": None}]
    ph.stable_prefix_len(h1, h2)
    ph.stable_prefix_len(h1, None)
    ph.check_ordering(h1 + h2)
    ph.hostile_reach(h1 + h2)
    ph.newest_person_instruction(h1 + h2)
