"""Four operator interviews on two models named the same contradiction:
the critical anchor says "Grep before Read", the chemistry module says
"typed tool BEFORE list/grep" for ORCA output. Both were right and the
model had to pick. The anchor now says both in one breath, in the same
token budget give or take four."""

from delfin.agent.prompt_loader import _CRITICAL_RULES


def test_the_write_rule_names_the_typed_tools_and_keeps_grep_before_read():
    rule = next(r for r in _CRITICAL_RULES["write"] if r.startswith("Grep before Read"))
    assert "extract_*" in rule and "ORCA output" in rule
    assert "relevant lines" in rule
    assert len(rule) <= 95, "the anchor is repeated at the end of every prompt; keep it short"
