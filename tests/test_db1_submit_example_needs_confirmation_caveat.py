"""A state-changing chain in the doc's own example must match its safety rule.

``dashboard_agent.md`` safety rule 1 says: for a state-changing command
(submit, recalc, cancel, delete) describe what you would do, ask, wait
for yes, THEN output the ACTION. The multi-step example in the /done
section shows the opposite shape: ACTION: /control set -> ACTION:
/tab submit -> ACTION: /done with no ask and no wait. A model that
copies the example undermines the rule the doc itself states.

The example is not wrong when the person already asked for the whole
chain — rule 1's "wait for yes" is then satisfied by the request. The
pinned invariant: the example must not be presented as a general
pattern without the confirmation caveat. The section holding it must
say which of its commands are state-changing and that safety rule 1
applies to them.
"""

import re
from pathlib import Path

_DOC = (Path(__file__).resolve().parents[1]
        / "delfin" / "agent" / "pack" / "agents" / "dashboard_agent.md")


def test_the_done_example_is_annotated_for_state_changing_commands():
    body = _DOC.read_text(encoding="utf-8")
    # The multi-step example with /control set -> /tab submit exists.
    assert "ACTION: /control" in body and "ACTION: /tab submit" in body
    m = re.search(
        r"Example — multi-step:.*?(?=\n## |\nIf you forget)", body, re.S)
    assert m, "multi-step example section not found in the doc"
    section = m.group(0)
    assert re.search(
        r"state[- ]changing|safety rule 1|confirmed", section, re.I), (
        "the multi-step example does not say that state-changing commands "
        "in it (submit) are subject to safety rule 1 — a model copying "
        "the example submits without the person's confirmation")


def test_safety_rule_1_names_its_commands():
    """Rule 1 must name the state-changing commands explicitly, so the
    two passages cannot drift apart silently."""
    body = _DOC.read_text(encoding="utf-8")
    rule1 = re.search(
        r"Never run a destructive action without asking\..*?output the `ACTION:`\.",
        body, re.S)
    assert rule1, "safety rule 1 not found"
    for cmd in ("Recalc", "submit", "cancel"):
        assert cmd in rule1.group(0), f"rule 1 does not name {cmd}"
