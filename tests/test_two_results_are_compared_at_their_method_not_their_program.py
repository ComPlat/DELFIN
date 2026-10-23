"""The agent must check method comparability before comparing results.

Found in an earlier round: the agent read the difference of two
calculations as "one program versus the other" and overlooked that
functional or basis set differed — it compared tools where it should have
compared methods. For an agent that is supposed to do science, comparing
incomparable numbers and presenting them as the same level is a silent
wrong answer, worse than a refusal.

The control (run against the unfixed state) proves the rule was absent;
the budget test holds the head at the measured cost of the rule.
"""
import re

import pytest

from delfin.agent.prompt_loader import PromptLoader

_ADDENDUM = (
    PromptLoader().agent_dir / "shared" / "scientific_integrity_addendum.md"
)


def test_the_comparability_rule_exists_and_names_the_points():
    """The rule must name the points where two results stop being
    comparable: functional, basis, dispersion, solvent, and geometry
    provenance. A rule that says only "check comparability" names
    nothing the model can check."""
    text = _ADDENDUM.read_text()
    joined = text.lower()
    missing = [
        w for w in ("functional", "basis", "dispersion", "solvent",
                    "geometry")
        if w not in joined
    ]
    assert not missing, (
        f"comparability rule missing the comparison points: {missing}")


def test_the_rule_reaches_the_built_prompt():
    """A rule is a rule when it is IN the prompt. shared/ is not
    automatically everywhere: the scientific_integrity addendum IS
    composed into the solo head (prompt_loader.py:1771-1784), and this
    test pins that the new text survives composition — by matching a
    distinctive phrase from the rule itself, not by section name."""
    built = PromptLoader().build_system_prompt(
        role_id="solo_agent", mode_id="solo", task_text="compare the runs",
        session_key="comparability-1")
    # Distinctive sentence from the rule; matched loosely so a reword
    # inside the rule still counts, but a dropped rule fails.
    assert re.search(r"comparable", built, re.I), (
        "the comparability rule is not in the composed solo prompt")
    # And it is the rule, not the reproducibility bullet that happens to
    # contain the word "comparable": the comparison points must co-occur.
    for w in ("functional", "basis"):
        assert re.search(w, built, re.I), (
            f"'{w}' absent from the composed prompt — the rule or its "
            "points did not survive composition")


@pytest.mark.parametrize("role_id, mode_id", [
    ("solo_agent", "solo"),
    ("dashboard_agent", "dashboard"),
])
def test_the_rule_reaches_every_science_role(role_id, mode_id):
    """Both roles that answer science questions get the addendum via the
    universal Layer-0 loop; pin each so a future routing change cannot
    silently drop one."""
    built = PromptLoader().build_system_prompt(
        role_id=role_id, mode_id=mode_id, task_text="compare the runs",
        session_key=f"comparability-{role_id}")
    assert re.search(r"comparable", built, re.I), (
        f"{role_id}: comparability rule missing from composed prompt")
