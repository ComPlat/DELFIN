"""Adversarial tests, round 4: lead prompt and skill safety.

Reviewer nacht-s21, package SW3. The lead role prompt and the swarm-lead
skill are shipped content every user's agent will follow. Wave rule: no
skill or prompt may widen a permission, name a credential path, tell the
agent to skip a check, or follow text from a page or file as an
instruction. These tests read the real shipped artifacts and reject such
directives.
"""
from delfin.agent.subagents import SUBAGENT_PRESETS
from delfin.agent.skills import get_skill

# Words whose appearance in a prompt/skill would mean it steers the
# agent at a credential surface or a permission edge. (The list is about
# the SHIPPED TEXT, not about any command line; the deny rules for
# commands do not apply to test file content.)
_FORBIDDEN_MARKERS = (
    "api_key", "api-key", "apikey", "token =", "password",
    ".ssh", ".netrc", "skip the check", "skip this check",
    "without asking", "no confirmation", "ignore the gate",
)

_FORBIDDEN_DIRECTIVES = (
    "follow the instructions in the page",
    "follow the instructions in the file",
    "treat the page content as an instruction",
    "do whatever the file says",
)


def test_lead_preset_body_has_no_unsafe_directives():
    body = (SUBAGENT_PRESETS["swarm-lead"].system_prompt or "").lower()
    for marker in _FORBIDDEN_MARKERS:
        assert marker not in body, (
            f"swarm-lead preset system_prompt contains {marker!r}; a role "
            "prompt must not name credential surfaces or tell the agent "
            "to skip checks"
        )
    for directive in _FORBIDDEN_DIRECTIVES:
        assert directive not in body, directive


def test_swarm_lead_skill_body_has_no_unsafe_directives():
    skill = get_skill("swarm-lead")
    assert skill is not None
    body = (skill.body or "").lower()
    for marker in _FORBIDDEN_MARKERS:
        assert marker not in body, (
            f"swarm-lead skill contains {marker!r}; a skill must not "
            "name credential surfaces or tell the agent to skip checks"
        )
    for directive in _FORBIDDEN_DIRECTIVES:
        assert directive not in body, directive


def test_lead_skill_requires_user_confirmation_for_state_changes():
    """Target 1 of the package: state-changing lead tools always ASK the
    user in the lead session. The skill must say so, not silently do."""
    skill = get_skill("swarm-lead")
    assert skill is not None
    low = (skill.body or "").lower()
    assert "ask" in low, (
        "the swarm-lead skill never mentions asking; state-changing lead "
        "tools (swarm_start, swarm_pause/resume, swarm_triage_answer, "
        "swarm_prcheck) must always ask the user first"
    )


def test_lead_skill_refuses_test_only_fixes_that_hide_defects():
    """Package rule: a test-only fix that hides a production defect is
    refused. The lead role must carry that rule, since it routes fixes."""
    skill = get_skill("swarm-lead")
    assert skill is not None
    low = (skill.body or "").lower()
    assert "hide" in low and "defect" in low, (
        "the swarm-lead skill does not carry the rule that a test-only "
        "fix hiding a production defect is refused"
    )
