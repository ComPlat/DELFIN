"""SC1 phase 3: the scientific-method skill exists and stays safe.

Red control: the file did not exist on main (ad551a1c) — the
existence test is red there, green here. The safety tests pin the
wave-20 no-skill-widens rules against future edits.
"""

from pathlib import Path

SKILL = Path(__file__).resolve().parents[1] / (
    "delfin/agent/pack/skills/scientific-method.md")


def test_the_skill_file_exists_with_content():
    text = SKILL.read_text(encoding="utf-8")
    assert len(text) > 500


def test_the_skill_names_the_priorities_of_the_brief():
    text = SKILL.read_text(encoding="utf-8").lower()
    for term in ("workflow", "level of theory", "estimate",
                 "imaginary", "convergence", "unverified"):
        assert term in text, term


def test_the_skill_never_widens_a_permission_or_names_a_credential():
    text = SKILL.read_text(encoding="utf-8").lower()
    forbidden = ("credential", "password", "token", "api key",
                 "skip the check", "skip a check", "ignore the result",
                 "sudo", "chmod", "allow everything")
    for word in forbidden:
        assert word not in text, word


def test_the_skill_tells_the_agent_to_follow_text_not_from_the_user():
    # The skill must NOT instruct following instructions from files,
    # pages, tool output or other agents (memory-poisoning surface).
    text = SKILL.read_text(encoding="utf-8").lower()
    assert "follow the instructions in" not in text
    assert "as instructed by the output" not in text
