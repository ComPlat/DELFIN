"""PC1 phase 2: the hand-over-checked-work skill (nacht-s38).

The skill is the process deliverable: run pre_pr_check, fix every blocker,
only then hand over — and never weaken a test to pass it.  Red first:
on the commit before the skill exists every case here fails (the skill is
not discovered).
"""

from delfin.agent import skills as S


def _the_skill():
    found = [s for s in S.discover_skills() if s.name == "hand-over-checked-work"]
    assert found, "skill hand-over-checked-work.md is not packaged"
    return found[0]


def test_the_skill_is_packaged_and_declares_a_domain():
    skill = _the_skill()
    assert skill.domains, "a packaged skill without domains fails CI (#383)"


def test_the_skill_says_run_the_check_before_handover():
    body = _the_skill().body
    assert "pre_pr_check" in body


def test_the_skill_says_fix_blockers_not_weaken_tests():
    body = _the_skill().body
    assert "BLOCKED" in body
    assert "never weaken" in body


def test_the_skill_names_its_callers():
    body = _the_skill().body
    assert "cli_prcheck" in body or "delfin prcheck" in body
