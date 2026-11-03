"""Adversarial tests for the GH1 skills (reviewer nacht-s2, phase 3).

Attack the skill TEXT the way an agent reads it: is the rule stated so
that it survives paraphrase, does it bind to the right verb (SUM vs
pick), does the negative example show what NOT to do, and does any
added line widen a permission or tell the agent to follow text from
another source? All rules under test come from commit 586fc8bc9.
"""
from __future__ import annotations

import re
from pathlib import Path

SKILLS = Path("delfin/agent/pack/skills")

MERGE = (SKILLS / "merge.md").read_text(encoding="utf-8")
PR = (SKILLS / "pr.md").read_text(encoding="utf-8")
NEW_BRANCH = (SKILLS / "new-branch.md").read_text(encoding="utf-8")
PUSH = (SKILLS / "push.md").read_text(encoding="utf-8")
PR_UG = (SKILLS / "pr-until-green.md").read_text(encoding="utf-8")


class TestCountSummingRule:
    def test_merge_names_the_verb_sum_and_bans_picking(self):
        # The rule must name SUM as the action and picking as the
        # forbidden alternative in the SAME block — a rule that only
        # says "additive" without the verb can be read as "pick one".
        step6 = MERGE[MERGE.index("COUNTER"):MERGE.index("7. Run")]
        assert re.search(r"\bSUM(MING)?\b", step6)
        assert re.search(r"never.*pick|pick.*never", step6, re.IGNORECASE)

    def test_merge_requires_sources_in_the_report(self):
        step6 = MERGE[MERGE.index("COUNTER"):MERGE.index("7. Run")]
        assert "sources" in step6.lower()


class TestBranchHygieneRule:
    def test_pr_step0_names_probe_files_and_scratch_branch(self):
        step0 = PR[PR.index("## Steps"):PR.index("1. GitHub CLI")]
        assert "test_zz_" in step0
        assert "scratch branch" in step0
        assert "main" in step0

    def test_new_branch_says_never_merge_the_reviewed_branch(self):
        assert "never merge the reviewed branch" in NEW_BRANCH
        assert "scratch branch" in NEW_BRANCH


class TestGuardFamiliesRule:
    def test_push_names_all_three_families(self):
        # The wave-23 lesson: guards are not all named *guard*.
        step4 = PUSH[PUSH.index("4. Before pushing"):PUSH.index("5. Push")]
        assert "test_*guard*" in step4
        assert "test_every_" in step4
        assert "test_delfins_own_" in step4

    def test_push_keeps_the_original_contract(self):
        # NOTHING NEW IS REFUSED: the original ruff+tests requirement
        # must still be there next to the extension.
        step4 = PUSH[PUSH.index("4. Before pushing"):PUSH.index("5. Push")]
        assert "ruff" in step4
        assert "tests that cover" in step4


class TestRedCIReproductionRule:
    def test_pr_until_green_bans_repair_from_log_alone(self):
        step3 = PR_UG[PR_UG.index("3. Reproduce"):PR_UG.index("4. Decide")]
        assert "locally FIRST" in step3
        assert "log text alone" in step3


class TestSkillSafety:
    """No skill text may widen a permission or outsource a decision."""

    def test_no_push_or_merge_to_main_instruction(self):
        for name, text in (("pr", PR), ("push", PUSH), ("merge", MERGE)):
            # Skills may describe the user's request, never instruct an
            # unprompted push of the base branch.
            assert not re.search(r"push (the )?main\b", text, re.IGNORECASE), \
                name

    def test_no_follow_text_from_outside_instruction(self):
        for name, text in (("pr", PR), ("merge", MERGE), ("push", PUSH),
                           ("pr-until-green", PR_UG), ("new-branch",
                                                       NEW_BRANCH)):
            assert not re.search(
                r"follow (the )?text|instructions in the (log|page|file)",
                text, re.IGNORECASE), name
