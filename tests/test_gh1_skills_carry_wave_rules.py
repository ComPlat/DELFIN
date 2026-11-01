"""The six git skills carry the wave-26 measured rules.

Every rule below encodes a REAL failure the wave measured (probe files
committed, cross-review merged into the builder branch, a count
conflict picked instead of summed, red CI repaired from the log text
alone, guards nobody ran locally). The skill is where an agent reads
the rule, so each rule must be there in plain words.
"""
import re
from pathlib import Path

_SKILLS = Path(__file__).resolve().parents[1] / "delfin" / "agent" / "pack" / "skills"


def _read(name: str) -> str:
    return (_SKILLS / name).read_text(encoding="utf-8")


def test_pr_skill_names_branch_hygiene():
    text = _read("pr.md")
    assert re.search(r"(?i)probe.{0,40}(never|not).{0,20}committed", text)
    assert re.search(r"(?i)scratch branch", text)
    assert re.search(r"(?i)merging the reviewed branch", text)


def test_pr_skill_says_main_stays_untouched():
    assert re.search(r"(?i)main.{0,30}is\s+left where", _read("pr.md"))


def test_merge_skill_says_counts_are_summed():
    text = _read("merge.md")
    assert re.search(
        r"(?i)(counter|counted|counts)[\s\S]{0,220}(SUMMING|additive)",
        text), text
    assert re.search(r"(?i)never by picking", text)


def test_new_branch_skill_says_cross_review_is_a_scratch_branch():
    text = _read("new-branch.md")
    assert re.search(r"(?i)scratch branch", text)
    assert re.search(r"(?i)never merge the reviewed branch", text)


def test_push_skill_names_the_repo_wide_guard_families():
    text = _read("push.md")
    for family in (r"test_\*guard\*", r"test_every_\*",
                   r"test_delfins_own_\*"):
        assert re.search(family, text), text


def test_pr_until_green_says_reproduce_before_repair():
    text = _read("pr-until-green.md")
    assert re.search(r"(?i)reproduce.{0,60}(locally|FIRST)", text)
    assert re.search(r"(?i)never repair from the\s+log text", text)
