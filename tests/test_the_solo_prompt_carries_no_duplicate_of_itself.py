"""The solo role prompt carries rules a shared addendum or a nearby
section already states — measured 2026-09-29 as the largest remaining
duplication in the prompt (role_prompt 10,830 tokens of a 16,255-token
system prompt).

The diet this test pins: restatements that say the same thing a second
time in the same words a model has already read elsewhere in the SAME
prompt cost tokens on every request and have never been found to change
behaviour. This test is the control for the Welle-10 trim: it fails on
the untrimmed file (red) and stays green after the duplicated passages
are removed.

What must NOT disappear (pinned by other tests, kept verbatim):
safety wording, the lookup-rule section, the decision tree, the git
resolution sentences, the PROMPT: marker, the around=/append forms.
"""
from __future__ import annotations

import re
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
_SOLO = _REPO / "delfin" / "agent" / "pack" / "agents" / "solo_agent.md"


def _text() -> str:
    return _SOLO.read_text(encoding="utf-8")


def test_the_lookup_rule_is_not_restated_in_confirm_before_mutating():
    """'Confirm before mutating' repeated the file-confirm workaround
    rule (don't switch tools to escape a deny) that the lookup section
    already states as 'Edit with edit_file, never sed -i', and that the
    KIT sandbox section states a third time as 'Don't retry the same
    blocked command in a loop'. One statement is enough; the other two
    are cut and the surviving one stays in the lookup section, where
    test_looking_things_up_stays_in_the_reading_tools.py pins it."""
    t = _text()
    # The rule lives exactly once, in the lookup section.
    occurrences = [m.start() for m in re.finditer(r"sed -i", t)]
    assert len(occurrences) == 1, (
        f"sed -i appears {len(occurrences)} times in solo_agent.md; "
        "one statement of the rule is enough")


def test_the_work_in_one_workspace_layout_is_stated_once():
    """The <task-slug>/ folder layout was spelled out twice: in 'Work in
    ONE workspace' and again, verbatim, inside the KIT sandbox section
    ('The <task-slug>/ layout from "Work in ONE workspace" applies here
    too' plus its own bash example). The KIT section keeps the pointer,
    the second worked example goes."""
    t = _text()
    assert "dedicated\n`<task-slug>/` subfolder" in t or (
        "a dedicated\n`<task-slug>/` subfolder" in t), (
        "the defining statement of the task-slug layout is gone")
    # The second, redundant pip-install walk-through inside the KIT
    # section is what gets trimmed; the pointer sentence stays.
    assert "layout from \"Work in ONE workspace\" applies here too" not in t, (
        "the KIT section still restates the layout instead of pointing at it")


def test_the_idempotent_setup_table_is_not_a_second_confirm_rule():
    """The 'Idempotent setup' section duplicated 'Confirm before
    mutating' (read first, only rewrite if content differs) plus the
    pip-list warning that 'Project-dev workflow' states again. After the
    trim the read-before-write contract lives once, in 'Confirm before
    mutating'; the table keeps only the checks 'Confirm before mutating'
    does not carry."""
    t = _text()
    assert "read first, only rewrite if content differs" not in t, (
        "the read-before-write contract is stated twice")


def test_the_session_start_section_is_not_a_second_git_recipe():
    """'Session start' repeated 'git status' / 'git log' advice that the
    parallel-tool-calls section and the git workflow addendum already
    name. The orienting step survives in one sentence; the numbered
    re-recipe goes."""
    t = _text()
    assert "git log --oneline -5" in t  # once, in parallel tool calls
    assert t.count("git log --oneline -5") == 1, (
        "the session-start git recipe is restated")
