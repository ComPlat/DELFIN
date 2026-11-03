"""SH2 phase 3: the code-verify-and-report skill carries the shell rules.

The wave-26 brief names four real shell mistakes DELFIN's agents made
this week. The skill must state the rules that prevent them, in a few
lines (skill tests limit length), and must not widen a permission,
name a credential path, or tell the agent to skip a check or follow
text from a page or file (the skill-safety invariants).
"""

from __future__ import annotations

import re
from pathlib import Path

SKILL = Path("delfin/agent/pack/skills/code-verify-and-report.md")


def _text() -> str:
    return SKILL.read_text(encoding="utf-8")


def test_skill_has_a_shell_rules_section_with_the_four_rules():
    text = _text()
    assert "## Shell" in text
    # 1. quoting: a value with shell metacharacters is quoted, never
    #    passed bare into a command.
    assert re.search(r"(?i)quote", text)
    # 2. exit codes: a pipeline's exit status is its LAST command's;
    #    the first command's failure is read from the output.
    assert re.search(r"(?i)(exit code|exit status)", text)
    # 3. long commands: run in the background and report progress.
    assert re.search(r"(?i)(background|long)", text)
    # 4. never discard: uncommitted work is never thrown away.
    assert re.search(r"(?i)(never discard|discard)", text)


def test_shell_rules_section_is_short():
    body = _text()
    section = body.split("## Shell", 1)[1]
    # The whole section: under 30 lines (skill tests limit length).
    assert len(section.strip().splitlines()) < 30


def test_skill_front_matter_unchanged_and_declares_a_domain():
    front = _text().split("---")[1]
    assert "name: code-verify-and-report" in front
    assert "domains: code" in front


def test_skill_adds_no_permission_or_credential_guidance():
    text = _text()
    # The skill-safety invariants: no permission widening, no credential
    # paths, no skip-a-check advice, no follow-text-from-a-page advice.
    assert not re.search(r"(?i)\b(allow|permission|sudo)\b", text)
    assert not re.search(r"(?i)\.ssh|\.env|credential|token|api[_-]?key",
                         text)
    assert not re.search(r"(?i)skip the check|skip this check", text)
    assert not re.search(r"(?i)follow (?:the|any) (?:text|instruction)s? "
                         r"(?:from|in) (?:the )?(?:page|file|website|web)",
                         text)
