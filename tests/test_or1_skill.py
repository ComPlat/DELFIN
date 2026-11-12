"""OR1 phase 2: the TS skill exists, is loadable and stays clean.

Red tests written BEFORE the skill file lands the workflow content.
"""
from __future__ import annotations

import re
from pathlib import Path

SKILL = Path("delfin/agent/pack/skills/orca-ts-workflow.md")

# Same SHAPE rules as tests/test_c1_skill_text_lint.py:23-35 — a skill
# must never widen a permission, name credentials, or leave the workspace.
FORBIDDEN = [
    r"(?i)\b(grant|widen|raise|bypass|skip)\b[^\n]{0,30}\b(permission|check|gate)\b",
    r"(?i)\ballow[_ -]?all\b",
    r"(?i)\b(sudo|chmod|chown)\b",
    r"(?i)\b(api[_ -]?key|token|password|secret)\b",
    r"(?i)(\.ssh|\.env|\.netrc|\.npmrc)\b",
    r"(?i)(etc/(passwd|sudoers)|/usr/(lib|share|bin)/\S+)",
    r"(?i)(^|\s)/home/\S+",
    r"(?i)\btmp/\S+",
]


def test_the_ts_skill_exists_and_names_the_delfin_steps():
    text = SKILL.read_text(encoding="utf-8")
    # the workflow rides DELFIN's own steps, nothing invented
    for needle in ("orca_ts_opt", "orca_ts_freq_check", "orca_ts_irc",
                   "xtb_opt", "hessian_mode"):
        assert needle in text, needle
    # the scientific invariants the package exists for
    for needle in ("exactly one imaginary", "IRC", "same level", "ZPE"):
        assert needle in text, needle


def test_the_ts_skill_carries_no_forbidden_shape():
    text = SKILL.read_text(encoding="utf-8")
    hits = [pat for pat in FORBIDDEN if re.search(pat, text)]
    assert hits == [], hits


def test_the_ts_skill_frontmatter_and_intent_line():
    text = SKILL.read_text(encoding="utf-8")
    assert text.startswith("---\n")
    assert "domains:" in text.split("---")[1]
