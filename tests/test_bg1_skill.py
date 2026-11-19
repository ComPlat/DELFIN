"""BG1 phase 3: the calc-budget skill is shipped and safe.

The skill is a DELFIN deliverable (package BG1 target 5).  It must exist as a
packaged skill, declare the "code" domain like the rest of the pack, and carry
no language that widens a permission, names a credential path, or tells the
agent to skip a check (wave-28 skill-safety rule).
"""

import re
from pathlib import Path

PACK = Path(__file__).resolve().parents[1] / "delfin" / "agent" / "pack" / "skills"
SKILL = PACK / "calc-budget.md"


def test_skill_exists_with_code_domain():
    assert SKILL.is_file(), "calc-budget.md must be a packaged skill"
    text = SKILL.read_text(encoding="utf-8")
    assert re.search(r"^domains:\s*code\b", text, re.M), \
        "skill must declare the code domain like the rest of the pack"
    assert re.search(r"^name:\s*calc-budget\b", text, re.M), \
        "skill must declare its own name"


def test_skill_text_is_safe_for_the_agent():
    # No permission-widening, credential-path or skip-a-check language -- the
    # same invariant test_gh1_adversarial_skills enforces over the whole pack.
    text = SKILL.read_text(encoding="utf-8")
    forbidden = [r"\bpermissions?\b", r"\bcredentials?\b", r"skip\s+che", r"~/\."]
    for pattern in forbidden:
        assert not re.search(pattern, text, re.I), \
            f"skill text must not contain {pattern!r}"
