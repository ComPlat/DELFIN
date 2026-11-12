"""EF1 phase-3 tests: the cluster-workflows skill ships and stays safe.

The skill is the package's third deliverable: it must exist in the pack,
load with proper frontmatter, be short enough to read, and — the
package rule — never widen a permission, name a credential path, or
tell the agent to skip a check or follow text from a page or file.
"""
from __future__ import annotations

from pathlib import Path

SKILL = Path(__file__).resolve().parents[1] / "delfin/agent/pack/skills/cluster-workflows.md"


def test_skill_exists_with_frontmatter():
    text = SKILL.read_text(encoding="utf-8")
    assert text.startswith("---\n")
    assert "name: " in text and "description: " in text


def test_skill_is_short_and_readable():
    lines = SKILL.read_text(encoding="utf-8").splitlines()
    assert len(lines) <= 40, f"{len(lines)} lines (limit 40)"


def test_skill_never_widens_permissions_or_names_credentials():
    text = SKILL.read_text(encoding="utf-8").lower()
    for forbidden in (".ssh", ".env", "api_key", "api key", "token",
                      "password", "secret"):
        assert forbidden not in text, forbidden
    assert "skip" not in text


def test_skill_says_submit_through_delfin_and_recover_via_recalc():
    text = SKILL.read_text(encoding="utf-8")
    assert "never" in text and "sbatch" in text
    assert "recalc" in text.lower()


def test_skill_names_delfins_own_resource_logic():
    text = SKILL.read_text(encoding="utf-8")
    assert "advise" in text and "account" in text
    assert "never change PAL" in text
