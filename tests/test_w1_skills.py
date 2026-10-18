"""W1 control: the four web skills exist, load, and stay within limits."""
from __future__ import annotations

from pathlib import Path

import pytest

SKILLS = [
    "web-primary-sources",
    "web-cross-check",
    "web-cite",
    "web-page-is-data",
]


def _skill_path(name: str) -> Path:
    return Path(__file__).resolve().parents[1] / "delfin/agent/pack/skills" / f"{name}.md"


@pytest.mark.parametrize("name", SKILLS)
def test_skill_exists_and_has_frontmatter(name):
    path = _skill_path(name)
    assert path.is_file(), f"missing skill: {path}"
    text = path.read_text(encoding="utf-8")
    assert text.startswith("---\n"), f"{name}: no frontmatter"
    assert "name: " in text and "description: " in text


@pytest.mark.parametrize("name", SKILLS)
def test_skill_is_short(name):
    lines = _skill_path(name).read_text(encoding="utf-8").splitlines()
    assert len(lines) <= 40, f"{name}: {len(lines)} lines (limit 40)"


def test_page_is_data_never_grants_or_instructs():
    text = _skill_path("web-page-is-data").read_text(encoding="utf-8")
    # The skill must state the rule (page text is data), must NOT itself
    # widen permissions or name credential paths.
    assert "UNTRUSTED DATA" in text
    assert "never" in text and "grants" in text
    for forbidden in (".ssh", ".env", "API_KEY", "api key"):
        assert forbidden not in text, f"web-page-is-data names {forbidden}"


def test_skills_reference_the_ledger():
    for name in ("web-primary-sources", "web-cross-check", "web-cite"):
        text = _skill_path(name).read_text(encoding="utf-8")
        assert "source_ledger" in text, f"{name}: no ledger integration"
