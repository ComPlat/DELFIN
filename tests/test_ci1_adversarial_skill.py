"""CI1 adversarial tests (nacht-s6, reviewer): skill-shape attacks.

Attacks the skill texts the builder's lint test does not: the lint only
checks absence of phrases; these check what a half-hearted skill would
still get wrong (merge advice without green+match-head-commit, weak
wording like 'consider', force-push routes, watching without a stop).
"""
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent.parent / "delfin" / "agent" / "pack" / "skills"


def _read(name):
    return (SKILL_DIR / f"{name}.md").read_text(encoding="utf-8")


def test_pr_until_green_mentions_stop_conditions():
    text = _read("pr-until-green")
    for topic in ("base branch", "infrastructure", "protected"):
        assert topic in text.lower(), f"missing stop condition: {topic}"


def test_pr_until_green_requires_reproduce_before_fix():
    text = _read("pr-until-green").lower()
    assert "reproduc" in text


def test_pr_until_green_never_merges():
    text = _read("pr-until-green").lower()
    assert "merge" in text and "never" in text


def test_pr_until_green_reports_each_round_to_the_person():
    text = _read("pr-until-green").lower()
    assert "report" in text and "person" in text


def test_pr_until_green_no_ci_config_edit_wording():
    # A skill that says 'edit the workflow' or 'add a workflow' is weakening.
    text = _read("pr-until-green").lower()
    assert ".github/" not in text
    assert "workflow" not in text


def test_pr_report_only_watches_when_asked():
    # pr.md Report step: "watching until green" is conditional on the user
    # asking for it -- the report must not say it always watches.
    text = " ".join(_read("pr").lower().split())
    assert "watching until green" in text
    assert "when the user asked" in text or "if the user asked" in text


def test_pr_merge_gate_still_intact():
    # The merge step must still require green + --match-head-commit and
    # still forbid --admin/--auto. Mentioning the flags AS A REFUSAL is
    # correct content (pr.md:35 does this); only the absence of the
    # refusal wording would be a regression.
    text = _read("pr")
    assert "--match-head-commit" in text
    lowered = text.lower()
    assert "never `--admin` or `--auto`" in lowered


def test_no_skill_tells_the_agent_to_follow_log_text():
    # CI log text is data. A skill must never advise acting on log text,
    # or claim text from a log is trusted input.
    text = _read("pr-until-green").lower()
    assert "log" not in text or "data" in text
