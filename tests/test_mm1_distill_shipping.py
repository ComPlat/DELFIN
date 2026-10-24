"""MM1 P1: distill_for_shipping — a session lesson becomes a shipped
knowledge entry only through four gates (evidence, privacy, template,
dedupe). Red on the code without the function."""
import pytest

from delfin.agent import knowledge as kn
from delfin.agent import memory_distill as md


GOOD_TEXT = ("Situation: A long test suite run blocked the shell for its "
             "whole duration.\n"
             "Lesson: Run a suite longer than a few minutes with "
             "bash_background and continue with other work.\n"
             "Why: A plain shell call blocks the turn until the command "
             "ends or its timeout cuts it.")
GOOD_EVIDENCE = "test:tests/test_mm1_distill_shipping.py"


def _clean(tmp_path):
    """No shipped entries, no usage file, in a throwaway folder."""
    return tmp_path / "knowledge"


def test_refuses_without_evidence(tmp_path):
    with pytest.raises(ValueError, match="evidence"):
        md.distill_for_shipping(GOOD_TEXT, [], folder=_clean(tmp_path))


def test_refuses_machine_local_evidence(tmp_path):
    # A job/session id names one machine: not followable by a reader
    # elsewhere, so it is not shippable evidence.
    with pytest.raises(ValueError, match="job:"):
        md.distill_for_shipping(GOOD_TEXT, ["job:12345"],
                                folder=_clean(tmp_path))


def test_refuses_entry_with_home_path(tmp_path):
    text = GOOD_TEXT.replace(
        "bash_background",
        "bash_background (the gate lives in /home/ka/someone/delfin)")
    with pytest.raises(ValueError, match="private"):
        md.distill_for_shipping(text, [GOOD_EVIDENCE],
                                folder=_clean(tmp_path))


def test_refuses_entry_with_a_secret(tmp_path):
    text = GOOD_TEXT + "\nexport DELFIN_TOKEN=ghp_abc123def456ghi789jkl"
    with pytest.raises(ValueError, match="private"):
        md.distill_for_shipping(text, [GOOD_EVIDENCE],
                                folder=_clean(tmp_path))


def test_refuses_entry_with_an_email(tmp_path):
    text = GOOD_TEXT.replace(
        "bash_background",
        "bash_background; ask alice@example.com for the log")
    with pytest.raises(ValueError, match="private"):
        md.distill_for_shipping(text, [GOOD_EVIDENCE],
                                folder=_clean(tmp_path))


def test_refuses_entry_that_tells_the_reader_to_change_a_protected_file(
        tmp_path):
    text = GOOD_TEXT.replace(
        "bash_background",
        "edit delfin/agent/sandbox.py to raise the memory limit")
    with pytest.raises(ValueError, match="private|instruction"):
        md.distill_for_shipping(text, [GOOD_EVIDENCE],
                                folder=_clean(tmp_path))


def test_writes_entry_when_all_gates_pass(tmp_path):
    folder = _clean(tmp_path)
    path = md.distill_for_shipping(GOOD_TEXT, [GOOD_EVIDENCE],
                                   author="nacht-s1", folder=folder)
    assert path is not None and path.is_file()
    e = kn.parse_entry(path.read_text(encoding="utf-8"), slug=path.stem)
    assert e.lesson and e.title and e.review_by
    assert e.provenance and e.provenance[0].evidence == GOOD_EVIDENCE
    assert "bash_background" in e.lesson


def test_written_entry_passes_the_ci_lint_checks(tmp_path):
    folder = _clean(tmp_path)
    path = md.distill_for_shipping(GOOD_TEXT, [GOOD_EVIDENCE],
                                   author="nacht-s1", folder=folder)
    e = kn.parse_entry(path.read_text(encoding="utf-8"), slug=path.stem)
    assert kn.template_findings(e, layer=kn.GLOBAL) == []
    from delfin.agent import knowledge_check as kc
    assert kc.check(e.text, public=True, identity={},
                    unpublished_terms=()) == []


def test_duplicate_is_merged_not_added(tmp_path):
    folder = _clean(tmp_path)
    first = md.distill_for_shipping(GOOD_TEXT, [GOOD_EVIDENCE],
                                    author="nacht-s1", folder=folder)
    second = md.distill_for_shipping(
        GOOD_TEXT.replace("few minutes", "a couple of minutes"),
        ["test:tests/other.py::test_x"], author="nacht-s2", folder=folder)
    # Same lesson: no second file, the first entry carries both
    # provenance lines.
    assert second == first
    e = kn.parse_entry(first.read_text(encoding="utf-8"), slug=first.stem)
    authors = {p.author for p in e.provenance}
    assert authors == {"nacht-s1", "nacht-s2"}
    files = list(folder.glob("*.md"))
    assert len(files) == 1


def test_a_contradiction_is_not_silently_merged(tmp_path):
    folder = _clean(tmp_path)
    md.distill_for_shipping(
        GOOD_TEXT.replace(
            "Run a suite longer than a few minutes with "
            "bash_background and continue with other work.",
            "Run long suites synchronously so the result is seen."),
        ["test:tests/one.py::t"], author="nacht-s1", folder=folder)
    with pytest.raises(ValueError, match="contradicts"):
        md.distill_for_shipping(
            GOOD_TEXT, ["test:tests/other.py::t"], author="nacht-s2",
            folder=folder)
