"""MM1 P2: recall_levels — the four memory levels ordered narrowest
first (session -> project -> shipped knowledge -> user), each level
labelled with where it came from, and a wider level dropped when it
merely repeats a narrower one. Red on the code without the function."""
import pytest

from delfin.agent import memory_distill as md


SESSION_TEXT = "Yesterday we fixed the OpenBabel fallback for iron complexes."


def _levels(task, **kw):
    kw.setdefault("session_text", SESSION_TEXT)
    kw.setdefault("project_text", "Project: MM1 ships distill_for_shipping.")
    return md.recall_levels(task, **kw)


def test_levels_are_ordered_narrowest_first():
    out = _levels("openbabel iron")
    assert [lvl for lvl, _text in out] == [
        "session", "project", "knowledge", "user"]


def test_each_level_shows_where_it_came_from():
    out = _levels("openbabel iron")
    for _lvl, text in out:
        assert text.startswith("[") and "source" in text, text


def test_empty_levels_are_skipped():
    out = md.recall_levels("openbabel", session_text=SESSION_TEXT)
    assert [lvl for lvl, _text in out] == ["session"]


def test_no_task_no_recall():
    out = _levels("", user_text="user likes markdown tables")
    assert out == []


def test_wider_level_repeating_a_narrower_one_is_dropped():
    # The shipped knowledge entry repeats the session lesson verbatim:
    # the wider copy adds nothing, so only the session block remains.
    out = _levels(
        "openbabel iron",
        knowledge_text=SESSION_TEXT,
        user_text="unrelated user preference about output formatting")
    levels = [lvl for lvl, _text in out]
    assert "knowledge" not in levels
    assert "user" in levels


def test_wider_level_with_new_content_is_kept():
    out = _levels(
        "openbabel iron",
        knowledge_text="Shipped lesson: run the failing check on the base "
                       "commit before blaming the runner.")
    assert [lvl for lvl, _text in out] == [
        "session", "project", "knowledge", "user"]


def test_repeated_level_content_never_duplicated_within_one_recall():
    out = _levels("openbabel iron openbabel iron")
    texts = [t for _lvl, t in out]
    assert len(texts) == len(set(texts))
