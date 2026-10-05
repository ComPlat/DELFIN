"""A session shows a name, not the first 40 characters of its prompt.

A session's ``title`` is the text the user opened with. Measured over
the 52 sessions in ``~/.delfin/agent_sessions`` on 2026-10-05: **38 are
unusable as a label** -- multi-line, markdown headings, code fences, cut
off mid-word. Two places show them as names:

* the Resume dropdown, where each row reads
  ``"Du bist heute nicht nur Assistent, sondern auch Tester dieses W — ~/x"``;
* the new-session hint that warns another session is already working in
  this repository, which joins up to three of them into one sentence.

So the control that exists to stop two sessions sharing a checkout
([[two worktrees]]) was unreadable exactly when it mattered.

Truncation is not the fix: the first 40 characters of a handover
briefing are "# Übergabe: Sechs Arbeitszweige prüfen, z". The first LINE,
stripped of markup and cut on a word boundary, is a name.
"""

from __future__ import annotations

from delfin.dashboard import agent_sessions as AS


def test_a_plain_short_title_is_left_alone():
    assert AS.session_label("Fix the failing test", "abc123") == (
        "Fix the failing test")


def test_only_the_first_line_is_used():
    assert AS.session_label(
        "Rename the branch\n\nThen push it and open a PR", "x") == (
        "Rename the branch")


def test_a_markdown_heading_loses_its_marks():
    assert AS.session_label("# Übergabe: Sechs Arbeitszweige", "x") == (
        "Übergabe: Sechs Arbeitszweige")
    assert AS.session_label("## M3 — Stumme Hänger sichtbar machen", "x") == (
        "M3 — Stumme Hänger sichtbar machen")


def test_fences_and_emphasis_are_dropped():
    assert "`" not in AS.session_label("Fix `api_client.py` now", "x")
    assert "*" not in AS.session_label("**Wichtig**: push the branch", "x")
    # A title that STARTS with a fence has its first real line taken.
    assert AS.session_label("```\ngit push\n```", "x") == "git push"


def test_a_long_title_is_cut_on_a_word_boundary():
    long = ("Du bist heute nicht nur Assistent, sondern auch Tester "
            "dieses Werkzeugkastens.")
    got = AS.session_label(long, "x")
    assert len(got) <= 44, got
    assert got.endswith("…")
    assert not got[:-1].rstrip().endswith(" ")
    # No word is broken: every word in the label is a word of the title.
    words = set(long.replace(",", " ").replace(".", " ").split())
    assert all(w in words for w in got[:-1].split()), got


def test_an_empty_or_markup_only_title_falls_back_to_the_id():
    for title in ("", "   ", "\n\n", "```", "###", None):
        got = AS.session_label(title, "9d19b79871ff4a0e")
        assert got, f"no label for {title!r}"
        assert "9d19b798" in got, got


def test_whitespace_inside_the_line_is_collapsed():
    assert AS.session_label("push    the     branch", "x") == (
        "push the branch")


# -- the two places that must use it ---------------------------------------

def test_the_resume_dropdown_and_the_hint_both_use_it():
    """Asserted on the source: both sites build their text inside
    create_tab's closures, which need a built dashboard to reach."""
    import pathlib

    src = pathlib.Path(AS.__file__).read_text(encoding="utf-8")
    assert src.count("session_label(") >= 3, (
        "the definition plus both call sites")
    # And neither may fall back to the raw title with a slice.
    assert 'str(row.get("title") or "Untitled")[:40]' not in src
    assert 'str(r.get("title") or "a session")' not in src


def test_an_underscore_is_part_of_a_name_not_emphasis():
    """Found by running this over the real 52 titles: stripping `_` as
    markdown emphasis turned "test_calc.py schlaegt fehl" into
    "testcalc.py", a file that does not exist. An emphasis mark lost is
    cosmetic; a filename altered is a wrong answer."""
    got = AS.session_label("test_calc.py schlaegt fehl. Finde den Fehler", "x")
    assert got.startswith("test_calc.py"), got
    assert "_" in AS.session_label("fix user_project_workspace now", "x")
    # A tilde is a home path here far more often than strikethrough.
    assert "~" in AS.session_label("read ~/notes.txt", "x")


def test_a_long_path_does_not_collapse_the_label():
    """Also from the real titles: "Bau in
    tests/fixtures/user_project_workspace/ ein kleines Modul ..." has no
    space inside the budget, so the word boundary left "Bau in" -- which
    names no session. Below half the budget the character cut says more."""
    got = AS.session_label(
        "Bau in tests/fixtures/user_project_workspace/ ein kleines Modul "
        "tagstat.py, das Zeilen zaehlt", "x")
    assert len(got) > AS._LABEL_MAX // 2, got
    assert "tests/fixtures" in got, got
