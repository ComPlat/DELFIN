"""Skill discovery and prompt expansion, against the shipped API.

This file was skipped whole for seven weeks. The skip reason named the
mismatch exactly: it called ``parse_skill_text`` / ``find_skill`` /
``format_skill_message`` and a ``Skill`` carrying ``title`` /
``source_path`` / ``slash_command``, none of which exist. The shipped
module parses front-matter (``_parse_frontmatter``), looks up by name
(``get_skill``) and renders an injectable block
(``render_skill_invocation``), and its ``Skill`` carries ``description``
and a ``source`` Path.

Twenty tests were therefore not running over the layer that decides
which playbook text is injected into a prompt. The coverage was written;
only the names were wrong.

Deliberately NOT duplicated here: what ``test_skills_system.py`` already
covers (project/user/pack precedence, the ``skill`` tool's dispatch) and
what ``test_skills_match_the_workspace.py`` covers (domain filtering).
What is left is the parsing and rendering underneath both, which had no
tests of its own.
"""

from __future__ import annotations

import pathlib

import pytest

from delfin.agent import skills


@pytest.fixture
def only_this_workspace(tmp_path, monkeypatch):
    """A workspace whose skills are the ONLY ones discoverable.

    ``discover_skills`` always also scans the packaged directory and the
    user's own, so without this every "no skills" assertion would be a
    statement about the 13 shipped playbooks instead.
    """
    empty = tmp_path / "no_pack"
    empty.mkdir()
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setattr(skills, "_PACK_SKILLS_DIR", empty)
    monkeypatch.setattr(skills.Path, "home", classmethod(lambda cls: home))
    ws = tmp_path / "projekt"
    (ws / ".delfin" / "skills").mkdir(parents=True)
    return ws


def _write(ws: pathlib.Path, name: str, text: str) -> pathlib.Path:
    p = ws / ".delfin" / "skills" / f"{name}.md"
    p.write_text(text, encoding="utf-8")
    return p


def _write_folder_skill(ws: pathlib.Path, dirname: str,
                        text: str) -> pathlib.Path:
    d = ws / ".delfin" / "skills" / dirname
    d.mkdir(parents=True)
    p = d / "SKILL.md"
    p.write_text(text, encoding="utf-8")
    return p


# ---------------------------------------------------------------------------
# Front-matter
# ---------------------------------------------------------------------------

def test_front_matter_is_split_from_the_body():
    meta, body = skills._parse_frontmatter(
        "---\nname: tune\ndescription: Tune it\n---\n\n# Heading\nbody line\n")
    assert meta == {"name": "tune", "description": "Tune it"}
    assert body == "# Heading\nbody line\n"


def test_a_file_without_front_matter_is_all_body():
    """The loader is permissive on purpose: a user's own note-shaped skill
    keeps working."""
    meta, body = skills._parse_frontmatter("# Heading\nbody\n")
    assert meta == {}
    assert body == "# Heading\nbody\n"


def test_an_unterminated_front_matter_block_is_not_eaten():
    """Without the closing fence the whole file is body -- losing it
    silently would delete the playbook the user wrote."""
    text = "---\nname: tune\n\n# Heading\nbody\n"
    meta, body = skills._parse_frontmatter(text)
    assert meta == {}
    assert body == text


def test_quotes_and_comments_in_front_matter():
    meta, _ = skills._parse_frontmatter(
        '---\n# a comment\nname: "tune"\ndescription: \'Tune it\'\n\n'
        'not-a-pair\n---\nbody\n')
    assert meta == {"name": "tune", "description": "Tune it"}


def test_a_single_line_file_is_not_front_matter():
    assert skills._parse_frontmatter("---") == ({}, "---")


# ---------------------------------------------------------------------------
# domains
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("raw,expected", [
    ("office", ("office",)),
    ("office, code", ("office", "code")),
    ("[office, code]", ("office", "code")),
    ("office; code", ("office", "code")),
    ('"Office"', ("office",)),
    ("", ()),
    (None, ()),
    (["office"], ()),
])
def test_domains_are_read_the_way_people_write_them(raw, expected):
    assert skills._parse_domains(raw) == expected


def test_an_unreadable_domain_means_everywhere_not_nowhere():
    """A typo must not be able to hide a skill from every session."""
    sk = skills.Skill(name="x", description="", body="b",
                      source=pathlib.Path("x.md"),
                      domains=skills._parse_domains(object()))
    assert sk.applies_to("office")
    assert sk.applies_to("code")


# ---------------------------------------------------------------------------
# The description a listing shows
# ---------------------------------------------------------------------------

def test_the_first_heading_stands_in_for_a_missing_description(
        only_this_workspace):
    _write(only_this_workspace, "tune", "# Tune CONTROL\n\nbody line\n")
    sk = skills.get_skill("tune", only_this_workspace)
    assert sk is not None
    assert sk.description == "Tune CONTROL"


def test_a_declared_description_wins_over_the_heading(only_this_workspace):
    _write(only_this_workspace, "tune",
           "---\ndescription: One line.\n---\n# Tune CONTROL\nbody\n")
    sk = skills.get_skill("tune", only_this_workspace)
    assert sk is not None
    assert sk.description == "One line."


def test_a_file_with_neither_still_loads(only_this_workspace):
    """An empty description is a worse listing, not a lost skill."""
    _write(only_this_workspace, "tune", "just the body, no heading\n")
    sk = skills.get_skill("tune", only_this_workspace)
    assert sk is not None
    assert sk.description == ""
    assert sk.body == "just the body, no heading"


@pytest.mark.parametrize("heading,expected", [
    ("# Title", "Title"),
    ("### Deep title", "Deep title"),
    ("   # Indented", "Indented"),
])
def test_heading_markers_are_stripped(heading, expected):
    assert skills._first_heading(f"{heading}\nbody") == expected


# ---------------------------------------------------------------------------
# What becomes a skill, and what does not
# ---------------------------------------------------------------------------

def test_the_name_comes_from_the_front_matter(only_this_workspace):
    _write(only_this_workspace, "datei",
           "---\nname: tune-control\n---\n# T\nbody\n")
    assert [s.name for s in skills.discover_skills(only_this_workspace)] == [
        "tune-control"]


def test_a_flat_file_is_named_after_itself(only_this_workspace):
    _write(only_this_workspace, "tune-control", "# T\nbody\n")
    assert [s.name for s in skills.discover_skills(only_this_workspace)] == [
        "tune-control"]


def test_a_folder_skill_is_named_after_its_folder(only_this_workspace):
    _write_folder_skill(only_this_workspace, "tune-control", "# T\nbody\n")
    assert [s.name for s in skills.discover_skills(only_this_workspace)] == [
        "tune-control"]


@pytest.mark.parametrize("bad", ["1tune", "-tune", "tune space", "tü"])
def test_a_name_that_is_not_a_name_is_refused(only_this_workspace, bad):
    """The name reaches a ``/<name>`` command and a tool argument."""
    _write(only_this_workspace, "datei", f"---\nname: {bad}\n---\nbody\n")
    assert skills.discover_skills(only_this_workspace) == []


def test_an_empty_file_is_not_a_skill(only_this_workspace):
    _write(only_this_workspace, "leer", "   \n\n  \n")
    assert skills.discover_skills(only_this_workspace) == []


def test_a_non_markdown_file_is_ignored(only_this_workspace):
    (only_this_workspace / ".delfin" / "skills" / "notes.txt").write_text(
        "# Not a skill\nbody\n", encoding="utf-8")
    assert skills.discover_skills(only_this_workspace) == []


def test_a_top_level_skill_md_is_not_a_skill(only_this_workspace):
    """``SKILL.md`` names the FOLDER form; at the top level it has no
    folder to be named after."""
    _write(only_this_workspace, "SKILL", "# T\nbody\n")
    assert skills.discover_skills(only_this_workspace) == []


def test_discovery_is_sorted_by_name(only_this_workspace):
    for n in ("gamma", "alpha", "beta"):
        _write(only_this_workspace, n, f"# {n}\nbody\n")
    assert [s.name for s in skills.discover_skills(only_this_workspace)] == [
        "alpha", "beta", "gamma"]


def test_the_source_path_is_recorded(only_this_workspace):
    p = _write(only_this_workspace, "tune", "# T\nbody\n")
    sk = skills.get_skill("tune", only_this_workspace)
    assert sk is not None
    assert sk.source == p


def test_the_summary_carries_what_a_listing_needs(only_this_workspace):
    _write(only_this_workspace, "tune",
           "---\ndescription: One line.\ndomains: office\n---\nbody\n")
    sk = skills.get_skill("tune", only_this_workspace)
    assert sk is not None
    assert sk.to_summary() == {
        "name": "tune",
        "description": "One line.",
        "source": str(sk.source),
        "domains": ["office"],
    }


# ---------------------------------------------------------------------------
# Lookup
# ---------------------------------------------------------------------------

def test_an_unknown_name_is_not_a_skill(only_this_workspace):
    _write(only_this_workspace, "tune", "# T\nbody\n")
    assert skills.get_skill("nichts", only_this_workspace) is None


def test_lookup_is_exact(only_this_workspace):
    """``/tune`` and ``/Tune`` are different commands; guessing between
    them would silently run a playbook the user did not name."""
    _write(only_this_workspace, "tune", "# T\nbody\n")
    assert skills.get_skill("TUNE", only_this_workspace) is None
    assert skills.get_skill("tune", only_this_workspace) is not None


def test_an_empty_name_matches_nothing(only_this_workspace):
    _write(only_this_workspace, "tune", "# T\nbody\n")
    assert skills.get_skill("", only_this_workspace) is None
    assert skills.get_skill("   ", only_this_workspace) is None


# ---------------------------------------------------------------------------
# Expansion -- what actually reaches the prompt
# ---------------------------------------------------------------------------

def _skill(body: str = "Do the thing.") -> skills.Skill:
    return skills.Skill(name="tune", description="d", body=body,
                        source=pathlib.Path("tune.md"))


def test_the_rendered_block_names_the_skill_and_is_delimited():
    out = skills.render_skill_invocation(_skill())
    assert out.startswith("--- Skill: tune ---")
    assert out.rstrip().endswith("--- end skill ---")
    assert "Do the thing." in out


def test_the_body_is_carried_whole():
    """The body IS the instruction; a renderer that trimmed it would drop
    steps out of the middle of a playbook."""
    body = "1. erste Zeile\n\n2. zweite Zeile\n   eingerückt\n"
    assert body.strip() in skills.render_skill_invocation(_skill(body))


def test_arguments_are_labelled_when_given():
    out = skills.render_skill_invocation(_skill(), args="molekuel.xyz")
    assert "Arguments: molekuel.xyz" in out


def test_no_arguments_means_no_argument_line():
    """An empty ``Arguments:`` line reads as "none were given" and the
    model answers about it."""
    assert "Arguments:" not in skills.render_skill_invocation(_skill())
    assert "Arguments:" not in skills.render_skill_invocation(
        _skill(), args="   ")


def test_a_discovered_skill_renders_its_own_file(only_this_workspace):
    """End to end: the file on disk is what reaches the prompt."""
    _write(only_this_workspace, "tune",
           "---\nname: tune\n---\n# Tune CONTROL\n\nSetze SCF auf tight.\n")
    sk = skills.get_skill("tune", only_this_workspace)
    assert sk is not None
    out = skills.render_skill_invocation(sk, args="job1")
    assert "Setze SCF auf tight." in out
    assert "Arguments: job1" in out


# ---------------------------------------------------------------------------
# The shipped set is still reachable
# ---------------------------------------------------------------------------

def test_the_packaged_playbooks_are_discoverable():
    """Without the pack directory in the search order the shipped skills
    were undiscoverable by the tool that offers them."""
    names = {s.name for s in skills.discover_skills()}
    assert "tune-control" in names, sorted(names)
