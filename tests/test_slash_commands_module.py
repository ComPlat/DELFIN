"""Tests for delfin.agent.slash_commands — the markdown-template
slash-command system that powers /commands and the user's own /<name>
shortcuts."""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import slash_commands as sc


@pytest.fixture
def fake_home(monkeypatch, tmp_path):
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    return tmp_path


def _write(d: Path, name: str, body: str) -> Path:
    d.mkdir(parents=True, exist_ok=True)
    p = d / f"{name}.md"
    p.write_text(body, encoding="utf-8")
    return p


# --- discovery --------------------------------------------------------------

def test_discover_picks_up_user_global(fake_home):
    _write(fake_home / ".delfin" / "commands", "review",
           "Review PR $1 and report.")
    cmds = sc.discover_commands(None)
    names = [c.name for c in cmds]
    assert "review" in names


def test_workspace_override_wins(fake_home, tmp_path):
    """When the same command exists in user + workspace, workspace wins."""
    ws = tmp_path / "repo"
    _write(fake_home / ".delfin" / "commands", "x", "USER body")
    _write(ws / ".delfin" / "commands", "x", "PROJECT body")
    cmds = sc.discover_commands(ws)
    x = next(c for c in cmds if c.name == "x")
    assert "PROJECT" in x.body


def test_canonical_commands_dir_also_discovered(fake_home):
    """~/.claude/commands/ is read as a compat layer so users who had
    files there don't lose them — but is the lowest-priority layer."""
    _write(fake_home / ".claude" / "commands", "from-canonical",
           "compat body")
    cmds = sc.discover_commands(None)
    assert any(c.name == "from-canonical" for c in cmds)


def test_invalid_filename_skipped(fake_home):
    d = fake_home / ".delfin" / "commands"
    d.mkdir(parents=True)
    (d / "bad name!.md").write_text("body", encoding="utf-8")
    (d / "good.md").write_text("body", encoding="utf-8")
    cmds = sc.discover_commands(None)
    assert any(c.name == "good" for c in cmds)
    assert not any(c.name.startswith("bad") for c in cmds)


def test_empty_file_skipped(fake_home):
    _write(fake_home / ".delfin" / "commands", "empty", "   \n  \n")
    assert all(c.name != "empty" for c in sc.discover_commands(None))


# --- frontmatter ------------------------------------------------------------

def test_frontmatter_description_used(fake_home):
    _write(
        fake_home / ".delfin" / "commands", "review",
        "---\ndescription: Audit a PR for review\nargument-hint: <pr-number>\n---\n\nReview PR #$1.",
    )
    tpl = sc.get_command("review", None)
    assert tpl is not None
    assert tpl.description == "Audit a PR for review"
    assert tpl.argument_hint == "<pr-number>"
    assert tpl.body == "Review PR #$1."


def test_description_falls_back_to_first_line(fake_home):
    _write(fake_home / ".delfin" / "commands", "explain",
           "Explain the implementation of $ARGUMENTS in plain language.")
    tpl = sc.get_command("explain", None)
    assert tpl is not None
    assert "Explain the implementation" in tpl.description


# --- expand_template --------------------------------------------------------

def test_expand_substitutes_arguments():
    out = sc.expand_template("Review PR $ARGUMENTS for me.", "8423")
    assert out == "Review PR 8423 for me."


def test_expand_substitutes_positional_tokens():
    out = sc.expand_template("file=$1 line=$2", "foo.py 42")
    assert out == "file=foo.py line=42"


def test_expand_zero_is_full_args():
    """``$0`` mirrors POSIX ``$0`` semantics: the whole args string."""
    out = sc.expand_template("got: $0", "first second")
    assert out == "got: first second"


def test_expand_unknown_placeholder_left_intact():
    """Random ``$0.05`` or unknown vars must NOT vanish."""
    out = sc.expand_template("price is $0.05 and $foo bar", "")
    # ``$0`` → "", then ".05" remains; ``$foo`` is unknown → untouched
    assert ".05" in out
    assert "$foo" in out


def test_expand_at_alias_for_arguments():
    out = sc.expand_template("got $@", "1 2 3")
    assert out == "got 1 2 3"


def test_expand_missing_positional_becomes_empty():
    out = sc.expand_template("a=$1 b=$2 c=$3", "first")
    assert out == "a=first b= c="


def test_expand_braced_placeholders():
    out = sc.expand_template("hi ${ARGUMENTS}, bye", "world")
    assert out == "hi world, bye"


# --- render_command end-to-end --------------------------------------------

def test_render_command_via_name(fake_home):
    _write(fake_home / ".delfin" / "commands", "greet",
           "Hello $1, welcome to $2.")
    rendered = sc.render_command("greet", "Alice DELFIN", None)
    assert rendered == "Hello Alice, welcome to DELFIN."


def test_render_command_strips_leading_slash(fake_home):
    """Convenience: get_command should accept both 'name' and '/name'."""
    _write(fake_home / ".delfin" / "commands", "x", "body")
    assert sc.get_command("/x", None) is not None
    assert sc.get_command("x", None) is not None


def test_render_unknown_raises():
    with pytest.raises(KeyError):
        sc.render_command("does-not-exist", "")
