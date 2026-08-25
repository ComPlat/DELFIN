"""What the user sees before typing anything, and after every turn.

Taken from a real session on a cluster login node. The complaints were
not aesthetic: an absolute `$HOME` path pushed the part that identifies
the directory off the readable width, a status line printed `branch=`
with nothing after it, and a multi-line task report arrived as one
wrapped paragraph with the checkbox glyphs run together.

Each of those is a mechanism describing itself wrongly, which is the
same class of defect as a gate that does not gate — it just fails in the
reader rather than in the code.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import cli as agent_cli
from delfin.agent import repl_render as rr
from delfin.agent import status_line as sl


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

def test_home_is_written_as_a_tilde(monkeypatch, tmp_path):
    home = tmp_path / "home" / "someone"
    (home / "work" / "thing").mkdir(parents=True)
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))

    assert agent_cli._tilde(home / "work" / "thing") == "~/work/thing"
    assert agent_cli._tilde(home) == "~"


def test_a_path_outside_home_is_left_alone(monkeypatch, tmp_path):
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path / "h"))
    assert agent_cli._tilde("/opt/data") == "/opt/data"


def test_a_sibling_of_home_is_not_shortened(monkeypatch, tmp_path):
    """`/home/someone-else` starts with `/home/someone`, and is not in it."""
    home = tmp_path / "home" / "someone"
    home.mkdir(parents=True)
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))
    other = str(home) + "-else/project"
    assert agent_cli._tilde(other) == other


# ---------------------------------------------------------------------------
# The banner
# ---------------------------------------------------------------------------

class _Git:
    def __init__(self, is_repo=True, branch="main", dirty=(), unborn=False):
        self.is_repo = is_repo
        self.branch = branch
        self.dirty = dirty
        self.unborn = unborn


class _Report:
    def __init__(self, git):
        self.git = git
        self.granted_dirs = ()
        self.read_dirs = ()


class _Engine:
    client = type("C", (), {"model": "kit.qwen3.5-397b-A17b"})()
    provider = "kit"
    mode = "solo"
    session_id = "0123456789abcdef0123456789abcdef"
    kit_permissions = None


def _banner(git, **kw):
    return agent_cli._startup_banner(_Engine(), _Report(git),
                                     Path("/tmp/project"), **kw)


def test_a_directory_without_git_says_so_in_words():
    text = _banner(_Git(is_repo=False, branch=""))
    assert "not a git repository" in text
    assert "branch" not in text, (
        "an empty field reads as a lookup that failed, not as a fact")


def test_a_detached_head_is_not_reported_as_no_repository():
    text = _banner(_Git(is_repo=True, branch=""))
    assert "detached HEAD" in text
    assert "not a git repository" not in text


# ---------------------------------------------------------------------------
# What `git init && delfin-agent` actually reports — found by running it
# ---------------------------------------------------------------------------

def test_a_repository_with_no_commits_is_not_a_detached_head():
    """`git init` then straight in is an ordinary way to start.

    On that repository `git rev-parse --abbrev-ref HEAD` exits 128, so a
    branch read only that way comes back empty and the banner announced
    a detached HEAD — a different and alarming state.
    """
    text = _banner(_Git(branch="main", unborn=True))
    assert "no commits yet" in text
    assert "detached" not in text


def _real_git_info(path):
    from delfin.agent import launch_guard as lg
    return lg._git_info(Path(path))


def test_a_fresh_repository_reports_its_branch_and_its_unbornness(tmp_path):
    import subprocess
    subprocess.run(["git", "init", "-q", "."], cwd=tmp_path, check=True)
    info = _real_git_info(tmp_path)
    assert info.is_repo is True
    assert info.branch, "an unborn branch still has a name"
    assert info.unborn is True


def test_a_repository_with_a_commit_is_not_unborn(tmp_path):
    import subprocess
    run = lambda *a: subprocess.run(a, cwd=tmp_path, check=True,
                                    capture_output=True)
    run("git", "init", "-q", ".")
    run("git", "config", "user.email", "t@example.invalid")
    run("git", "config", "user.name", "t")
    (tmp_path / "f.txt").write_text("x")
    run("git", "add", "f.txt")
    run("git", "commit", "-qm", "first")
    info = _real_git_info(tmp_path)
    assert info.unborn is False
    assert info.branch


def test_the_agents_own_state_directory_is_not_the_users_loose_work(tmp_path):
    """`.delfin/` is created by the agent, in the workspace, on turn one.

    Counting it as uncommitted work told the user they had changes in a
    directory they had never touched — and it did so in the same sentence
    that warns their work and the agent's will be hard to tell apart.
    """
    import subprocess
    subprocess.run(["git", "init", "-q", "."], cwd=tmp_path, check=True)
    (tmp_path / ".delfin").mkdir()
    (tmp_path / ".delfin" / "session_tasks.json").write_text("{}")
    (tmp_path / "mine.txt").write_text("x")

    info = _real_git_info(tmp_path)
    assert "mine.txt" in info.dirty
    assert not [p for p in info.dirty if p.startswith(".delfin")], info.dirty


def test_the_session_id_is_shortened_to_what_a_person_types():
    text = _banner(_Git())
    assert "session    01234567" in text
    assert _Engine.session_id not in text, (
        "the full id made the widest line of the banner the least useful")


def test_the_hint_line_names_the_keys_that_exist():
    """The line used to name Ctrl+C alone.

    Esc ends a turn and Shift+Tab cycles the approval mode; both are
    implemented in the key layer, and neither was mentioned anywhere the
    user would look.
    """
    from delfin.agent import repl_keys as rk

    text = _banner(_Git())
    hint = text.splitlines()[-1]
    assert "esc" in hint
    assert "shift+tab" in hint
    assert "/help" in hint
    # The keys named are the ones the decoder actually produces.
    assert rk.INTERRUPT and rk.CYCLE_MODE


def test_a_granted_directory_appears_on_screen():
    """A grant nobody can see is a grant nobody can audit."""
    report = _Report(_Git())
    report.granted_dirs = (Path("/srv/lib-a"), Path("/srv/lib-b"))
    report.read_dirs = (Path("/srv/reference"),)
    text = agent_cli._startup_banner(_Engine(), report, Path("/tmp/project"))
    assert "/srv/lib-a" in text and "/srv/lib-b" in text
    assert "readable" in text and "/srv/reference" in text


# ---------------------------------------------------------------------------
# The status line
# ---------------------------------------------------------------------------

def test_the_default_line_omits_a_branch_it_does_not_have(tmp_path):
    ctx = sl.StatusContext(workspace=tmp_path, mode="plan", tokens=0)
    ctx.branch = ""
    out = sl.render_status_line(ctx)
    assert "branch=" not in out
    assert "mode=plan" in out


def test_the_default_line_still_names_a_branch_it_does_have(tmp_path):
    ctx = sl.StatusContext(workspace=tmp_path, mode="plan", tokens=12)
    ctx.branch = "feature/x"
    assert "branch=feature/x" in sl.render_status_line(ctx)


def test_a_users_own_template_still_gets_the_empty_string(tmp_path):
    """The truth goes to a template the user wrote; only OUR default hides it."""
    import json
    (tmp_path / ".delfin").mkdir()
    (tmp_path / ".delfin" / "settings.json").write_text(
        json.dumps({"statusLine": {"template": "[{branch}]"}}))
    ctx = sl.StatusContext(workspace=tmp_path, mode="plan", tokens=0)
    ctx.branch = ""
    assert sl.render_status_line(ctx) == "[]"


def test_configured_and_not_configured_are_distinguishable(tmp_path):
    import json
    assert sl.has_custom_status_line(tmp_path) is False
    (tmp_path / ".delfin").mkdir()
    (tmp_path / ".delfin" / "settings.json").write_text(
        json.dumps({"statusLine": {"template": "x"}}))
    assert sl.has_custom_status_line(tmp_path) is True


def test_the_terminal_prints_it_only_when_it_was_configured():
    """The docstring said "if they configured one" and the code always did.

    `render_status_line` falls back to a built-in template, so every turn
    ended with a line repeating the mode and branch the banner already
    carries.
    """
    import inspect
    from delfin.agent import repl
    src = inspect.getsource(repl.TerminalAgent._report_status)
    assert "has_custom_status_line" in src


# ---------------------------------------------------------------------------
# Notices
# ---------------------------------------------------------------------------

def test_a_multi_line_notice_keeps_its_lines():
    """The open-task report is a heading plus one line per item.

    Collapsing whitespace turned it into a paragraph in which the
    checkbox glyphs ran together and nothing could be counted.
    """
    text = ("⚠ This turn ended with 3 open task(s):\n"
            "  ☐ task 2 (id 7) refactor loader\n"
            "  ☐ task 3 (id 8) write tests\n"
            "  ⛔ task 6 (id 11) deploy — waiting on a key")
    out = rr.notice_line(text, theme=rr.Theme())
    lines = out.splitlines()
    assert len(lines) == 4, out
    assert lines[0].startswith("! ⚠")
    assert lines[1].strip().startswith("☐")
    assert lines[3].strip().startswith("⛔")


def test_a_single_line_notice_is_unchanged():
    assert rr.notice_line("retrying (2/3)", theme=rr.Theme()) == \
        "! retrying (2/3)"


def test_columns_inside_a_line_are_still_collapsed():
    """A notice built from tool output must not paint its own table."""
    out = rr.notice_line("a\t\t\tb     c", theme=rr.Theme())
    assert out == "! a b c"


def test_a_notice_cannot_move_the_cursor():
    out = rr.notice_line("done\rHACKED\n\x1b[2Jgone", theme=rr.Theme())
    assert "\r" not in out and "\x1b[" not in out


def test_an_empty_notice_prints_nothing():
    assert rr.notice_line("   \n\n  ", theme=rr.Theme()) == ""


# ---------------------------------------------------------------------------
# The flag whose consumer had no producer
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("value", ["auto", "always", "never"])
def test_colour_can_actually_be_chosen(value):
    args = agent_cli.build_parser().parse_args(["chat", "--color", value])
    assert args.color == value


def test_the_repl_reads_the_flag_the_parser_declares():
    import inspect
    from delfin.agent import repl
    assert "color" in inspect.signature(repl.ReplOptions).parameters or \
        "color" in repl.ReplOptions.__dataclass_fields__
