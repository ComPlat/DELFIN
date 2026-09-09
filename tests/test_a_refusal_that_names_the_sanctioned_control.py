"""It tried the control three times in sixteen seconds. All three refused.

The most expensive thing found tonight, and it was found by reading what
a live run left behind rather than by reading code. Asked whether the
last commit broke a red test -- the task written that afternoon to
measure exactly this -- kit.deepseek-v4-flash reached for the control
three times:

    18:35:17  cd ensemble_tools && git stash -u; git checkout <ref> -- .
    18:35:24  git checkout <ref> -- test_weights.py energies.py && pytest
    18:35:33  git worktree add /tmp/ctl-weights <ref>

Every one came back "not on the auto-allow list", with no alternative
named. It stopped trying and answered from the diff, and the rubric --
mine, written that afternoon -- failed it for not running a control.

The framework asked for the work in the integrity addendum, refused all
three spellings of it, and then scored it down for not doing it. That is
the shape this project keeps finding, and it is always more expensive
than it looks: 33% pass rate, attributed to the model.

The refusals themselves are RIGHT. `git stash` and `git checkout <ref> --
<paths>` overwrite the user's working tree, and this project's git rules
name both as destructive. `git worktree add` destroys nothing, but its
path argument is a whole checkout written wherever it is pointed, and
_bash_write_targets cannot see that path -- so it stays off the list too.

What was missing is the sentence saying what yes looks like. Same lesson
as test_a_blocked_interpreter_names_the_allowed_way, which was written
for the same failure one tool over: a gate that says no without naming
the sanctioned spelling costs the turn AND the behaviour.
"""

from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path

import pytest

import delfin.agent.api_client as A


@pytest.fixture
def repo():
    with tempfile.TemporaryDirectory(prefix="ctl-") as tmp:
        d = Path(tmp)
        subprocess.run(["git", "init", "-q"], cwd=str(d), check=True)
        yield d


def _err(cmd: str, repo: Path) -> str:
    perms = A.KitToolPermissions(mode="default", workspace=str(repo))
    out = A._doc_executor.execute(
        "bash", {"command": cmd, "description": "d"}, perms)
    return str(json.loads(out).get("error", ""))


# ---------------------------------------------------------------------------
# The three spellings from the recorded run
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cmd", [
    "git worktree add /tmp/ctl-weights 6714d47~1 2>&1 | tail -5",
    "cd ensemble_tools && git stash -u 2>/dev/null; "
    "git checkout 6714d47~1 -- . 2>&1; python -m pytest test_weights.py -q",
    "git checkout 6714d47~1 -- test_weights.py energies.py export.py "
    "&& python -m pytest test_weights.py -q",
])
def test_the_refusal_names_the_tool_that_does_it(cmd, repo):
    err = _err(cmd, repo)
    assert "not on the auto-allow list" in err
    assert "enter_worktree" in err
    assert "base_ref" in err


def test_it_says_why_that_is_not_a_workaround(repo):
    """Without this the hint reads as a loophole and the next reader
    deletes it -- which is exactly what happened to the last one."""
    err = _err("git worktree add /tmp/x HEAD~1", repo)
    assert "cannot destroy work" in err


def test_the_hint_does_not_fire_on_unrelated_git(repo):
    """A hint attached to everything is noise, and noise is what people
    learn to skip."""
    err = _err("git push origin main", repo)
    assert "not on the auto-allow list" in err
    assert "enter_worktree" not in err


# ---------------------------------------------------------------------------
# ...and the read-only half that was refused for no reason at all
# ---------------------------------------------------------------------------

def test_listing_worktrees_is_a_read(repo):
    """`git worktree list` changes nothing and answers the question the
    agent was asking -- do I already have a control tree? It was refused
    while `git stash list` beside it was allowed."""
    perms = A.KitToolPermissions(mode="default", workspace=str(repo))
    out = json.loads(A._doc_executor.execute(
        "bash", {"command": "git worktree list", "description": "d"}, perms))
    assert out.get("exit_code") == 0, out


@pytest.mark.parametrize("cmd", [
    "git worktree add /tmp/x HEAD",
    "git worktree remove /tmp/x",
    "git worktree prune",
])
def test_the_other_worktree_subcommands_are_not_allowed(cmd, repo):
    """add writes a whole checkout wherever it is pointed and the write
    gate cannot see that path; remove and prune delete one."""
    assert not A.KitToolPermissions(
        workspace=repo, mode="default").matches_bash_auto_allow(cmd), cmd


# ---------------------------------------------------------------------------
# Cleaning up after itself
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cmd", [
    "rm -f _tagreport_check.py",
    "rm _verify_export.py && ls -la",
])
def test_removing_its_own_scratch_file_names_undo_changes(cmd, repo):
    """Both from one recorded run: files the agent had written a minute
    earlier to check its own work. Both refusals left the scratch file in
    the user's directory. `rm` cannot tell whose file it is -- which is
    why it is not on the list -- and undo_changes can."""
    err = _err(cmd, repo)
    assert "undo_changes" in err
    assert "this session" in err


def test_a_recursive_delete_is_still_refused_outright(repo):
    """The deny-list runs first and must keep running first: this one is
    not a missing hint, it is the thing the deny-list exists for."""
    err = _err("rm -rf build", repo)
    assert "deny-pattern" in err
    assert "undo_changes" not in err


# ---------------------------------------------------------------------------
# Writing a file through the shell
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cmd", [
    "cat > run.py << 'EOF'\nprint(1)\nEOF",
    "tee out.txt",
    "printf 'x' > cfg.json",
])
def test_a_shell_write_names_write_file(cmd, repo):
    """`echo x > f` happens to run -- echo is on the list and the redirect
    goes through the write gate -- while `cat > f << EOF` and `tee f` did
    not, and the refusal named nothing. Both spellings were observed in
    one recorded run, writing a launcher the agent had just been asked to
    build."""
    err = _err(cmd, repo)
    assert "write_file" in err
    assert "pre-image" in err


def test_the_reason_is_the_journal_not_the_rule(repo):
    """write_file is the better spelling, not merely the allowed one. If
    the hint said only "that is not allowed" the next reader would widen
    the list instead."""
    err = _err("tee out.txt", repo)
    assert "undo_changes cannot take it back" in err


def test_a_refusal_that_writes_nothing_gets_no_write_hint(repo):
    err = _err("curl https://example.com", repo)
    assert "not on the auto-allow list" in err
    assert "write_file" not in err


def test_the_hints_do_not_collide(repo):
    """Four hints share one branch now. Each command must get its own."""
    assert "enter_worktree" in _err("git worktree add /tmp/x HEAD~1", repo)
    assert "undo_changes" in _err("rm -f scratch.py", repo)
    assert "write_file" in _err("tee out.txt", repo)
    assert "py_compile" in _err('python3 -c "import os; os.getcwd()"', repo)
