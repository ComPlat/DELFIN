"""A refusal reads the command, not the words inside it.

Driven on 2026-09-17 by three sessions building a feature together:

- The deny pattern for ending the machine matched the WORD anywhere in the
  line. It refused `grep -rn "shutdown\\|def stop" engine.py` and then a
  `git commit -m "... on shutdown"` — reading the commit message and costing
  that session its commit.
- A session working in its own git worktree touched the checkout the
  worktree belongs to and was told the path was "the archive of stored
  calculations" and to "COPY it into calc or agent_workspace": the refusal
  was right, the reason and the advice were not.
"""
import pytest

from delfin.agent.api_client import KitToolPermissions, _read_only_reason


@pytest.fixture
def perms(tmp_path):
    return KitToolPermissions(workspace=str(tmp_path))


@pytest.mark.parametrize("cmd", [
    "shutdown -h now", "sudo shutdown -r now", "/sbin/reboot", "echo hi; reboot",
    "make && poweroff", "$(reboot)", "`reboot`", "init 0", "init 6", "halt",
    "python -c \"import os; os.system('reboot')\"",
])
def test_ending_the_machine_is_still_refused(perms, cmd):
    assert perms.matches_bash_deny(cmd)


@pytest.mark.parametrize("cmd", [
    'grep -rn "shutdown\\|def stop|atexit" delfin/agent/engine.py',
    'git commit -m "Add session report written on shutdown"',
    'git add x && git commit -m "wire the report into the shutdown path"',
    'grep -n "shutdown" delfin/agent/cli.py',
    'echo refresh the report on shutdown',
])
def test_the_word_inside_a_string_is_not_a_command(perms, cmd):
    assert perms.matches_bash_deny(cmd) is None


def test_a_repository_checkout_says_so_and_not_calc(tmp_path):
    repo = tmp_path / "repo"
    (repo / ".git").mkdir(parents=True)
    worktree = tmp_path / "wt"
    worktree.mkdir()
    perms = KitToolPermissions(workspace=str(worktree))
    reason = _read_only_reason(repo / "delfin" / "x.py", perms)
    assert "repository checkout" in reason and str(repo) in reason
    assert "archive of stored calculations" not in reason
    assert "workspace directory" in reason          # what to do instead


def test_a_path_that_is_no_repository_keeps_the_archive_wording(tmp_path):
    plain = tmp_path / "calc" / "JOB-1"
    plain.mkdir(parents=True)
    perms = KitToolPermissions(workspace=str(tmp_path / "ws"))
    reason = _read_only_reason(plain / "out.txt", perms)
    assert "archive of stored calculations" in reason


def test_the_workspace_itself_is_not_called_a_foreign_checkout(tmp_path):
    ws = tmp_path / "ws"
    (ws / ".git").mkdir(parents=True)
    perms = KitToolPermissions(workspace=str(ws))
    assert "repository checkout" not in _read_only_reason(ws / "f.py", perms)
