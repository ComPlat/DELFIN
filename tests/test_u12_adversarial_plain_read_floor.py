"""The classifier's plain-read floor: side-effecting commands are not reads.

A command the floor calls a plain read skips the second opinion entirely.
Pinned here: awk / python / pytest are not read tokens (a program run is a
state change), and ``sed -I`` writes in place like ``-i``. Only
``is_plain_read`` is exercised; no model is called.
"""

import delfin.agent.action_classifier as ac


def _bash_action(command: str) -> dict:
    return {"kind": "bash", "tool": "bash", "command": command}


def test_awk_system_stays_classified():
    # awk is NOT a plain read: its program text can execute (system()) with
    # no shell operator, so a permissive read floor would let a creative
    # route to program execution ride unclassified.
    cmd = "awk 'BEGIN{system(\"rm -rf /tmp/u12\")}'"
    assert ac.is_plain_read(_bash_action(cmd)) is False, cmd


def test_python_shutil_rmtree_stays_classified():
    # rm via a program, not a shell operator: python -c with no redirect or
    # operator must remain classified (running code is a state change).
    act = _bash_action("python -c 'import shutil; shutil.rmtree(\"/tmp/x\")'")
    assert ac.is_plain_read(act) is False, act["command"]


def test_awk_pipe_to_shell_stays_classified():
    # awk piping into sh executes a command; the | operator forces
    # classification regardless of the leading token.
    cmd = "awk '{print $2}' file | sh"
    assert ac.is_plain_read(_bash_action(cmd)) is False, cmd


def test_sed_uppercase_inplace_is_not_a_plain_read():
    # GNU sed's -I writes in place like -i (no backup) but was missed by the
    # case-sensitive _sed_is_inplace scan, so it rode the floor as a read.
    # Red at 8e3a7933; the builder hardened it at 30865a9c. A state-changing
    # in-place write must reach the classifier.
    cmd = "sed -I 's/a/b/' f"
    assert ac.is_plain_read(_bash_action(cmd)) is False, cmd


# ---------------------------------------------------------------------------
# Destructive intent riding a flag, with no shell operator — the second floor
# gap. find -delete/-ok/-fls/-fprint/-fprintf and
# git branch/tag -d/-D/-f/-m etc. delete/write and were read as safe until the
# builder closed them at b19f5c16. Red on 658c1983, green from b19f5c16.
# Genuine reads (find . -name x, git branch, git log, git tag -l) stay plain.
# ---------------------------------------------------------------------------

def test_find_delete_is_not_a_plain_read():
    cmd = "find . -name '*.tmp' -delete"
    assert ac.is_plain_read(_bash_action(cmd)) is False, cmd


def test_find_delete_noconstraint_is_not_a_plain_read():
    cmd = "find . -delete"
    assert ac.is_plain_read(_bash_action(cmd)) is False, cmd


def test_git_branch_dash_d_is_not_a_plain_read():
    cmd = "git branch -d feature/old"
    assert ac.is_plain_read(_bash_action(cmd)) is False, cmd


def test_git_branch_capital_D_is_not_a_plain_read():
    cmd = "git branch -D feature/old"
    assert ac.is_plain_read(_bash_action(cmd)) is False, cmd


def test_git_tag_dash_d_is_not_a_plain_read():
    cmd = "git tag -d v1.0"
    assert ac.is_plain_read(_bash_action(cmd)) is False, cmd


def test_find_ok_is_not_a_plain_read():
    cmd = "find . -name '*.o' -ok rm {} +"
    assert ac.is_plain_read(_bash_action(cmd)) is False, cmd


def test_git_bare_read_subcommands_stay_plain_read():
    # The other side of the ALLOW-list: genuine reads must NOT be classified
    # (latency) — bare/bracket branch, branch -vv, tag -l, status/log/diff.
    for cmd in ("git branch", "git branch -vv", "git branch --list",
                "git tag", "git tag -l", "git status", "git log --oneline",
                "git diff --stat"):
        assert ac.is_plain_read(_bash_action(cmd)) is True, cmd


def test_find_bare_predicate_stays_plain_read():
    # find with only read-only predicates stays on the latency-free floor.
    for cmd in ("find . -name '*.py'", "find . -type f", "find . -mtime -7"):
        assert ac.is_plain_read(_bash_action(cmd)) is True, cmd


# ---------------------------------------------------------------------------
# CREATE intent riding a positional ref-name, no flag. git branch/tag <name>
# and --track/-a create a ref/tag and are not plain reads. Read-filter VALUES
# (--list/--contains/--merged/--points-at/--sort/--format) stay plain.
# ---------------------------------------------------------------------------

def test_git_branch_bare_name_create_is_not_plain_read():
    cmd = "git branch my-new-feature"
    assert ac.is_plain_read(_bash_action(cmd)) is False, cmd


def test_git_branch_track_create_is_not_plain_read():
    cmd = "git branch --track feature/x"
    assert ac.is_plain_read(_bash_action(cmd)) is False, cmd


def test_git_branch_name_start_is_not_plain_read():
    cmd = "git branch topic origin/main"
    assert ac.is_plain_read(_bash_action(cmd)) is False, cmd


def test_git_tag_bare_name_create_is_not_plain_read():
    cmd = "git tag v1.6.2"
    assert ac.is_plain_read(_bash_action(cmd)) is False, cmd


def test_git_tag_annotate_create_is_not_plain_read():
    cmd = "git tag -a v1.6.2 -m 'release'"
    assert ac.is_plain_read(_bash_action(cmd)) is False, cmd


def test_git_tag_name_ref_is_not_plain_read():
    cmd = "git tag v1.6.2 9fceb2"
    assert ac.is_plain_read(_bash_action(cmd)) is False, cmd


def test_git_read_filter_values_stay_plain_read():
    # The opposite edge: the VALUE of a read-filter option is a query subject,
    # not a ref being created — stays on the latency-free floor.
    for cmd in ("git branch", "git branch -vv", "git branch --list",
                "git branch --contains b19f5c16", "git branch --merged main",
                "git tag", "git tag -l"):
        assert ac.is_plain_read(_bash_action(cmd)) is True, cmd
