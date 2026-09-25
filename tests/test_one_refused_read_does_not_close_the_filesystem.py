"""One refused read does not close the filesystem.

A refused read is remembered and holds against every later command, and
it used to be matched as a plain substring of the command line. Measured
in a supervised run: an agent asked to read '/' (it meant its own
tests/), was refused, and from then on every command with a slash in it
was blocked -- its own test runner included. A refused '/tests' blocked
every path containing '/tests', and a refused file directly in the home
directory extended to the home directory, which holds the workspace.

The refusal must still hold -- the other half of every case below.
"""

from __future__ import annotations

from pathlib import Path

from delfin.agent.api_client import _bash_reads_denied_path as refused


def test_a_refused_root_blocks_only_the_root():
    denied = {"/"}
    assert refused("ls /", denied)
    for cmd in ("/opt/tools/gate tests/test_x.py | tail -5",
                "cat /etc/hostname",
                "ls tests/"):
        assert refused(cmd, denied) == "", cmd


def test_a_refused_top_level_directory_is_a_path_not_a_substring():
    denied = {"/tests"}
    assert refused("ls /tests", denied)
    assert refused("cat /tests/a.py", denied)
    assert refused("cat /repo/delfin/tests/test_x.py", denied) == ""


def test_a_refused_file_in_the_home_does_not_close_the_home():
    home = str(Path.home())
    denied = {home + "/.bashrc"}
    assert refused("cat ~/.bashrc", denied)
    assert refused("cat $HOME/.bashrc", denied)
    assert refused(f"ls {home}/project/tests", denied) == "", (
        "the home is where the workspace lives; it is never opened by a "
        "read grant, so a refusal does not close it either")


def test_a_deep_refusal_still_holds_in_every_spelling():
    denied = {"/data/secret/parser.py"}
    for cmd in ("cat /data/secret/parser.py",
                "xargs -a/data/secret/parser.py echo",
                "python3 -c \"print(open('/data/secret/parser.py').read())\"",
                "cat file:///data/secret/parser.py",
                "ls -la /data/secret",
                "grep -rn x /data/secret/",
                "cat /data/secret/*"):
        assert refused(cmd, denied), cmd


def test_a_neighbour_with_a_longer_name_is_not_the_refused_one():
    denied = {"/data/secret/parser.py"}
    assert refused("cat /data/secret2/x", denied) == ""
    assert refused("cat /data/public/notes.txt", denied) == ""


def test_a_refused_directory_that_holds_the_workspace_is_refused_as_itself():
    # Night run 2026-09-25: the home directory was refused, and every
    # command naming a path under it -- the workspace, the session's own
    # test gate -- was blocked from then on.
    home = "/data/home/user"
    ws = home + "/software/repo/.delfin/worktrees/s6"
    denied = {home}
    roots = [ws]
    assert refused(f"ls {home}", denied, roots=roots)
    assert refused(f"ls {home}/", denied, roots=roots)
    for cmd in (f"{home}/software/repo/.delfin/gate-tools/gate tests/t.py",
                f"cat {ws}/README.md",
                "ls tests/"):
        assert refused(cmd, denied, roots=roots) == "", cmd


def test_a_refused_directory_without_a_granted_root_still_covers_its_contents():
    denied = {"/data/secret"}
    assert refused("cat /data/secret/key.txt", denied, roots=["/work/ws"])


def test_a_refused_read_of_a_program_does_not_stop_running_it():
    # Night run 2026-09-25: refusing to let a session READ its test gate
    # blocked every later RUN of the gate.
    gate = "/data/repo/.delfin/gate-tools/gate"
    denied = {gate}
    assert refused(f"{gate} tests/test_x.py -q", denied) == ""
    assert refused(f"cat {gate}", denied)
    assert refused(f"bash {gate}", denied)
    assert refused(f"ls -l; sed -n 1,5p {gate}", denied)
