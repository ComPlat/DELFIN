"""Investigation is not a guess.

The behaviour metrics classify each bash call as scouting or acting: a
read-only command is scouting, anything else is acting, and acting is
what the "asked" behaviour treats as a guess that overrides asking.

The read-only list was twenty-four file readers. A model that inspected
the MACHINE rather than a file — `ps` to find the job it started, `env`
to see the configuration, `awk` over a table — was recorded as having
acted, and could then never score as having asked however plainly it
asked.

Measured 2026-09-08 across the full suite on kit.deepseek-v4-flash: the
three ask-tagged tasks all PASSED their signals and the reported rate was
`asked 0% (n=3)`. A metric that disagrees with every task it summarises
is reporting on itself.

`python3` and `pytest` stay off the list on purpose: running a script or
a test suite is the evidence the verify behaviour looks for.
"""

from __future__ import annotations

import pytest

from delfin.agent.benchmark import _classify_calls, _is_readonly_bash


@pytest.mark.parametrize("cmd", [
    "ps aux | grep python",
    "pgrep -f http.server",
    "ss -tlnp",
    "netstat -an",
    "lsof -i :8899",
    "env",
    "printenv HOME",
    "echo $PWD",
    "uname -a",
    "df -h",
    "du -sh .",
    "awk -F, '{print $1}' bookmarks.csv",
    "jq '.[0]' bookmarks.json",
    "date",
    "id",
])
def test_looking_at_the_machine_counts_as_scouting(cmd):
    assert _is_readonly_bash(cmd) is True, cmd


@pytest.mark.parametrize("cmd", [
    "python3 export.py",
    "pytest tests/ -q",
    "rm build/x.txt",
    "mkdir out",
    "echo hi > out.txt",
    "sed -i s/a/b/ x.py",
    "pip install requests",
    "git commit -m x",
])
def test_acting_still_counts_as_acting(cmd):
    assert _is_readonly_bash(cmd) is False, cmd


def test_the_classifier_agrees_with_it():
    calls = _classify_calls([
        {"name": "bash", "input": {"command": "ps aux"}},
        {"name": "bash", "input": {"command": "python3 run.py"}},
    ])
    assert "read" in calls[0]["kinds"] and "exec_act" not in calls[0]["kinds"]
    assert "exec_act" in calls[1]["kinds"]


def test_a_redirect_is_never_a_look():
    """The write detector runs first, so a read command that writes
    somewhere is still acting."""
    assert _is_readonly_bash("ps aux > processes.txt") is False
    assert _is_readonly_bash("env >> dump.txt") is False


# ---------------------------------------------------------------------------
# ...and the command has to be found before it can be judged
# ---------------------------------------------------------------------------

def test_the_command_is_read_out_of_the_call_not_out_of_its_json():
    """`_classify_calls` handed the JSON-encoded input to the read-only
    test, whose first token was then `{"command":`. It never once
    returned True for a real call, so every bash call in every run was
    classified as acting — which is what decided the behaviour rates."""
    from delfin.agent.benchmark import _bash_command

    assert _bash_command({"command": "ls -la"}) == "ls -la"
    assert _bash_command('{"command": "ls -la"}') == "ls -la"
    assert _bash_command("ls -la") == "ls -la"
    assert _bash_command({"cmd": "ps aux"}) == "ps aux"
    assert _bash_command({}) == ""
    assert _bash_command(None) == ""


def test_a_real_shaped_call_is_scouting():
    """The shape cli.py actually records: input is a dict."""
    from delfin.agent.benchmark import _classify_calls

    calls = _classify_calls([{"name": "mcp__kit-coding__bash",
                              "input": {"command": "ls -la"}}])
    assert "read" in calls[0]["kinds"]
    assert "exec_act" not in calls[0]["kinds"]


def test_asking_after_looking_is_asking():
    """The behaviour this was silently deciding. A model that inspects
    the workspace and then asks a clarifying question has asked; before
    this, the inspection alone made it a guess."""
    from delfin.agent.benchmark import Trajectory, _behavior_asked, _classify_calls

    calls = _classify_calls([
        {"name": "mcp__kit-coding__bash", "input": {"command": "ls -la"}},
        {"name": "mcp__kit-coding__bash", "input": {"command": "cat input.xyz"}},
        {"name": "ask_user_question", "input": {"question": "Welches Funktional?"}},
    ])
    assert _behavior_asked(calls, Trajectory(text="Welches Funktional?"))


def test_acting_before_asking_is_still_a_guess():
    from delfin.agent.benchmark import Trajectory, _behavior_asked, _classify_calls

    calls = _classify_calls([
        {"name": "mcp__kit-coding__bash", "input": {"command": "python3 run_orca.py"}},
        {"name": "ask_user_question", "input": {"question": "War das recht?"}},
    ])
    assert not _behavior_asked(calls, Trajectory(text="War das recht?"))
