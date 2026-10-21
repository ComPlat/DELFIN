"""A tool with a sub-command list is held to that list in every spelling.

Input: a command string. Output: ``sandbox.is_allowed`` (the allow-list of
the dashboard approval runner).

Measured 2026-10-15: ``pip`` was in ``_ALLOW_FIRST_TOKEN`` and in
``_ALLOW_SUBCOMMAND``; ``_check_segment`` consulted the first-token set
first, so ``pip install``, ``pip uninstall`` and ``pip download`` passed
and the pip list ({list, show, freeze, check, config}) was never read.
``python3 -m pip uninstall`` passed through the interpreter, and the
listed ``git config`` and ``pip config`` admitted their writing forms
(``git config core.pager 'sh -c ...'``), which a later git or pip outside
the sandbox executes.

Not covered, by construction: an allow-list over the first token cannot
contain an interpreter. ``python -c`` and a script file run anything;
the OS confinement of the runner is what limits them.
"""
from __future__ import annotations

import sys

import pytest

from delfin.agent import sandbox as SB


@pytest.fixture(autouse=True)
def _session_is_this_interpreter(monkeypatch):
    # The install-target rule runs before the list; with the running
    # interpreter as the session venv it lets installs into it through,
    # so the cases below test the list, not that rule.
    monkeypatch.setattr(SB, "_session_python", lambda: sys.executable)


PY = sys.executable

REFUSED = [
    "pip install requests",
    "pip uninstall -y numpy",
    "pip download requests",
    "pip3 install requests",
    "pip3.11 install requests",
    f"{PY} -m pip install requests",
    f"{PY} -m pip uninstall -y numpy",
    f"{PY} -I -m pip install requests",
    f"{PY} -mpip install requests",
    f"{PY} -m pip._internal install requests",
    "python3 -m pip uninstall -y numpy",
    "env pip install requests",
    "env FOO=1 pip install requests",
    "env -i PATH=/usr/bin pip install requests",
    "env -S 'pip install requests'",
    "nohup pip install requests",
    "echo requests | xargs pip install",
    "xargs -n 1 pip install < reqs.txt",
    "xargs -I{} pip install {} < reqs.txt",
    "pip config set global.index-url http://mirror.invalid/simple",
    "pip config unset global.index-url",
    "pip config edit",
    "git config core.pager 'sh -c id'",
    "git config core.hooksPath /tmp/hooks",
    "git config --add core.fsmonitor ./x",
    "git config --unset user.name",
    "git config set core.pager less",
    "git config --file other.cfg a.b c",
    "git config -e",
    "git remote set-url origin https://example.invalid/x.git",
    "git remote add mirror https://example.invalid/x.git",
    "git remote -v add mirror https://example.invalid/x.git",
    "env git push origin main",
]

ALLOWED = [
    "pip list",
    "pip show numpy",
    "pip freeze",
    "pip check",
    "pip config list",
    "pip config get global.index-url",
    "pip3 list",
    f"{PY} -m pip list",
    "python -m pip show numpy",
    "python -m pytest -q tests/",
    "python script.py",
    "env",
    "env FOO=1 python script.py",
    "nohup python run.py",
    "ls | xargs grep foo",
    "ls | xargs -n 1 wc -l",
    "ls | xargs",
    "git status",
    "git config user.name",
    "git config --get user.name",
    "git config --list",
    "git config --get-regexp '^remote'",
    "git config get user.name",
    "git remote",
    "git remote -v",
    "git remote get-url origin",
    "git remote show origin",
]


@pytest.mark.parametrize("cmd", REFUSED)
def test_a_listed_tool_runs_only_its_listed_subcommands(cmd):
    res = SB.is_allowed(cmd)
    assert not res.allowed, (cmd, res.reason)


@pytest.mark.parametrize("cmd", ALLOWED)
def test_the_read_forms_stay_allowed(cmd):
    res = SB.is_allowed(cmd)
    assert res.allowed, (cmd, res.reason)


def test_no_tool_is_in_both_lists():
    assert not set(SB._ALLOW_FIRST_TOKEN) & set(SB._ALLOW_SUBCOMMAND)


def test_a_wrapper_does_not_hide_a_tree_walk():
    assert not SB.is_allowed("env find ~ -name x").allowed
